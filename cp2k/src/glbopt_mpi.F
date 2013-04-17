!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_mpi
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_para_env,                     ONLY: cp_para_env_create
  USE f77_interface,                   ONLY: create_force_env
  USE glbopt_master,                   ONLY: glbopt_master_type,&
                                             glbopt_master_init,&
                                             glbopt_master_steer,&
                                             glbopt_master_finalize
  USE glbopt_types,                    ONLY: glbopt_command_type,&
                                             glbopt_report_type
  USE glbopt_walker,                   ONLY: glbopt_walker_init,&
                                             glbopt_walker_finalize,&
                                             glbopt_walker_execute_command,&
                                             glbopt_walker_assemble_report
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length
  USE message_passing,                 ONLY: mp_comm_free,&
                                             mp_comm_split,&
                                             mp_comm_split_direct,&
                                             mp_environ, mp_sum
  USE timings,                         ONLY: timeset,&
                                             timestop

  USE physcon,                         ONLY: angstrom
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: integer_t, real_t
  USE input_constants,                 ONLY: low_print_level, add_last_numeric
  USE message_passing,                 ONLY: mp_any_source,mp_any_tag,&
                                             mp_recv,&
                                             mp_send, mp_abort
  USE machine,                         ONLY: default_output_unit
  USE message_passing,                 ONLY: mp_bcast, mp_sync
  USE cp_files,                        ONLY: close_file,&
                                             open_file
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_mpi'

 PUBLIC :: glbopt_mpi_type, glbopt_mpi_init, glbopt_mpi_finalize
 PUBLIC :: glbopt_mpi_send_report, glbopt_mpi_recv_report
 PUBLIC :: glbopt_mpi_send_command, glbopt_mpi_recv_command


 TYPE glbopt_mpi_type
    TYPE(cp_para_env_type), POINTER          :: world => Null()
    TYPE(cp_para_env_type), POINTER          :: walker  => Null()
    TYPE(cp_para_env_type), POINTER          :: master  => Null()
    INTEGER, DIMENSION(:), ALLOCATABLE       :: wid2group
 END TYPE glbopt_mpi_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_mpi_init(world_para_env, glbopt_mpi, n_walkers, walker_id, iw, error)
    TYPE(cp_para_env_type), POINTER          :: world_para_env
    TYPE(glbopt_mpi_type)                    :: glbopt_mpi
    INTEGER, INTENT(IN)                      :: n_walkers
    INTEGER, INTENT(OUT)                     :: walker_id
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(LEN=default_path_length)       :: output_path
    INTEGER                                  :: output_unit
    INTEGER :: n_groups_created, &
      pe_per_walker, subgroup, subgroup_rank, subgroup_size, walker_group
    LOGICAL                                  :: im_the_master
    INTEGER, DIMENSION(:), POINTER           :: group_distribution_p
    INTEGER, DIMENSION(0:world_para_env%num_pe-2), TARGET :: group_distribution


    ! ====== Setup of MPI-Groups ======
    walker_id = -1
    glbopt_mpi%world => world_para_env

    IF (MOD(glbopt_mpi%world%num_pe-1, n_walkers) /= 0) &
       STOP "number of processors-1 is not divisible by n_walkers."
    IF (glbopt_mpi%world%num_pe < n_walkers + 1) &
       STOP "There are not enough processes for n_walkers + 1. Aborting."

    pe_per_walker = (glbopt_mpi%world%num_pe-1)/n_walkers

    IF(iw>0) THEN
       WRITE(iw,*) "GLOBAL_OPT| Number of mpi ranks: ", glbopt_mpi%world%num_pe
       WRITE(iw,*) "GLOBAL_OPT| Number of walkers: ", n_walkers
    ENDIF

    ! the last task becomes the master. Preseves node-alignment of other tasks.
    im_the_master = (glbopt_mpi%world%mepos == glbopt_mpi%world%num_pe-1)

    ! First split split para_env%group into a master- and a walkers-groups...
    IF (im_the_master) THEN
       CALL mp_comm_split_direct(glbopt_mpi%world%group, subgroup, 1)
       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
       !WRITE (*,*) "Hello, this is a Master - there are ",subgroup_size, " of us."
       IF(subgroup_size/=1) STOP "glbopt: mp_comm_split_direct failed (master)"
    ELSE
       CALL mp_comm_split_direct(glbopt_mpi%world%group, subgroup, 2)
       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
       !WRITE (*,*) "Hello, this is a Walker - there are ",subgroup_size, " of us."
       IF(subgroup_size/=glbopt_mpi%world%num_pe-1) STOP "glbopt: mp_comm_split_direct failed (walker)"
    ENDIF

    CALL cp_para_env_create(glbopt_mpi%master, group=subgroup, error=error)

    ALLOCATE(glbopt_mpi%wid2group(n_walkers))
    glbopt_mpi%wid2group = 0

    IF(.NOT. im_the_master) THEN
       ! ...then split walkers-group into n_walkers groups - one for each walker.
       group_distribution_p => group_distribution
       CALL mp_comm_split(subgroup, walker_group, n_groups_created, group_distribution_p, n_subgroups=n_walkers)
       walker_id = group_distribution(subgroup_rank) + 1 ! shall start by 1
       IF(n_groups_created/=n_walkers) STOP "glbopt: mp_comm_split failed."
       CALL cp_para_env_create(glbopt_mpi%walker, group=walker_group, error=error)

       ! collect world-ranks of each walker groups rank-0 node
       IF(glbopt_mpi%walker%mepos == 0) &
          glbopt_mpi%wid2group(walker_id) = glbopt_mpi%world%mepos

    ENDIF

    CALL mp_sum(glbopt_mpi%wid2group, glbopt_mpi%world%group)
    !WRITE (*,*), "wid2group table: ",glbopt_mpi%wid2group


    CALL logger_init(glbopt_mpi, error)

  END SUBROUTINE glbopt_mpi_init

 ! *****************************************************************************
 ! *****************************************************************************
  SUBROUTINE logger_init(glbopt_mpi, error)
    TYPE(glbopt_mpi_type)                    :: glbopt_mpi
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_logger_type), POINTER            :: logger, new_logger
    CHARACTER(LEN=default_path_length)       :: output_path
    INTEGER                                  :: output_unit

    ! broadcast output_path to all ranks
    logger => cp_error_get_logger(error)
    output_unit = logger%default_local_unit_nr
    output_path = ""
    IF(glbopt_mpi%world%mepos == glbopt_mpi%world%source) THEN
       output_path = "__STD_OUT__"
       IF(output_unit /= default_output_unit)&
          INQUIRE(unit=output_unit, name=output_path)
       CALL FLUSH(output_unit) !TODO: really needed?
    END IF

    CALL mp_bcast(output_path, glbopt_mpi%world%source, glbopt_mpi%world%group)

    ! restore default logger on master's rank-0
    IF(.NOT. ASSOCIATED(glbopt_mpi%master)) RETURN

    output_unit = -1
    IF(glbopt_mpi%master%source == glbopt_mpi%master%mepos) THEN
      ! open output_unit according to output_path
       output_unit = default_output_unit
       IF (output_path /= "__STD_OUT__") &
          CALL open_file(file_name=output_path,file_status="UNKNOWN",&
             file_action="WRITE", file_position="APPEND", unit_number=output_unit)
    ENDIF

    ! create default logger from output_unit
    NULLIFY(new_logger)
    CALL cp_logger_create(new_logger, para_env=glbopt_mpi%master,&
       default_global_unit_nr=output_unit, close_global_unit_on_dealloc=.FALSE.)
    CALL cp_add_default_logger(new_logger)
    CALL cp_logger_release(new_logger)
    CALL cp_error_set(error, logger=new_logger)
    !TODO also change the loggers of the walkers?
  END SUBROUTINE logger_init

 ! *****************************************************************************
 ! *****************************************************************************
  SUBROUTINE logger_finalize(glbopt_mpi, error)
    TYPE(glbopt_mpi_type)                    :: glbopt_mpi
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(cp_logger_type), POINTER            :: logger, old_logger
    INTEGER                                  :: output_unit

    IF(.NOT.ASSOCIATED(glbopt_mpi%master)) RETURN
    logger => cp_error_get_logger(error)
    output_unit = logger%default_local_unit_nr
    IF(output_unit > 0) CALL FLUSH(output_unit) !TODO: really needed?

    CALL cp_rm_default_logger() !pops the top-most logger
    old_logger => cp_error_get_logger(error)
    CALL cp_error_set(error, logger=old_logger)
  END SUBROUTINE logger_finalize


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_mpi_finalize(glbopt_mpi, error)
    TYPE(glbopt_mpi_type)                    :: glbopt_mpi
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL mp_sync(glbopt_mpi%world%group)
    CALL logger_finalize(glbopt_mpi, error)

    If(ASSOCIATED(glbopt_mpi%walker)) CALL mp_comm_free(glbopt_mpi%walker%group)
    CALL mp_comm_free(glbopt_mpi%master%group)
    NULLIFY(glbopt_mpi%walker, glbopt_mpi%master)
    DEALLOCATE(glbopt_mpi%wid2group)
  END SUBROUTINE glbopt_mpi_finalize


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE glbopt_mpi_send_report(glbopt_mpi, report)
     TYPE(glbopt_mpi_type)                    :: glbopt_mpi
     TYPE(glbopt_report_type)                 :: report

     INTEGER :: tag, dest

     ! Only rank-0 of walker group sends it's report
     IF(glbopt_mpi%walker%source /= glbopt_mpi%walker%mepos) RETURN

     dest = glbopt_mpi%world%num_pe-1
     tag = 42
     !TODO: this is inefficient - could pack everything into one buffer
     CALL mp_send(report%walker_id,       dest, tag, glbopt_mpi%world%group)
     CALL mp_send(report%Epot,            dest, tag, glbopt_mpi%world%group)
     CALL mp_send(SIZE(report%positions), dest, tag, glbopt_mpi%world%group)
     CALL mp_send(report%positions,       dest, tag, glbopt_mpi%world%group)
  END SUBROUTINE glbopt_mpi_send_report

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE glbopt_mpi_recv_report(glbopt_mpi, report)
     TYPE(glbopt_mpi_type)                    :: glbopt_mpi
     TYPE(glbopt_report_type)                 :: report
     INTEGER :: src, tag, s

     tag=42
     src=mp_any_source

     CALL mp_recv(report%walker_id, src, tag, glbopt_mpi%world%group)

     IF(src /= glbopt_mpi%wid2group(report%walker_id)) &
        CALL mp_abort("wid2group table corrupted")

     CALL mp_recv(report%Epot,      src, tag, glbopt_mpi%world%group)
     CALL mp_recv(s,                src, tag, glbopt_mpi%world%group)
     ! TODO atom numbers do not change, make this more efficient
     IF(ALLOCATED(report%positions)) DEALLOCATE(report%positions)
     ALLOCATE(report%positions(s))
     CALL mp_recv(report%positions, src, tag, glbopt_mpi%world%group)
  END SUBROUTINE glbopt_mpi_recv_report

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE glbopt_mpi_send_command(glbopt_mpi, walker_id, cmd)
     TYPE(glbopt_mpi_type)                    :: glbopt_mpi
     INTEGER       :: walker_id
     TYPE(glbopt_command_type)                 :: cmd

     INTEGER :: tag, dest
     tag = 42
     dest = glbopt_mpi%wid2group(walker_id)

     !TODO: this is inefficient - could pack everything into one buffer
     CALL mp_send(cmd%cmd_id,          dest, tag, glbopt_mpi%world%group)
     CALL mp_send(cmd%temperature,     dest, tag, glbopt_mpi%world%group)
     CALL mp_send(SIZE(cmd%positions), dest, tag, glbopt_mpi%world%group)
     CALL mp_send(cmd%positions,       dest, tag, glbopt_mpi%world%group)
  END SUBROUTINE glbopt_mpi_send_command

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE glbopt_mpi_recv_command(glbopt_mpi, cmd)
     TYPE(glbopt_mpi_type)                    :: glbopt_mpi
     TYPE(glbopt_command_type)                :: cmd

     INTEGER :: tag, src, s

     IF(ALLOCATED(cmd%positions)) DEALLOCATE(cmd%positions)

     ! This is a two step communication schema.
     ! First: The rank-0 of the walker groups receives the command from the master.
     IF(glbopt_mpi%walker%source == glbopt_mpi%walker%mepos) THEN 
        tag = 42
        src = glbopt_mpi%world%num_pe-1
        CALL mp_recv(cmd%cmd_id,      src, tag, glbopt_mpi%world%group)
        CALL mp_recv(cmd%temperature, src, tag, glbopt_mpi%world%group)
        CALL mp_recv(s,               src, tag, glbopt_mpi%world%group)
        ! TODO atom numbers do not change, make this more efficient
        ALLOCATE(cmd%positions(s))
        CALL mp_recv(cmd%positions,   src, tag, glbopt_mpi%world%group)
     ENDIF

     ! Second: The command is broadcasted within the walker group.
     CALL mp_bcast(cmd%cmd_id,      glbopt_mpi%walker%source, glbopt_mpi%walker%group)
     CALL mp_bcast(cmd%temperature, glbopt_mpi%walker%source, glbopt_mpi%walker%group)
     CALL mp_bcast(s,               glbopt_mpi%walker%source, glbopt_mpi%walker%group)
     IF(glbopt_mpi%walker%source /= glbopt_mpi%walker%mepos) ALLOCATE(cmd%positions(s))
     CALL mp_bcast(cmd%positions,   glbopt_mpi%walker%source, glbopt_mpi%walker%group)
  END SUBROUTINE glbopt_mpi_recv_command


END MODULE glbopt_mpi

