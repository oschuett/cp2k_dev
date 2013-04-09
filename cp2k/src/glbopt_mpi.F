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
  USE glbopt_types,                    ONLY: glbopt_config_type,&
                                             glbopt_command_type,&
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

    INTEGER :: n_groups_created, &
      pe_per_walker, subgroup, subgroup_rank, subgroup_size, walker_group
    LOGICAL                                  :: im_the_master
    INTEGER, DIMENSION(:), POINTER           :: group_distribution_p
    INTEGER, DIMENSION(0:world_para_env%num_pe-2), TARGET :: group_distribution

    walker_id = -1
    glbopt_mpi%world => world_para_env

    ! ====== Setup of MPI-Groups ======
    IF (MOD(glbopt_mpi%world%num_pe-1, n_walkers) /= 0) &
       STOP "number of processors-1 is not divisible by n_walkers."
    IF (glbopt_mpi%world%num_pe < n_walkers + 1) &
       STOP "There are not enough processes for n_walkers + 1. Aborting."

    pe_per_walker = (glbopt_mpi%world%num_pe-1)/n_walkers

    WRITE(*,*) "!!! Number of mpi ranks: ", glbopt_mpi%world%num_pe
    WRITE(*,*) "!!! Number of walkers: ", n_walkers


    ! the last task becomes the master. Preseves node-alignment of other tasks.
    im_the_master = (glbopt_mpi%world%mepos == glbopt_mpi%world%num_pe-1)

    ! First split split para_env%group into a master- and a walkers-groups...
    IF (im_the_master) THEN
       CALL mp_comm_split_direct(glbopt_mpi%world%group, subgroup, 1)
       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
       WRITE (*,*) "Hello, this is a Master - there are ",subgroup_size, " of us."
       IF(subgroup_size/=1) STOP "glbopt: mp_comm_split_direct failed (master)"
    ELSE
       CALL mp_comm_split_direct(glbopt_mpi%world%group, subgroup, 2)
       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
       WRITE (*,*) "Hello, this is a Walker - there are ",subgroup_size, " of us."
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
    WRITE (*,*), "wid2group table: ",glbopt_mpi%wid2group

  END SUBROUTINE glbopt_mpi_init



! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_mpi_finalize(glbopt_mpi, error)
    TYPE(glbopt_mpi_type)                    :: glbopt_mpi
    TYPE(cp_error_type), INTENT(inout)       :: error

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
     CALL mp_send(cmd%id,              dest, tag, glbopt_mpi%world%group)
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
     tag = 42
     src = glbopt_mpi%world%num_pe-1

     CALL mp_recv(cmd%id,          src, tag, glbopt_mpi%world%group)
     CALL mp_recv(cmd%temperature, src, tag, glbopt_mpi%world%group)
     CALL mp_recv(s,               src, tag, glbopt_mpi%world%group)
     ! TODO atom numbers do not change, make this more efficient
     IF(ALLOCATED(cmd%positions)) DEALLOCATE(cmd%positions)
     ALLOCATE(cmd%positions(s))
     CALL mp_recv(cmd%positions,   src, tag, glbopt_mpi%world%group)
  END SUBROUTINE glbopt_mpi_recv_command

  
!! *****************************************************************************
!! *****************************************************************************
!   SUBROUTINE driver_master_walkers(glbopt_cfg, root_section, para_env, input_path, output_unit, error)
!    TYPE(glbopt_config_type)                 :: glbopt_cfg
!    TYPE(section_vals_type), POINTER         :: root_section
!    TYPE(cp_para_env_type), POINTER          :: para_env
!    CHARACTER(LEN=*), INTENT(IN)             :: input_path
!    INTEGER, INTENT(IN)                      :: output_unit
!    TYPE(cp_error_type), INTENT(inout)       :: error
!
!    INTEGER :: f_env_id, my_walker_id, n_groups_created, &
!      pe_per_walker, subgroup, subgroup_rank, subgroup_size, walker_group
!    LOGICAL                                  :: im_the_master
!    INTEGER, DIMENSION(:), POINTER           :: group_distribution_p
!    INTEGER, &
!      DIMENSION(0:para_env%num_pe-2), TARGET :: group_distribution
!
!
!    ! ====== Setup of MPI-Groups ======
!    IF (MOD(para_env%num_pe-1, glbopt_cfg%n_walkers) /= 0) &
!       STOP "number of processors-1 is not divisible by n_walkers."
!    IF (para_env%num_pe < glbopt_cfg%n_walkers + 1) &
!       STOP "There are not enough processes for n_walkers + 1. Aborting."
!
!    pe_per_walker = (para_env%num_pe-1)/glbopt_cfg%n_walkers
!
!    WRITE(*,*) "!!! Number of mpi ranks: ",para_env%num_pe
!    WRITE(*,*) "!!! Number of walkers: ", glbopt_cfg%n_walkers
!    WRITE(*,*) "!!! Number of colvars: ", glbopt_cfg%n_colvars
!
!    ! the last task becomes the master. Keeps node-alignment of other tasks.
!    im_the_master = (para_env%mepos == para_env%num_pe-1)
!
!    ! First split split para_env%group into a master- and a walkers-groups...
!    IF (im_the_master) THEN
!       CALL mp_comm_split_direct(para_env%group, subgroup, 1)
!       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
!       !WRITE (*,*) "Hello, this is a Master - there are ",subgroup_size, " of us."
!       IF(subgroup_size/=1) STOP "glbopt: mp_comm_split_direct failed (master)"
!    ELSE
!       CALL mp_comm_split_direct(para_env%group, subgroup, 2)
!       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
!       !WRITE (*,*) "Hello, this is a Walker - there are ",subgroup_size, " of us."
!       IF(subgroup_size/=para_env%num_pe-1) STOP "glbopt: mp_comm_split_direct failed (walker)"
!    ENDIF
!
!
!    IF (im_the_master) THEN
!       CALL glbopt_drive_master(glbopt_cfg, root_section, para_env, error)
!    ELSE
!       ! ...then split walkers-group into n_walkers groups - one for each walker.
!       group_distribution_p => group_distribution
!       CALL mp_comm_split(subgroup, walker_group, n_groups_created, group_distribution_p, n_subgroups=glbopt_cfg%n_walkers)
!       my_walker_id = group_distribution(subgroup_rank)
!       IF(n_groups_created/=glbopt_cfg%n_walkers) STOP "glbopt: mp_comm_split failed."
!       ! f_env_id = walker_create_f_env(root_section, input_path, walker_group, my_walker_id, error)
!       CALL glbopt_drive_walker(glbopt_cfg, root_section, input_path, output_unit, para_env, walker_group, my_walker_id, error)
!       CALL mp_comm_free(walker_group)
!    ENDIF
!
!    CALL mp_comm_free(subgroup)
!
!  END SUBROUTINE driver_master_walkers
!
!
!! *****************************************************************************
!  FUNCTION walker_create_f_env(root_section, input_path, mpi_comm, walker_id, error) RESULT(f_env_id)
!    TYPE(section_vals_type), POINTER         :: root_section
!    CHARACTER(LEN=*), INTENT(IN)             :: input_path
!    INTEGER, INTENT(IN)                      :: mpi_comm, walker_id
!    TYPE(cp_error_type), INTENT(inout)       :: error
!    INTEGER                                  :: ierr, f_env_id
!
!    CHARACTER(len=default_path_length)       :: output_path
!    CHARACTER(len=default_string_length)     :: new_project_name, &
!                                                project_name, walker_name
!
!! change the PROJECT_NAME to avoid clashes during backup of RESTART-files
!    IF(walker_id > 999) STOP "walker_create_f_env: Did not expect so many walkers."
!    CALL section_vals_val_get(root_section, "GLOBAL%PROJECT_NAME", c_val=project_name,error=error)
!    WRITE(walker_name,"('WALKER',I3.3)"), walker_id
!    !walker_name = "walker"//ADJUSTL(cp_to_string(walker_id))
!    new_project_name = TRIM(project_name)//"-"//walker_name
!    CALL section_vals_val_set(root_section, "GLOBAL%PROJECT_NAME", c_val=new_project_name,error=error)
!
!    output_path = TRIM(new_project_name)//".out"
!    WRITE (*,*) output_path
!
!    CALL create_force_env(f_env_id, &
!                          input_path=input_path, &
!                          input=root_section, &
!                          output_path=output_path, &
!                          mpi_comm=mpi_comm,&
!                          ierr=ierr)
!
!    IF(ierr/=0) STOP "glbopt: failed to create force_env"
!
!  END FUNCTION walker_create_f_env
!
!
END MODULE glbopt_mpi

