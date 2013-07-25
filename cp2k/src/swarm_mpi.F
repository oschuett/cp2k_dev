!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE swarm_mpi
  USE cp_files,                        ONLY: open_file
  USE cp_para_env,                     ONLY: cp_para_env_create
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_path_length
  USE machine,                         ONLY: default_output_unit
  USE message_passing,                 ONLY: &
       mp_abort, mp_any_source, mp_bcast, mp_comm_free, mp_comm_split, &
       mp_comm_split_direct, mp_environ, mp_recv, mp_send, mp_sum, mp_sync

  USE swarm_message,                  ONLY:  swarm_message_type,&
                                             swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_mpi_send,&
                                             swarm_message_mpi_recv,&
                                             swarm_message_mpi_bcast
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_mpi'

 PUBLIC :: swarm_mpi_type, swarm_mpi_init, swarm_mpi_finalize
 PUBLIC :: swarm_mpi_send_report, swarm_mpi_recv_report
 PUBLIC :: swarm_mpi_send_command, swarm_mpi_recv_command


 TYPE swarm_mpi_type
    TYPE(cp_para_env_type), POINTER          :: world => Null()
    TYPE(cp_para_env_type), POINTER          :: worker  => Null()
    TYPE(cp_para_env_type), POINTER          :: master  => Null()
    INTEGER, DIMENSION(:), ALLOCATABLE       :: wid2group
 END TYPE swarm_mpi_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_mpi_init(world_para_env, swarm_mpi, n_workers, worker_id, iw, error)
    TYPE(cp_para_env_type), POINTER          :: world_para_env
    TYPE(swarm_mpi_type)                    :: swarm_mpi
    INTEGER, INTENT(IN)                      :: n_workers
    INTEGER, INTENT(OUT)                     :: worker_id
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: n_groups_created, &
                                                pe_per_worker, subgroup, &
                                                subgroup_rank, subgroup_size, &
                                                worker_group
    LOGICAL                                  :: im_the_master
    INTEGER, DIMENSION(:), POINTER           :: group_distribution_p
    INTEGER, &
      DIMENSION(0:world_para_env%num_pe-2), &
      TARGET                                 :: group_distribution

! ====== Setup of MPI-Groups ======

    worker_id = -1
    swarm_mpi%world => world_para_env

    IF (MOD(swarm_mpi%world%num_pe-1, n_workers) /= 0) &
       STOP "number of processors-1 is not divisible by n_workers."
    IF (swarm_mpi%world%num_pe < n_workers + 1) &
       STOP "There are not enough processes for n_workers + 1. Aborting."

    pe_per_worker = (swarm_mpi%world%num_pe-1)/n_workers

    IF(iw>0) THEN
       WRITE(iw,'(A,45X,I8)') " SWARM| Number of mpi ranks", swarm_mpi%world%num_pe
       WRITE(iw,'(A,47X,I8)') " SWARM| Number of workers", n_workers
    ENDIF

    ! the last task becomes the master. Preseves node-alignment of other tasks.
    im_the_master = (swarm_mpi%world%mepos == swarm_mpi%world%num_pe-1)

    ! First split split para_env%group into a master- and a workers-groups...
    IF (im_the_master) THEN
       CALL mp_comm_split_direct(swarm_mpi%world%group, subgroup, 1)
       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
       IF(subgroup_size/=1) STOP "swarm: mp_comm_split_direct failed (master)"
       CALL cp_para_env_create(swarm_mpi%master, group=subgroup, error=error)
       !WRITE (*,*) "this is a master ", swarm_mpi%master%mepos, swarm_mpi%master%num_pe
    ELSE
       CALL mp_comm_split_direct(swarm_mpi%world%group, subgroup, 2)
       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
       !WRITE (*,*) "Hello, this is a Worker - there are ",subgroup_size, " of us."
       IF(subgroup_size/=swarm_mpi%world%num_pe-1) STOP "swarm: mp_comm_split_direct failed (worker)"
    ENDIF

    ALLOCATE(swarm_mpi%wid2group(n_workers))
    swarm_mpi%wid2group = 0

    IF(.NOT. im_the_master) THEN
       ! ...then split workers-group into n_workers groups - one for each worker.
       group_distribution_p => group_distribution
       CALL mp_comm_split(subgroup, worker_group, n_groups_created, group_distribution_p, n_subgroups=n_workers)
       worker_id = group_distribution(subgroup_rank) + 1 ! shall start by 1
       IF(n_groups_created/=n_workers) STOP "swarm: mp_comm_split failed."
       CALL cp_para_env_create(swarm_mpi%worker, group=worker_group, error=error)

       !WRITE (*,*) "this is worker ", worker_id, swarm_mpi%worker%mepos, swarm_mpi%worker%num_pe

       ! collect world-ranks of each worker groups rank-0 node
       IF(swarm_mpi%worker%mepos == 0) &
          swarm_mpi%wid2group(worker_id) = swarm_mpi%world%mepos

    ENDIF

    CALL mp_sum(swarm_mpi%wid2group, swarm_mpi%world%group)
    !WRITE (*,*), "wid2group table: ",swarm_mpi%wid2group


    CALL logger_init(swarm_mpi, error)

  END SUBROUTINE swarm_mpi_init

 ! *****************************************************************************
 ! *****************************************************************************
  SUBROUTINE logger_init(swarm_mpi, error)
    TYPE(swarm_mpi_type)                    :: swarm_mpi
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=default_path_length)       :: output_path
    INTEGER                                  :: output_unit
    TYPE(cp_logger_type), POINTER            :: logger, new_logger

! broadcast output_path to all ranks

    logger => cp_error_get_logger(error)
    output_unit = logger%default_local_unit_nr
    output_path = ""
    IF(swarm_mpi%world%mepos == swarm_mpi%world%source) THEN
       output_path = "__STD_OUT__"
       IF(output_unit /= default_output_unit)&
          INQUIRE(unit=output_unit, name=output_path)
       CALL FLUSH(output_unit) !TODO: really needed?
    END IF

    CALL mp_bcast(output_path, swarm_mpi%world%source, swarm_mpi%world%group)

    ! restore default logger on master's rank-0
    IF(.NOT. ASSOCIATED(swarm_mpi%master)) RETURN

    output_unit = -1
    IF(swarm_mpi%master%source == swarm_mpi%master%mepos) THEN
      ! open output_unit according to output_path
       output_unit = default_output_unit
       IF (output_path /= "__STD_OUT__") &
          CALL open_file(file_name=output_path,file_status="UNKNOWN",&
             file_action="WRITE", file_position="APPEND", unit_number=output_unit)
    ENDIF

    ! create default logger from output_unit
    NULLIFY(new_logger)
    CALL cp_logger_create(new_logger, para_env=swarm_mpi%master,&
       default_global_unit_nr=output_unit, close_global_unit_on_dealloc=.FALSE.,&
       template_logger=logger)
    CALL cp_add_default_logger(new_logger)
    CALL cp_logger_release(new_logger)
    CALL cp_error_set(error, logger=new_logger)
    !TODO also change the loggers of the workers?
  END SUBROUTINE logger_init

 ! *****************************************************************************
 ! *****************************************************************************
  SUBROUTINE logger_finalize(swarm_mpi, error)
    TYPE(swarm_mpi_type)                    :: swarm_mpi
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: output_unit
    TYPE(cp_logger_type), POINTER            :: logger, old_logger

    IF(.NOT.ASSOCIATED(swarm_mpi%master)) RETURN
    logger => cp_error_get_logger(error)
    output_unit = logger%default_local_unit_nr
    IF(output_unit > 0) CALL FLUSH(output_unit) !TODO: really needed?

    CALL cp_rm_default_logger() !pops the top-most logger
    old_logger => cp_error_get_logger(error)
    CALL cp_error_set(error, logger=old_logger)
  END SUBROUTINE logger_finalize


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_mpi_finalize(swarm_mpi, error)
    TYPE(swarm_mpi_type)                    :: swarm_mpi
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL mp_sync(swarm_mpi%world%group)
    CALL logger_finalize(swarm_mpi, error)

    IF(ASSOCIATED(swarm_mpi%worker)) CALL mp_comm_free(swarm_mpi%worker%group)
    IF(ASSOCIATED(swarm_mpi%master)) CALL mp_comm_free(swarm_mpi%master%group)
    NULLIFY(swarm_mpi%worker, swarm_mpi%master)
    DEALLOCATE(swarm_mpi%wid2group)
  END SUBROUTINE swarm_mpi_finalize


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_mpi_send_report(swarm_mpi, report)
    TYPE(swarm_mpi_type)                    :: swarm_mpi
    TYPE(swarm_message_type)                 :: report

    INTEGER                                  :: dest, tag

! Only rank-0 of worker group sends it's report

     IF(swarm_mpi%worker%source /= swarm_mpi%worker%mepos) RETURN

     dest = swarm_mpi%world%num_pe-1
     tag = 42
     CALL swarm_message_mpi_send(report, group=swarm_mpi%world%group, dest=dest, tag=tag)

     !!TODO: this is inefficient - could pack everything into one buffer
     !CALL mp_send(report%worker_id,       dest, tag, swarm_mpi%world%group)
     !CALL mp_send(report%iframe,          dest, tag, swarm_mpi%world%group)
     !CALL mp_send(report%Epot,            dest, tag, swarm_mpi%world%group)
     !CALL mp_send(SIZE(report%positions), dest, tag, swarm_mpi%world%group)
     !CALL mp_send(report%positions,       dest, tag, swarm_mpi%world%group)
  END SUBROUTINE swarm_mpi_send_report

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_mpi_recv_report(swarm_mpi, report)
     TYPE(swarm_mpi_type)                    :: swarm_mpi
     TYPE(swarm_message_type), INTENT(OUT)    :: report

     INTEGER                                  :: s, src, tag

     tag=42
     src=mp_any_source

     CALL swarm_message_mpi_recv(report, group=swarm_mpi%world%group, src=src, tag=tag)

     !CALL mp_recv(report%worker_id, src, tag, swarm_mpi%world%group)
     !
     !IF(src /= swarm_mpi%wid2group(report%worker_id)) &
     !   CALL mp_abort("wid2group table corrupted")
     !
     !CALL mp_recv(report%iframe,    src, tag, swarm_mpi%world%group)
     !CALL mp_recv(report%Epot,      src, tag, swarm_mpi%world%group)
     !CALL mp_recv(s,                src, tag, swarm_mpi%world%group)
     !ALLOCATE(report%positions(s))
     !CALL mp_recv(report%positions, src, tag, swarm_mpi%world%group)
  END SUBROUTINE swarm_mpi_recv_report

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_mpi_send_command(swarm_mpi, cmd)
    TYPE(swarm_mpi_type)                    :: swarm_mpi
    TYPE(swarm_message_type)                :: cmd
    INTEGER                                 :: worker_id, dest, tag

     CALL swarm_message_get(cmd, "worker_id", worker_id)
     tag = 42
     dest = swarm_mpi%wid2group(worker_id)

     CALL swarm_message_mpi_send(cmd, group=swarm_mpi%world%group, dest=dest, tag=tag)

     !!TODO: this is inefficient - could pack everything into one buffer
     !CALL mp_send(cmd%cmd_id,          dest, tag, swarm_mpi%world%group)
     !CALL mp_send(cmd%iframe,          dest, tag, swarm_mpi%world%group)
     !CALL mp_send(cmd%temperature,     dest, tag, swarm_mpi%world%group)
     !IF(ALLOCATED(cmd%positions)) THEN
     !  CALL mp_send(SIZE(cmd%positions), dest, tag, swarm_mpi%world%group)
     !  CALL mp_send(cmd%positions,       dest, tag, swarm_mpi%world%group)
     !ELSE
     !  CALL mp_send(0,                   dest, tag, swarm_mpi%world%group)
     !ENDIF
  END SUBROUTINE swarm_mpi_send_command

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_mpi_recv_command(swarm_mpi, cmd)
    TYPE(swarm_mpi_type)                    :: swarm_mpi
    TYPE(swarm_message_type), INTENT(OUT)   :: cmd

    INTEGER                                  :: s, src, tag


! This is a two step communication schema.
! First: The rank-0 of the worker groups receives the command from the master.

     
     IF(swarm_mpi%worker%source == swarm_mpi%worker%mepos) THEN
        src = swarm_mpi%world%num_pe-1!
        tag = 42
        CALL swarm_message_mpi_recv(cmd, group=swarm_mpi%world%group, src=src, tag=tag)


!        CALL mp_recv(cmd%cmd_id,      src, tag, swarm_mpi%world%group)
!        CALL mp_recv(cmd%iframe,      src, tag, swarm_mpi%world%group)
!        CALL mp_recv(cmd%temperature, src, tag, swarm_mpi%world%group)
!        CALL mp_recv(s,               src, tag, swarm_mpi%world%group)
!        IF(s > 0) THEN
!          ALLOCATE(cmd%positions(s))
!          CALL mp_recv(cmd%positions,   src, tag, swarm_mpi%world%group)
!        END IF
     ENDIF
!

      CALL swarm_message_mpi_bcast(cmd, src=swarm_mpi%worker%source, group=swarm_mpi%worker%group)

!     ! Second: The command is broadcasted within the worker group.
!     CALL mp_bcast(cmd%cmd_id,      swarm_mpi%worker%source, swarm_mpi%worker%group)
!     CALL mp_bcast(cmd%iframe,      swarm_mpi%worker%source, swarm_mpi%worker%group)
!     CALL mp_bcast(cmd%temperature, swarm_mpi%worker%source, swarm_mpi%worker%group)
!     CALL mp_bcast(s,               swarm_mpi%worker%source, swarm_mpi%worker%group)
!     IF(s > 0) THEN
!       IF(swarm_mpi%worker%source/=swarm_mpi%worker%mepos) ALLOCATE(cmd%positions(s))
!       CALL mp_bcast(cmd%positions,   swarm_mpi%worker%source, swarm_mpi%worker%group)
!     ENDIF
  END SUBROUTINE swarm_mpi_recv_command


END MODULE swarm_mpi

