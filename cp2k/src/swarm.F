!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE swarm
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE swarm_input
  USE swarm_master,                   ONLY: swarm_master_finalize,&
                                             swarm_master_init,&
                                             swarm_master_steer,&
                                             swarm_master_type
  USE swarm_mpi,                      ONLY: swarm_mpi_finalize,&
                                             swarm_mpi_init,&
                                             swarm_mpi_recv_command,&
                                             swarm_mpi_recv_report,&
                                             swarm_mpi_send_command,&
                                             swarm_mpi_send_report,&
                                             swarm_mpi_type
  USE swarm_message,                  ONLY:  swarm_message_type,&
                                             swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_mpi_send,&
                                             swarm_message_mpi_recv

  USE swarm_worker,                   ONLY:  swarm_worker_execute,&
                                             swarm_worker_finalize,&
                                             swarm_worker_init,&
                                             swarm_worker_create_empty_report,&
                                             swarm_worker_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE kinds,                           ONLY: default_string_length

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm'

 PUBLIC :: run_swarm
 
  CONTAINS 
 
 

! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE run_swarm(root_section, para_env, globenv, input_path, error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'run_SWARM', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, n_workers
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,root_section,&
          "SWARM%PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    IF(iw > 0) WRITE(iw,"(A)") " SWARM| Ready to roll :-)"

    CALL section_vals_val_get(root_section,"SWARM%NUMBER_OF_WORKERS",&
       i_val=n_workers,error=error)

    IF(n_workers == 1) THEN
       IF(iw > 0) WRITE(iw,"(A)") " SWARM| Running in single worker mode."
       CALL swarm_serial_driver(root_section, input_path, para_env, globenv, error)
    ELSE
       IF(iw > 0) WRITE(iw,"(A)") " SWARM| Running in master / workers mode."
       !printkey iw passed on for output from swarm_mpi_init()
       CALL swarm_parallel_driver(n_workers, root_section, input_path, para_env, globenv, iw, error)
    ENDIF

    CALL timestop(handle)
   END SUBROUTINE run_swarm

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_serial_driver(root_section, input_path, para_env, globenv, error)
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    LOGICAL                                  :: should_stop
    TYPE(swarm_message_type)                :: cmd, report
    TYPE(swarm_master_type)                 :: master
    TYPE(swarm_worker_type)                 :: worker

       CALL swarm_master_init(master, para_env, globenv, root_section, n_workers=1, error=error)
       CALL swarm_worker_init(worker, para_env, root_section, input_path, &
                                   worker_id=1, error=error)

       CALL swarm_worker_create_empty_report(worker, report)
       should_stop = .FALSE.
       DO WHILE(.NOT. should_stop)
          CALL swarm_master_steer(master, report, cmd)
          CALL swarm_worker_execute(worker, cmd, report, should_stop)
       END DO

       CALL swarm_worker_finalize(worker)
       CALL swarm_master_finalize(master)

    END SUBROUTINE swarm_serial_driver



! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_parallel_driver(n_workers, root_section, input_path, para_env, globenv, iw, error)
    INTEGER, INTENT(IN)                      :: n_workers
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i_shutdowns, worker_id
    LOGICAL                                  :: should_stop
    CHARACTER(len=default_string_length)     :: command
    TYPE(swarm_message_type)                :: cmd, report
    TYPE(swarm_master_type)                 :: master
    TYPE(swarm_mpi_type)                    :: swarm_mpi
    TYPE(swarm_worker_type)                 :: worker

       CALL swarm_mpi_init(para_env, swarm_mpi, n_workers, worker_id, iw, error=error)

       IF(ASSOCIATED(swarm_mpi%worker)) THEN ! I'm a worker
          CALL swarm_worker_init(worker, swarm_mpi%worker,&
                     root_section, input_path, worker_id=worker_id, error=error)

          CALL swarm_worker_create_empty_report(worker, report)
          should_stop = .FALSE.
          DO WHILE(.NOT. should_stop)
             CALL swarm_mpi_send_report(swarm_mpi, report)
             CALL swarm_mpi_recv_command(swarm_mpi, cmd)
             CALL swarm_worker_execute(worker, cmd, report, should_stop)
          END DO
          CALL swarm_worker_finalize(worker)

       !------------------------------------------------------------------------
       ELSE  ! I'm the master
          CALL swarm_master_init(master, swarm_mpi%master, globenv, root_section,&
                 n_workers, error)
          i_shutdowns = 0
          DO WHILE(i_shutdowns < n_workers)
             CALL swarm_mpi_recv_report(swarm_mpi, report)
             CALL swarm_master_steer(master, report, cmd)
             CALL swarm_message_get(report, "worker_id", worker_id)
             CALL swarm_message_get(cmd, "command", command)
             CALL swarm_mpi_send_command(swarm_mpi, worker_id, cmd)
             IF(TRIM(command) == "shutdown") i_shutdowns = i_shutdowns + 1
          END DO
          CALL swarm_master_finalize(master)
       END IF

       CALL swarm_mpi_finalize(swarm_mpi, error)

   END SUBROUTINE swarm_parallel_driver


END MODULE swarm

