!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE glbopt_input
  USE glbopt_master,                   ONLY: glbopt_master_finalize,&
                                             glbopt_master_init,&
                                             glbopt_master_steer,&
                                             glbopt_master_type
  USE glbopt_mpi,                      ONLY: glbopt_mpi_finalize,&
                                             glbopt_mpi_init,&
                                             glbopt_mpi_recv_command,&
                                             glbopt_mpi_recv_report,&
                                             glbopt_mpi_send_command,&
                                             glbopt_mpi_send_report,&
                                             glbopt_mpi_type
  USE glbopt_types,                    ONLY: GLBOPT_CMD_SHUTDOWN,&
                                             glbopt_command_type,&
                                             glbopt_report_type
  USE glbopt_walker,                   ONLY: glbopt_walker_assemble_report,&
                                             glbopt_walker_execute_command,&
                                             glbopt_walker_finalize,&
                                             glbopt_walker_init,&
                                             glbopt_walker_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt'

 PUBLIC :: run_global_opt
 
  CONTAINS 
 
 

! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE run_global_opt(root_section, para_env, globenv, input_path, error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'run_global_opt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, n_walkers
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,root_section,&
          "GLOBAL_OPT%PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    IF(iw > 0) WRITE(iw,"(A)") " GLOBAL_OPT| Ready to roll :-)"

    CALL section_vals_val_get(root_section,"GLOBAL_OPT%NUMBER_OF_WALKERS",&
       i_val=n_walkers,error=error)

    IF(n_walkers == 1) THEN
       IF(iw > 0) WRITE(iw,"(A)") " GLOBAL_OPT| Running in single walker mode."
       CALL glbopt_serial_driver(root_section, input_path, para_env, globenv, error)
    ELSE
       IF(iw > 0) WRITE(iw,"(A)") " GLOBAL_OPT| Running in master / walkers mode."
       !printkey iw passed on for output from glbopt_mpi_init()
       CALL glbopt_parallel_driver(n_walkers, root_section, input_path, para_env, globenv, iw, error)
    ENDIF

    CALL timestop(handle)
   END SUBROUTINE run_global_opt

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_serial_driver(root_section, input_path, para_env, globenv, error)
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    LOGICAL                                  :: should_stop
    TYPE(glbopt_command_type)                :: cmd
    TYPE(glbopt_master_type)                 :: master
    TYPE(glbopt_report_type)                 :: report
    TYPE(glbopt_walker_type)                 :: walker

       CALL glbopt_master_init(master, para_env, globenv, root_section, n_walkers=1, error=error)
       CALL glbopt_walker_init(walker, para_env, root_section, input_path, &
                                   walker_id=1, error=error)

       should_stop = .FALSE.
       DO WHILE(.NOT. should_stop)
          CALL glbopt_walker_assemble_report(walker, report)
          CALL glbopt_master_steer(master, report, cmd)
          CALL glbopt_walker_execute_command(walker, cmd, should_stop)
       END DO

       CALL glbopt_walker_finalize(walker)
       CALL glbopt_master_finalize(master)

    END SUBROUTINE glbopt_serial_driver


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_parallel_driver(n_walkers, root_section, input_path, para_env, globenv, iw, error)
    INTEGER, INTENT(IN)                      :: n_walkers
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i_shutdowns, walker_id
    LOGICAL                                  :: should_stop
    TYPE(glbopt_command_type)                :: cmd
    TYPE(glbopt_master_type)                 :: master
    TYPE(glbopt_mpi_type)                    :: glbopt_mpi
    TYPE(glbopt_report_type)                 :: report
    TYPE(glbopt_walker_type)                 :: walker

       CALL glbopt_mpi_init(para_env, glbopt_mpi, n_walkers, walker_id, iw, error=error)

       IF(ASSOCIATED(glbopt_mpi%walker)) THEN ! I'm a walker
          CALL glbopt_walker_init(walker, glbopt_mpi%walker,&
                     root_section, input_path, walker_id=walker_id, error=error)

          should_stop = .FALSE.
          DO WHILE(.NOT. should_stop)
             CALL glbopt_walker_assemble_report(walker, report)
             CALL glbopt_mpi_send_report(glbopt_mpi, report)
             CALL glbopt_mpi_recv_command(glbopt_mpi, cmd)
             CALL glbopt_walker_execute_command(walker, cmd, should_stop)
          END DO
          CALL glbopt_walker_finalize(walker)

       !------------------------------------------------------------------------
       ELSE  ! I'm the master
          CALL glbopt_master_init(master, glbopt_mpi%master, globenv, root_section,&
                 n_walkers, error)
          i_shutdowns = 0
          DO WHILE(i_shutdowns < n_walkers)
             CALL glbopt_mpi_recv_report(glbopt_mpi, report)
             CALL glbopt_master_steer(master, report, cmd)
             CALL glbopt_mpi_send_command(glbopt_mpi, report%walker_id, cmd)
             IF(cmd%cmd_id == GLBOPT_CMD_SHUTDOWN) i_shutdowns = i_shutdowns + 1
          END DO
          CALL glbopt_master_finalize(master)
       END IF

       CALL glbopt_mpi_finalize(glbopt_mpi, error)

   END SUBROUTINE glbopt_parallel_driver


END MODULE glbopt

