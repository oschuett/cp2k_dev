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
  USE cp_para_types,                   ONLY: cp_para_env_type
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
                                             glbopt_walker_assemble_report,&
                                             glbopt_walker_type
  USE glbopt_mpi,                      ONLY: glbopt_mpi_init,&
                                             glbopt_mpi_finalize,&
                                             glbopt_mpi_type
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
                                             mp_environ
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
  USE glbopt_input
  USE glbopt_types,                    ONLY: GLBOPT_CMD_SHUTDOWN
  USE glbopt_mpi,                      ONLY: glbopt_mpi_type,&
                                             glbopt_mpi_init,&
                                             glbopt_mpi_finalize,&
                                             glbopt_mpi_send_report,&
                                             glbopt_mpi_recv_report,&
                                             glbopt_mpi_send_command,&
                                             glbopt_mpi_recv_command

  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE global_types,                    ONLY: global_environment_type

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

    TYPE(cp_logger_type), POINTER            :: logger

    INTEGER ::  handle, iw, n_walkers
    CHARACTER(len=*), PARAMETER :: routineN = 'run_global_opt', &
      routineP = moduleN//':'//routineN
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

       TYPE(glbopt_master_type)                 :: master
       TYPE(glbopt_walker_type)                 :: walker
       TYPE(glbopt_command_type)                :: cmd
       TYPE(glbopt_report_type)                 :: report
       INTEGER                                  :: i
       LOGICAL                                  :: should_stop

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

       TYPE(cp_para_env_type), POINTER          :: walker_para_env
       TYPE(glbopt_master_type)                 :: master
       TYPE(glbopt_walker_type)                 :: walker
       TYPE(glbopt_command_type)                :: cmd
       TYPE(glbopt_report_type)                 :: report
       INTEGER                                  :: i_shutdowns, walker_id
       LOGICAL                                  :: should_stop

       TYPE(glbopt_mpi_type)                    :: glbopt_mpi


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

