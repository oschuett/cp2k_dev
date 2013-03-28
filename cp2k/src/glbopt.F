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
  USE glbopt_types,                    ONLY: glbopt_config_type,&
                                             glbopt_command_type,&
                                             glbopt_report_type,&
                                             glbopt_walker_type
  USE glbopt_walker,                   ONLY: glbopt_walker_init,&
                                             glbopt_walker_finalize,&
                                             glbopt_walker_execute_command,&
                                             glbopt_walker_assemble_report
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


#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt'

 PUBLIC :: run_global_opt, glbopt_config_type
 
  CONTAINS 
 
 

! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE run_global_opt(root_section, para_env, input_path, output_unit, error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(glbopt_config_type)                 :: glbopt_cfg

    TYPE(section_vals_type), POINTER         :: colvar_section

    CALL section_vals_val_get(root_section,"GLOBAL_OPT%NUMBER_OF_WALKERS",&
       i_val=glbopt_cfg%n_walkers,error=error)
    CALL section_vals_val_get(root_section,"GLOBAL_OPT%EMIN",&
       r_val=glbopt_cfg%Emin,error=error)
    !TODO: read from input or make educated guess
    ! TODO: more hierarchic naming scheme eg. bloom_* md_* ...
    glbopt_cfg%n_hashes = 10
    glbopt_cfg%bloom_size = 2**10
    glbopt_cfg%epot_ndown = 2
    glbopt_cfg%epot_nup = 2
    glbopt_cfg%md_mincount_max = 3
    glbopt_cfg%hop_count_max = 100
 
    glbopt_cfg%md_max_dist   = 3.0/angstrom
    glbopt_cfg%gopt_max_dist = 3.0/angstrom
    glbopt_cfg%colvar_precision = 1.0e-6

    ! make sure there are colvars
    colvar_section => section_vals_get_subs_vals(root_section,"FORCE_EVAL%SUBSYS%COLVAR",error=error)
    CALL section_vals_get(colvar_section,n_repetition=glbopt_cfg%n_colvars,error=error)
    IF(glbopt_cfg%n_colvars == 0) STOP "No colvars found. Aborting."

    CALL glbopt_driver(glbopt_cfg, root_section, input_path, output_unit, para_env, error)

  END SUBROUTINE run_global_opt

! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE glbopt_driver(glbopt_cfg, root_section, input_path, output_unit, para_env, error)
    TYPE(glbopt_config_type)                 :: glbopt_cfg
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER :: f_env_id, handle
    CHARACTER(len=*), PARAMETER :: routineN = 'glbopt_driver', &
      routineP = moduleN//':'//routineN
    CALL timeset(routineN,handle)

    IF(output_unit > 0) WRITE(output_unit,"(A)") " GLOBAL_OPT| Reporting for duty :-)"

    IF(glbopt_cfg%n_walkers == 1) THEN
       IF(output_unit > 0) WRITE(output_unit,"(A)") " GLOBAL_OPT| Running in single walker mode."
       CALL glbopt_serial_driver(glbopt_cfg, root_section, input_path, output_unit, para_env, error)
    ELSE
       IF(output_unit > 0) WRITE(output_unit,"(A)") " GLOBAL_OPT| Running in master / walkers mode."
       STOP "not implemented yet"
       !CALL driver_master_walkers(glbopt_cfg, root_section, para_env, input_path, output_unit, error)
    ENDIF

    CALL timestop(handle)
   END SUBROUTINE glbopt_driver

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_serial_driver(glbopt_cfg, root_section, input_path, output_unit, para_env, error)
       TYPE(glbopt_config_type)                 :: glbopt_cfg
       TYPE(section_vals_type), POINTER         :: root_section
       CHARACTER(LEN=*), INTENT(IN)             :: input_path
       INTEGER, INTENT(IN)                      :: output_unit
       TYPE(cp_para_env_type), POINTER          :: para_env
       TYPE(cp_error_type), INTENT(inout)       :: error

       TYPE(glbopt_master_type)                 :: master
       TYPE(glbopt_walker_type)                 :: walker
       TYPE(glbopt_command_type)                :: cmd
       TYPE(glbopt_report_type)                 :: report
       INTEGER                                  :: i
       LOGICAL                                  :: should_stop

       CALL glbopt_master_init(master, glbopt_cfg, root_section, input_path, output_unit, para_env, error)
       CALL glbopt_walker_init(walker, glbopt_cfg, root_section, input_path, output_unit, para_env,&
                                  walker_group=para_env%group, walker_id=1, error=error)

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
   SUBROUTINE glbopt_parallel_driver(glbopt_cfg, root_section, input_path, output_unit, para_env, error)
       TYPE(glbopt_config_type)                 :: glbopt_cfg
       TYPE(section_vals_type), POINTER         :: root_section
       CHARACTER(LEN=*), INTENT(IN)             :: input_path
       INTEGER, INTENT(IN)                      :: output_unit
       TYPE(cp_para_env_type), POINTER          :: para_env
       TYPE(cp_error_type), INTENT(inout)       :: error

       TYPE(cp_para_env_type), POINTER          :: walker_para_env
       TYPE(glbopt_master_type)                 :: master
       TYPE(glbopt_walker_type)                 :: walker
       TYPE(glbopt_command_type)                :: cmd
       TYPE(glbopt_report_type)                 :: report
       INTEGER                                  :: i_shutdowns
       LOGICAL                                  :: should_stop

       TYPE(glbopt_mpi_type)                    :: glbopt_mpi

       INTEGER :: n_walkers
       n_walkers = glbopt_cfg%n_walkers

       CALL glbopt_mpi_init(para_env, glbopt_mpi, n_walkers, output_unit, error)

       IF(ASSOCIATED(glbopt_mpi%walker)) THEN ! I'm a walker
          CALL glbopt_walker_init(walker, glbopt_cfg, root_section, input_path, output_unit, para_env,&
                                  walker_group=glbopt_mpi%walker%group, walker_id=1, error=error)
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
          CALL glbopt_master_init(master, glbopt_cfg, root_section, input_path, output_unit, para_env, error)
          i_shutdowns = 0
          DO WHILE(i_shutdowns < n_walkers)
             CALL glbopt_mpi_recv_report(glbopt_mpi, report)
             CALL glbopt_master_steer(master, report, cmd)
             CALL glbopt_mpi_send_command(glbopt_mpi, cmd)
             IF(cmd%id == GLBOPT_CMD_SHUTDOWN) i_shutdowns = i_shutdowns + 1
          END DO
          CALL glbopt_master_finalize(master)
       END IF

       CALL glbopt_mpi_finalize(glbopt_mpi, output_unit, error)

   END SUBROUTINE glbopt_parallel_driver


END MODULE glbopt

