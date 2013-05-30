!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE swarm_worker
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type,&
                                             pack_subsys_particles
  USE f77_interface,                   ONLY: create_force_env,&
                                             destroy_force_env,&
                                             f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE geo_opt,                         ONLY: cp_geo_opt
  USE swarm_types,                    ONLY:  swarm_cmd_shutdown,&
                                             swarm_command_type,&
                                             swarm_report_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE md_run,                          ONLY: qs_mol_dyn
  USE mdctrl_types,                    ONLY: mdctrl_type
  USE parallel_rng_types,              ONLY: reset_to_next_rng_substream
  USE physcon,                         ONLY: angstrom,&
                                             femtoseconds,&
                                             kelvin
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE glbopt_worker,                   ONLY:  glbopt_worker_init,&
                                              glbopt_worker_finalize,&
                                              glbopt_worker_execute,&
                                              glbopt_worker_type

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_worker'

 PUBLIC :: swarm_worker_init, swarm_worker_finalize
 PUBLIC :: swarm_worker_execute
 PUBLIC :: swarm_worker_type
 PUBLIC :: swarm_worker_create_empty_report

  TYPE swarm_worker_type
   PRIVATE
   INTEGER                                  :: id
   INTEGER                                  :: iw
   TYPE(cp_error_type)                      :: error
   TYPE(glbopt_worker_type)                 :: glbopt
 END TYPE swarm_worker_type

 CONTAINS




! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_worker_init(worker, para_env, root_section,&
                 input_path, worker_id, error)
    TYPE(swarm_worker_type), INTENT(INOUT)   :: worker
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(in)                      :: worker_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_logger_type), POINTER            :: logger

    worker%error = error
    worker%id = worker_id
    ! getting an output unit for logging
    ! TODO: this seems to be too early, because the the force env has not been created
    logger => cp_error_get_logger(worker%error)
    worker%iw = cp_print_key_unit_nr(logger,root_section,&
          "SWARM%PRINT%WORKER_RUN_INFO",extension=".workerLog", error=worker%error)

    CALL glbopt_worker_init(worker%glbopt, para_env, root_section,&
                 input_path, worker_id, worker%iw, error)

  END SUBROUTINE swarm_worker_init



! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_worker_finalize(worker)
    TYPE(swarm_worker_type), INTENT(INOUT)  :: worker

    CALL glbopt_worker_finalize(worker%glbopt)

   END SUBROUTINE swarm_worker_finalize


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_worker_execute(worker, cmd, report, should_stop)
     TYPE(swarm_worker_type), INTENT(INOUT)  :: worker
     TYPE(swarm_command_type), INTENT(IN)    :: cmd
     TYPE(swarm_report_type), INTENT(OUT)    :: report
     LOGICAL, INTENT(INOUT)                  :: should_stop

     IF(cmd%cmd_id == SWARM_CMD_SHUTDOWN) THEN
        IF(worker%iw>0) WRITE(worker%iw,*) "SWARM| Received shutdown command, quitting."
        should_stop = .TRUE.
     ELSE
        CALL glbopt_worker_execute(worker%glbopt, cmd, report)
     ENDIF

   END SUBROUTINE swarm_worker_execute

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_worker_create_empty_report(worker, report)
       TYPE(swarm_worker_type), INTENT(INOUT)               :: worker
       TYPE(swarm_report_type), INTENT(OUT)                 :: report

       report%worker_id = worker%id
       ALLOCATE(report%positions(3))
       report%positions(:) = 0.0

   END SUBROUTINE swarm_worker_create_empty_report
END MODULE swarm_worker

