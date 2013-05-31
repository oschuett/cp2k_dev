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
  USE glbopt_worker,                   ONLY: glbopt_worker_init,&
                                             glbopt_worker_finalize,&
                                             glbopt_worker_execute,&
                                             glbopt_worker_type
 USE paramopt_worker,                  ONLY: paramopt_worker_init,&
                                             paramopt_worker_finalize,&
                                             paramopt_worker_execute,&
                                             paramopt_worker_type
  USE swarm_message,                   ONLY: swarm_message_type,&
                                             swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_length
  USE input_constants,                 ONLY: swarm_do_glbopt,&
                                             swarm_do_paramopt
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_worker'

 PUBLIC :: swarm_worker_init, swarm_worker_finalize
 PUBLIC :: swarm_worker_execute
 PUBLIC :: swarm_worker_type
 PUBLIC :: swarm_worker_create_initial_report

  TYPE swarm_worker_type
   PRIVATE
   INTEGER                                  :: id = -1
   INTEGER                                  :: iw = -1
   TYPE(cp_error_type)                      :: error
   INTEGER                                  :: behavior = -1
   TYPE(glbopt_worker_type), POINTER        :: glbopt => Null()
   TYPE(paramopt_worker_type), POINTER      :: paramopt => Null()
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

    CALL section_vals_val_get(root_section,"SWARM%BEHAVIOR", i_val=worker%behavior, error=error)

    SELECT CASE(worker%behavior)
      CASE(swarm_do_glbopt)
         ALLOCATE(worker%glbopt)
         CALL glbopt_worker_init(worker%glbopt, para_env, root_section,&
                 input_path, worker_id, worker%iw, error)
      CASE(swarm_do_paramopt)
         ALLOCATE(worker%paramopt)
         CALL paramopt_worker_init(worker%paramopt, para_env, root_section,&
                 input_path, worker_id, worker%iw, error)
      CASE DEFAULT
         STOP "swarm_worker_init: got unkown behavior"
    END SELECT

  END SUBROUTINE swarm_worker_init



! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_worker_finalize(worker)
    TYPE(swarm_worker_type), INTENT(INOUT)  :: worker

     SELECT CASE(worker%behavior)
      CASE(swarm_do_glbopt)
         CALL glbopt_worker_finalize(worker%glbopt)
         DEALLOCATE(worker%glbopt)
      CASE(swarm_do_paramopt)
         CALL paramopt_worker_finalize(worker%paramopt)
         DEALLOCATE(worker%glbopt)
      CASE DEFAULT
         STOP "swarm_worker_finalize: got unkown behavior"
    END SELECT

   END SUBROUTINE swarm_worker_finalize


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_worker_execute(worker, cmd, report, should_stop)
     TYPE(swarm_worker_type), INTENT(INOUT)  :: worker
     TYPE(swarm_message_type), INTENT(IN)    :: cmd
     TYPE(swarm_message_type), INTENT(OUT)   :: report
     LOGICAL, INTENT(INOUT)                  :: should_stop

     CHARACTER(LEN=default_string_length)    :: command

     CALL swarm_message_get(cmd, "command", command)
     CALL swarm_message_add(report, "worker_id", worker%id)

     IF(TRIM(command) == "shutdown") THEN
        IF(worker%iw>0) WRITE(worker%iw,*) "SWARM| Received shutdown command, quitting."
        should_stop = .TRUE.
     ELSE
        SELECT CASE(worker%behavior)
           CASE(swarm_do_glbopt)
             CALL glbopt_worker_execute(worker%glbopt, cmd, report)
           CASE(swarm_do_paramopt)
              CALL paramopt_worker_execute(worker%paramopt, cmd, report)
           CASE DEFAULT
              STOP "swarm_worker_execute: got unkown behavior"
        END SELECT
     ENDIF

   END SUBROUTINE swarm_worker_execute

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE swarm_worker_create_initial_report(worker, report)
       TYPE(swarm_worker_type), INTENT(INOUT)                :: worker
       TYPE(swarm_message_type), INTENT(OUT)                 :: report
       REAL(KIND=dp), DIMENSION(:), POINTER :: fake_positions

       CALL swarm_message_add(report, "worker_id", worker%id)
       CALL swarm_message_add(report, "initial_hello", 1)
       !CALL swarm_message_add(report, "iframe", 1)
       !
       !ALLOCATE(fake_positions(3))
       !fake_positions = 0.0
       !CALL swarm_message_add(report, "positions", fake_positions)
       !WRITE (*,*) "swarm_worker_create_empty_report: ugly hack using fake positions"

   END SUBROUTINE swarm_worker_create_initial_report
END MODULE swarm_worker

