!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_master
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE glbopt_types,                    ONLY: GLBOPT_CMD_SHUTDOWN,&
                                             GLBOPT_CMD_MDGOPT,&
                                             glbopt_command_type,&
                                             glbopt_report_type
  USE glbopt_history,                  ONLY: glbopt_history_type,&
                                             glbopt_history_init,&
                                             glbopt_history_finalize,&
                                             glbopt_history_add,&
                                             glbopt_history_lookup
  USE hash_functions,                  ONLY: b3hs_hash_key_jenkins
  USE kinds,                           ONLY: dp,&
                                             int_4
  USE message_passing,                 ONLY: mp_any_source,mp_any_tag,&
                                             mp_recv,&
                                             mp_send
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_type
  USE physcon,                         ONLY: kelvin, femtoseconds

  USE glbopt_minhop,                   ONLY: glbopt_minhop_type,&
                                             glbopt_minhop_init,&
                                             glbopt_minhop_finalize,&
                                             glbopt_minhop_steer
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
 USE timings_types,                   ONLY: timer_env_type
 USE timings,                         ONLY: add_timer_env, rm_timer_env, timer_env_create

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_master'


 PUBLIC :: glbopt_master_type
 PUBLIC :: glbopt_master_init, glbopt_master_finalize
 PUBLIC :: glbopt_master_steer


 TYPE glbopt_master_type
   PRIVATE
   TYPE(glbopt_minhop_type), ALLOCATABLE :: minhop
   REAL(KIND=dp)                                       :: Elowest = Huge(1.0_dp)
   REAL(KIND=dp)                                       :: Emin = Tiny(1.0_dp)
   INTEGER                                             :: iw = 0
   INTEGER                                             :: i_iteration = 0
   INTEGER                                             :: max_iter = 0
   LOGICAL                                             :: should_stop = .FALSE.
   INTEGER                                             :: n_walkers = -1
 END TYPE glbopt_master_type

 CONTAINS



! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_init(master, para_env, root_section, n_walkers, error)
    TYPE(glbopt_master_type)                 :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: n_walkers
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: output_unit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(timer_env_type), POINTER            :: timer_env
    TYPE(section_vals_type), POINTER         :: glbopt_section

    NULLIFY(timer_env)
    CALL timer_env_create(timer_env)
    CALL add_timer_env(timer_env)
    !CALL add_mp_perf_env(f_env%mp_perf_env) TODO: this as well?

    glbopt_section => section_vals_get_subs_vals(root_section, "GLOBAL_OPT", error=error)

    logger => cp_error_get_logger(error)
    master%n_walkers = n_walkers
    master%iw = cp_print_key_unit_nr(logger, glbopt_section,&
          "PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    CALL section_vals_val_get(glbopt_section,"E_MIN", r_val=master%Emin, error=error)
    CALL section_vals_val_get(glbopt_section,"MAX_ITER", i_val=master%max_iter, error=error)

    ALLOCATE(master%minhop)
    CALL glbopt_minhop_init(master%minhop, glbopt_section, n_walkers, master%iw, error)

 END SUBROUTINE glbopt_master_init


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_steer(master, report, cmd)
    TYPE(glbopt_master_type)   :: master
    TYPE(glbopt_report_type)   :: report
    TYPE(glbopt_command_type)  :: cmd

    master%Elowest = MIN(master%Elowest, report%Epot)
    master%i_iteration = master%i_iteration + 1
    IF(master%iw>0) WRITE(master%iw,'(A,I3,A,E15.8,A,E15.8)') &
       " GLOBAL_OPT| Walker ", report%walker_id, " reported minima with Epot=",&
       report%Epot," Elowest=",master%Elowest

    IF(master%Elowest<master%Emin .AND. .NOT. master%should_stop) THEN
       CALL print(master, "Reached Emin > Epot. Quitting.")
       master%should_stop = .TRUE.
    ENDIF

    IF(master%i_iteration>=master%max_iter .AND. .NOT. master%should_stop) THEN
       CALL print(master, "Reached MAX_ITER. Quitting.")
       master%should_stop = .TRUE.
    ENDIF

    CALL glbopt_minhop_steer(master%minhop, report, cmd)

    IF(master%should_stop) THEN
       !CALL print(master, "Sending shutdown command")
       cmd%id = GLBOPT_CMD_SHUTDOWN  !overwrite command
    ENDIF

 END SUBROUTINE


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_finalize(master)
    TYPE(glbopt_master_type)                 :: master

    CALL glbopt_minhop_finalize(master%minhop)
    DEALLOCATE(master%minhop)

    CALL rm_timer_env() !pops the top-most timer
 END SUBROUTINE glbopt_master_finalize

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE print(master, message)
      TYPE(glbopt_master_type), INTENT(INOUT)         :: master
      CHARACTER(LEN=*)                         :: message
      IF (master%iw>0) &
        WRITE (master%iw,*) "GLOBAL_OPT| "//TRIM(message)
   END SUBROUTINE print

END MODULE glbopt_master

