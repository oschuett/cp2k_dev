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
  USE glbopt_types,                    ONLY: glbopt_config_type,&
                                             GLBOPT_CMD_SHUTDOWN,&
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
                                             glbopt_minhop_steer,&
                                             glbopt_minhop_register_input
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type


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
   INTEGER                                             :: max_iteration = 0
   LOGICAL                                             :: should_stop = .FALSE.
 END TYPE glbopt_master_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_steer(master, report, cmd)
    TYPE(glbopt_master_type)   :: master
    TYPE(glbopt_report_type)   :: report
    TYPE(glbopt_command_type)  :: cmd

    master%Elowest = MIN(master%Elowest, report%Epot)
    master%i_iteration = master%i_iteration + 1

    IF(master%Elowest<master%Emin .AND. .NOT. master%should_stop) THEN
       CALL print(master, "Reached Emin > Epot. Quitting.")
       master%should_stop = .TRUE.
    END IF

    IF(master%i_iteration>=master%max_iteration .AND. .NOT. master%should_stop) THEN
       CALL print(master, "Reached max iterations. Quitting.")
       master%should_stop = .TRUE.
    END IF

    CALL glbopt_minhop_steer(master%minhop, report, cmd)

    IF(master%should_stop) cmd%id = GLBOPT_CMD_SHUTDOWN  !overwrite command

 END SUBROUTINE


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_init(master, cfg, root_section, input_path, output_unit, para_env, error)
    TYPE(glbopt_master_type)                 :: master
    TYPE(glbopt_config_type)                 :: cfg
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_logger_type), POINTER            :: logger

    ! getting an output unit for logging
    logger => cp_error_get_logger(error)
    master%iw = cp_print_key_unit_nr(logger,root_section,&
          "GLOBAL_OPT%PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    ALLOCATE(master%minhop)
    CALL glbopt_minhop_init(master%minhop, cfg, root_section, input_path, output_unit, para_env, error)
    master%Emin = cfg%Emin !TODO: read directly from glbopt_section
    master%max_iteration = 3
 END SUBROUTINE glbopt_master_init

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_finalize(master)
    TYPE(glbopt_master_type)                 :: master

    CALL glbopt_minhop_finalize(master%minhop)
    DEALLOCATE(master%minhop)
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

