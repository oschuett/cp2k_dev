!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_minhop
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE glbopt_types,                    ONLY: glbopt_cmd_shutdown,&
                                             glbopt_cmd_mdgopt,&
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
  !USE input_constants

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
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_minhop'


 PUBLIC :: glbopt_minhop_type
 PUBLIC :: glbopt_minhop_init, glbopt_minhop_finalize
 PUBLIC :: glbopt_minhop_steer


 TYPE glbopt_minhop_type
   PRIVATE
   TYPE(glbopt_history_type), DIMENSION(:),ALLOCATABLE :: history
   TYPE(glbopt_report_type), DIMENSION(:), ALLOCATABLE :: prev_report
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: temperature
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: Eaccept
   REAL(KIND=dp)                                       :: beta1 = 0
   REAL(KIND=dp)                                       :: beta2 = 0
   REAL(KIND=dp)                                       :: beta3 = 0
   REAL(KIND=dp)                                       :: Eaccept0 = 0
   REAL(KIND=dp)                                       :: temp_init = 0
   REAL(KIND=dp)                                       :: temp_throttle = 0
   REAL(KIND=dp)                                       :: alpha1 = -0
   REAL(KIND=dp)                                       :: alpha2 = 0
   INTEGER                                             :: iterations = 0
   INTEGER                                             :: n_accepted = 0
   INTEGER                                             :: iw = 0
   INTEGER                                             :: n_walkers =0
   LOGICAL                                             :: share_history = .FALSE.
 END TYPE glbopt_minhop_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_minhop_init(minhop, glbopt_section, n_walkers, iw, error)
    TYPE(glbopt_minhop_type)                 :: minhop
    TYPE(section_vals_type), POINTER         :: glbopt_section
    INTEGER, INTENT(IN)                      :: n_walkers
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, n_histories
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: minhop_section
    REAL(kind=dp)                            :: temp_in_kelvin

    ! read input
    minhop_section => section_vals_get_subs_vals(glbopt_section, "MINIMA_HOPPING", error=error)
    CALL section_vals_val_get(minhop_section,"BETA_1",  r_val=minhop%beta1,error=error)
    CALL section_vals_val_get(minhop_section,"BETA_2",  r_val=minhop%beta2,error=error)
    CALL section_vals_val_get(minhop_section,"BETA_3",  r_val=minhop%beta3,error=error)
    CALL section_vals_val_get(minhop_section,"ALPHA_1", r_val=minhop%alpha1,error=error)
    CALL section_vals_val_get(minhop_section,"ALPHA_2", r_val=minhop%alpha2,error=error)
    CALL section_vals_val_get(minhop_section,"E_ACCEPT_INIT", r_val=minhop%Eaccept0,error=error)
    CALL section_vals_val_get(minhop_section,"TEMPERATURE_INIT", r_val=temp_in_kelvin,error=error)
    minhop%temp_init = temp_in_kelvin / kelvin
    CALL section_vals_val_get(minhop_section,"TEMPERATURE_THROTTLE", r_val=temp_in_kelvin,error=error)
    minhop%temp_throttle = temp_in_kelvin / kelvin
    CALL section_vals_val_get(minhop_section,"SHARE_HISTORY", l_val=minhop%share_history,error=error)

    ! allocate history / histories
    n_histories = n_walkers
    IF(minhop%share_history) n_histories = 1
    ALLOCATE(minhop%history(n_histories))
    DO i=1, n_histories
      CALL glbopt_history_init(minhop%history(i))
    END DO

    ALLOCATE(minhop%prev_report(n_walkers))
    ALLOCATE(minhop%temperature(n_walkers))
    ALLOCATE(minhop%Eaccept(n_walkers))

    minhop%temperature(:) = minhop%temp_init
    minhop%Eaccept(:)     = minhop%Eaccept0
    minhop%n_walkers = n_walkers
    minhop%iw = iw
    minhop%iterations = 0
    minhop%n_accepted = 0

    IF(minhop%iw>0) &
      WRITE(minhop%iw, '(A,43X,L1)') " GLOBAL_OPT| Using shared histories:", minhop%share_history
 END SUBROUTINE glbopt_minhop_init

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_minhop_steer(minhop, report, cmd)
    TYPE(glbopt_minhop_type)           :: minhop
    TYPE(glbopt_report_type)           :: report
    TYPE(glbopt_command_type)          :: cmd

    INTEGER :: wid, hid
    LOGICAL :: minima_known


    wid = report%walker_id

    IF(minhop%prev_report(wid)%walker_id<0) minhop%prev_report(wid)=report !init

    hid = wid ! history_id = walker_id unless ....
    IF(minhop%share_history) hid = 1 !...there is ONE shared history.

    CALL glbopt_history_lookup(minhop%history(hid), report, minima_known)

    IF(minima_known) THEN ! known minima restart MD
       IF(minhop%iw>0) WRITE(minhop%iw,*) "GLOBAL_OPT| Conclusion: old minima"
       minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta1

    ELSE IF (report%Epot - minhop%prev_report(wid)%Epot > minhop%Eaccept(wid)) THEN
       IF(minhop%iw>0) WRITE(minhop%iw,*)&
          "GLOBAL_OPT| Conclusion: new minima, but got rejected due to Eaccept"
       minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta3
       minhop%Eaccept(wid)     = minhop%Eaccept(wid) * minhop%alpha2

    ELSE ! Hurra, we accept a new minimum
      IF(minhop%iw>0) WRITE(minhop%iw,*) "GLOBAL_OPT| Conclusion: new minima got accepted"
      minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta3
      minhop%Eaccept(wid) = minhop%Eaccept(wid) * minhop%alpha1
      CALL glbopt_history_add(minhop%history(hid), report)
      minhop%prev_report(wid) = report
      minhop%n_accepted = minhop%n_accepted + 1
    END IF

    IF(minhop%temperature(wid) > minhop%temp_throttle) THEN
       minhop%temperature(wid) = minhop%temp_throttle
      IF(minhop%iw > 0) WRITE(minhop%iw,*) "GLOBAL_OPT| Throttleing Temperature"
    END IF

    cmd%id              = GLBOPT_CMD_MDGOPT
    cmd%positions       = minhop%prev_report(wid)%positions
    cmd%temperature     = minhop%temperature(wid)

    minhop%iterations = minhop%iterations + 1
    IF(minhop%iw > 0) &
      WRITE(minhop%iw,'(A,I10,A,I10)') " GLOBAL_OPT| Minima Hopping Interation ",&
        minhop%iterations, " accepted number of minima: ",minhop%n_accepted

    !IF(minhop%iw > 0) &
    !WRITE (minhop%iw,'(A,A,A,I9,A,A)') ' GLOBAL_OPT| ',REPEAT("*", 15),' ',&
    !     minhop%iterations,'. Minima Hopping Iteration ',REPEAT("*", 15)

 END SUBROUTINE glbopt_minhop_steer



! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_minhop_finalize(minhop)
    TYPE(glbopt_minhop_type)                 :: minhop
    INTEGER                                  :: i

    DO i=1, SIZE(minhop%history)
      CALL glbopt_history_finalize(minhop%history(i))
    END DO
 END SUBROUTINE glbopt_minhop_finalize



END MODULE glbopt_minhop

