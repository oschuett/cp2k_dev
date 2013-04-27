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
  USE glbopt_history,                  ONLY: glbopt_history_add,&
                                             glbopt_history_finalize,&
                                             glbopt_history_init,&
                                             glbopt_history_lookup,&
                                             glbopt_history_type
  USE glbopt_types,                    ONLY: glbopt_cmd_mdgopt,&
                                             glbopt_command_type,&
                                             glbopt_report_type
  USE glbopt_utils,                    ONLY: particles_connected
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: angstrom,&
                                             kelvin
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
   REAL(KIND=dp)                                       :: beta2 = 0
   REAL(KIND=dp)                                       :: beta3 = 0
   REAL(KIND=dp)                                       :: Eaccept0 = 0
   REAL(KIND=dp)                                       :: temp_init = 0
   REAL(KIND=dp)                                       :: temp_max = 0
   REAL(KIND=dp)                                       :: temp_min = 0
   REAL(KIND=dp)                                       :: alpha1 = 0
   REAL(KIND=dp)                                       :: alpha2 = 0
   REAL(KIND=dp)                                       :: max_dist = 0
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
    INTEGER, INTENT(IN)                      :: n_walkers, iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, n_histories
    REAL(kind=dp)                            :: temp_in_kelvin, dist_in_angstrom
    TYPE(section_vals_type), POINTER         :: minhop_section

! read input

    minhop_section => section_vals_get_subs_vals(glbopt_section, "MINIMA_HOPPING", error=error)
    CALL section_vals_val_get(minhop_section,"BETA_2",  r_val=minhop%beta2,error=error)
    CALL section_vals_val_get(minhop_section,"BETA_3",  r_val=minhop%beta3,error=error)
    CALL section_vals_val_get(minhop_section,"ALPHA_1", r_val=minhop%alpha1,error=error)
    CALL section_vals_val_get(minhop_section,"ALPHA_2", r_val=minhop%alpha2,error=error)
    CALL section_vals_val_get(minhop_section,"E_ACCEPT_INIT", r_val=minhop%Eaccept0,error=error)
    CALL section_vals_val_get(minhop_section,"TEMPERATURE_INIT", r_val=temp_in_kelvin,error=error)
    minhop%temp_init = temp_in_kelvin / kelvin
    CALL section_vals_val_get(minhop_section,"TEMPERATURE_MAX", r_val=temp_in_kelvin,error=error)
    minhop%temp_max = temp_in_kelvin / kelvin
    CALL section_vals_val_get(minhop_section,"TEMPERATURE_MIN", r_val=temp_in_kelvin,error=error)
    minhop%temp_min = temp_in_kelvin / kelvin
    CALL section_vals_val_get(minhop_section,"SHARE_HISTORY", l_val=minhop%share_history,error=error)
    CALL section_vals_val_get(minhop_section,"MAX_DIST", r_val=dist_in_angstrom, error=error)
    minhop%max_dist = dist_in_angstrom/angstrom

    IF(minhop%temp_max < minhop%temp_init .OR. minhop%temp_init < minhop%temp_min) &
       STOP "minhop: temp_max < temp_init < temp_min not fulfilled"

    ! allocate history / histories
    n_histories = n_walkers
    IF(minhop%share_history) n_histories = 1
    ALLOCATE(minhop%history(n_histories))
    DO i=1, n_histories
      CALL glbopt_history_init(minhop%history(i), iw=iw)
    END DO

    ALLOCATE(minhop%prev_report(n_walkers))
    ALLOCATE(minhop%temperature(n_walkers))
    ALLOCATE(minhop%Eaccept(n_walkers))

    minhop%temperature(:) = minhop%temp_init
    minhop%Eaccept(:)     = minhop%Eaccept0
    minhop%n_walkers = n_walkers
    minhop%iw = iw
    minhop%n_accepted = 0

    IF(minhop%iw>0) &
      WRITE(minhop%iw, '(A,36X,L1)') " GLOBAL_OPT|MINHOP| Using shared histories:", minhop%share_history
 END SUBROUTINE glbopt_minhop_init

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_minhop_steer(minhop, report, cmd)
    TYPE(glbopt_minhop_type)                 :: minhop
    TYPE(glbopt_report_type)                 :: report
    TYPE(glbopt_command_type)                :: cmd

    INTEGER                                  :: hid, wid
    LOGICAL                                  :: minima_known

    wid = report%walker_id

    IF(minhop%prev_report(wid)%walker_id<0) minhop%prev_report(wid)=report !init

    hid = wid ! history_id = walker_id unless ....
    IF(minhop%share_history) hid = 1 !...there is ONE shared history.

    CALL glbopt_history_lookup(minhop%history(hid), report, minima_known)

    ! Stop if cluster starts to splits into multiple fragments.
    IF(.NOT. particles_connected(report%positions, minhop%max_dist)) THEN
       IF(minhop%iw > 0) WRITE (minhop%iw,"(A)") " GLOBAL_OPT|MINHOP| Conclusion: rejecting minima because of MAX_DIST."
       minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta3  !decreasing temperature

    ELSE IF(minima_known) THEN ! known minima restart MD
       IF(minhop%iw>0) WRITE(minhop%iw,'(A,40X,A)') " GLOBAL_OPT|MINHOP| Conclusion","old minima"
       minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta2 !increasing temperature

    ELSE IF (report%Epot - minhop%prev_report(wid)%Epot > minhop%Eaccept(wid)) THEN
       IF(minhop%iw>0) WRITE(minhop%iw,'(A,7X,A)')&
          " GLOBAL_OPT|MINHOP| Conclusion","new minima, but got rejected due to Eaccept"
       minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta3 !decreasing temperature
       minhop%Eaccept(wid)     = minhop%Eaccept(wid) * minhop%alpha2 !increasing temperature

    ELSE ! Hurra, we accept a new minimum
      IF(minhop%iw>0) WRITE(minhop%iw,'(A,27X,A)') " GLOBAL_OPT|MINHOP| Conclusion","new minima got accepted"
      minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta3 !decreasing temperature
      minhop%Eaccept(wid) = minhop%Eaccept(wid) * minhop%alpha1 !decreasing temperature
      CALL glbopt_history_add(minhop%history(hid), report)
      minhop%prev_report(wid) = report
      minhop%n_accepted = minhop%n_accepted + 1
    END IF


    IF(minhop%temperature(wid) > minhop%temp_max) THEN
       minhop%temperature(wid) = minhop%temp_max
      IF(minhop%iw > 0) WRITE(minhop%iw,*) " GLOBAL_OPT|MINHOP| Throttleing to Max Temperature"
    END IF

    IF(minhop%temperature(wid) < minhop%temp_min) THEN
       minhop%temperature(wid) = minhop%temp_min
      IF(minhop%iw > 0) WRITE(minhop%iw,*) " GLOBAL_OPT|MINHOP| Throttleing to Min Temperature"
    END IF

    IF(minhop%iw>0) THEN
        WRITE (minhop%iw,'(A,4X,E20.10)')&
         " GLOBAL_OPT|MINHOP| Walker's acceptance Energy [Hartree]",minhop%Eaccept(wid)
        WRITE (minhop%iw,'(A,11X,F20.3)')&
         " GLOBAL_OPT|MINHOP| Walker's temperature [Kelvin]",minhop%temperature(wid)*kelvin
    END IF

    cmd%cmd_id          = GLBOPT_CMD_MDGOPT
    cmd%positions       = minhop%prev_report(wid)%positions
    cmd%temperature     = minhop%temperature(wid)

    IF(minhop%iw > 0) &
      WRITE(minhop%iw,'(A,19X,I10)')&
       " GLOBAL_OPT|MINHOP| Total number of accepted minima",minhop%n_accepted

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

