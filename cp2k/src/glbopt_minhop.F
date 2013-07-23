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
  USE glbopt_history,                  ONLY: history_add,&
                                             history_finalize,&
                                             history_init,&
                                             history_lookup,&
                                             history_type
  USE swarm_message,                   ONLY: swarm_message_type,&
                                             swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_haskey,&
                                             swarm_message_length
  USE glbopt_utils,                    ONLY: particles_connected
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: angstrom,&
                                             kelvin

  USE swarm_message,                   ONLY: swarm_message_type,&
                                             swarm_message_add,&
                                             swarm_message_get
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_minhop'


 PUBLIC :: minhop_type
 PUBLIC :: minhop_init, minhop_finalize
 PUBLIC :: minhop_steer


 TYPE minhop_type
   PRIVATE
   TYPE(history_type), DIMENSION(:),ALLOCATABLE :: history
   TYPE(swarm_message_type), DIMENSION(:), ALLOCATABLE :: prev_report
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
   INTEGER                                             :: n_workers =0
   LOGICAL                                             :: share_history = .FALSE.
 END TYPE minhop_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE minhop_init(minhop, glbopt_section, n_workers, iw, error)
    TYPE(minhop_type)                        :: minhop
    TYPE(section_vals_type), POINTER         :: glbopt_section
    INTEGER, INTENT(IN)                      :: n_workers, iw
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
    n_histories = n_workers
    IF(minhop%share_history) n_histories = 1
    ALLOCATE(minhop%history(n_histories))
    DO i=1, n_histories
      CALL history_init(minhop%history(i), iw=iw)
    END DO

    ALLOCATE(minhop%prev_report(n_workers))
    ALLOCATE(minhop%temperature(n_workers))
    ALLOCATE(minhop%Eaccept(n_workers))

    minhop%temperature(:) = minhop%temp_init
    minhop%Eaccept(:)     = minhop%Eaccept0
    minhop%n_workers = n_workers
    minhop%iw = iw
    minhop%n_accepted = 0

    IF(minhop%iw>0) &
      WRITE(minhop%iw, '(A,47X,L1)') " MINHOP| Using shared histories:", minhop%share_history
 END SUBROUTINE minhop_init


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE minhop_steer(minhop, report, cmd)
    TYPE(minhop_type)                       :: minhop
    TYPE(swarm_message_type)                :: report
    TYPE(swarm_message_type)                :: cmd

    INTEGER                                 :: hid, wid
    LOGICAL                                 :: minima_known

    REAL(KIND=dp)                           :: report_Epot, prev_Epot
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE:: report_positions, prev_positions
    INTEGER                                 :: iframe


    CALL swarm_message_get(report, "worker_id", wid)

    IF(swarm_message_haskey(report, "initial_hello")) THEN
       CALL swarm_message_add(cmd, "command", "md_and_gopt")
       CALL swarm_message_add(cmd, "iframe", 1)
       CALL swarm_message_add(cmd, "temperature", minhop%temperature(wid))
       IF(minhop%iw > 0) WRITE(minhop%iw,'(A,30X,I10)') &
                     " MINHOP| Sending initial command to worker", wid
       RETURN
    ENDIF

    IF(swarm_message_length(minhop%prev_report(wid))==0) &
              minhop%prev_report(wid)=report !init (first real report)

    CALL swarm_message_get(report, "Epot", report_Epot)
    CALL swarm_message_get(report, "positions", report_positions)
    CALL swarm_message_get(minhop%prev_report(wid), "Epot", prev_Epot)

    hid = wid ! history_id = worker_id unless ....
    IF(minhop%share_history) hid = 1 !...there is ONE shared history.

    CALL history_lookup(minhop%history(hid), report_Epot, minima_known)

    ! Stop if cluster starts to splits into multiple fragments.
    IF(.NOT. particles_connected(report_positions, minhop%max_dist)) THEN
       IF(minhop%iw > 0) WRITE (minhop%iw,"(A)") " MINHOP| Conclusion: rejecting minima because of MAX_DIST."
       minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta3  !decreasing temperature

    ELSE IF(minima_known) THEN ! known minima restart MD
       IF(minhop%iw>0) WRITE(minhop%iw,'(A,51X,A)') " MINHOP| Conclusion","old minima"
       minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta2 !increasing temperature

    ELSE IF (report_Epot - prev_Epot > minhop%Eaccept(wid)) THEN
       IF(minhop%iw>0) WRITE(minhop%iw,'(A,18X,A)')&
          " MINHOP| Conclusion","new minima, but got rejected due to Eaccept"
       minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta3 !decreasing temperature
       minhop%Eaccept(wid)     = minhop%Eaccept(wid) * minhop%alpha2 !increasing temperature

    ELSE ! Hurra, we accept a new minimum
      IF(minhop%iw>0) WRITE(minhop%iw,'(A,38X,A)') " MINHOP| Conclusion","new minima got accepted"
      minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta3 !decreasing temperature
      minhop%Eaccept(wid) = minhop%Eaccept(wid) * minhop%alpha1 !decreasing temperature
      CALL history_add(minhop%history(hid), report_Epot)
      minhop%prev_report(wid) = report
      minhop%n_accepted = minhop%n_accepted + 1
    END IF


    IF(minhop%temperature(wid) > minhop%temp_max) THEN
       minhop%temperature(wid) = minhop%temp_max
      IF(minhop%iw > 0) WRITE(minhop%iw,*) " MINHOP| Throttleing to Max Temperature"
    END IF

    IF(minhop%temperature(wid) < minhop%temp_min) THEN
       minhop%temperature(wid) = minhop%temp_min
      IF(minhop%iw > 0) WRITE(minhop%iw,*) " MINHOP| Throttleing to Min Temperature"
    END IF

    IF(minhop%iw>0) THEN
        WRITE (minhop%iw,'(A,15X,E20.10)')&
         " MINHOP| Worker's acceptance Energy [Hartree]",minhop%Eaccept(wid)
        WRITE (minhop%iw,'(A,22X,F20.3)')&
         " MINHOP| Worker's temperature [Kelvin]",minhop%temperature(wid)*kelvin
    END IF


    !TODO write copy routine
    CALL swarm_message_get(report, "iframe", iframe)
    CALL swarm_message_add(cmd, "iframe", iframe)
    CALL swarm_message_add(cmd, "command", "md_and_gopt")
    CALL swarm_message_get(minhop%prev_report(wid), "positions", prev_positions)
    CALL swarm_message_add(cmd, "positions", prev_positions)
    CALL swarm_message_add(cmd, "temperature", minhop%temperature(wid))

    IF(minhop%iw > 0) &
      WRITE(minhop%iw,'(A,30X,I10)')&
       " MINHOP| Total number of accepted minima",minhop%n_accepted

 END SUBROUTINE minhop_steer


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE minhop_finalize(minhop)
    TYPE(minhop_type)                 :: minhop

    INTEGER                                  :: i

    DO i=1, SIZE(minhop%history)
      CALL history_finalize(minhop%history(i))
    END DO
 END SUBROUTINE minhop_finalize



END MODULE glbopt_minhop

