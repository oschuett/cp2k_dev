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


#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_minhop'


 PUBLIC :: glbopt_minhop_type
 PUBLIC :: glbopt_minhop_init, glbopt_minhop_finalize
 PUBLIC :: glbopt_minhop_steer, glbopt_minhop_register_input


 TYPE glbopt_minhop_type
   PRIVATE
   TYPE(glbopt_history_type), DIMENSION(:),ALLOCATABLE :: history
   TYPE(glbopt_report_type), DIMENSION(:), ALLOCATABLE :: prev_report
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: temperature
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: Eaccept
   REAL(KIND=dp)                                       :: beta1 = 1.1
   !REAL(KIND=dp)                                      :: beta2 1.1 !we do not differentiate for this case
   REAL(KIND=dp)                                       :: beta3 = 1.0/1.1
   REAL(KIND=dp)                                       :: Eaccept0 = 0.005
 !  REAL(KIND=dp)                                       :: temperature0 = 100*kelvin
   REAL(KIND=dp)                                       :: alpha1 = 0.98
   REAL(KIND=dp)                                       :: alpha2 = 1/0.98
   INTEGER                                             :: iterations = 0
   INTEGER                                             :: n_accepted = 0
   INTEGER                                             :: iw
   INTEGER                                             :: n_walkers
   LOGICAL                                             :: share_history = .TRUE.
 END TYPE glbopt_minhop_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_minhop_steer(minhop, report, cmd)
    TYPE(glbopt_minhop_type)           :: minhop
    TYPE(glbopt_report_type)           :: report
    TYPE(glbopt_command_type)          :: cmd

    INTEGER :: wid, hid
    LOGICAL :: minima_known


    wid = report%walker_id

    hid = wid ! history_id = walker_id unless ....
    IF(minhop%share_history) hid = 1 !...there is ONE shared history.

    CALL glbopt_history_lookup(minhop%history(hid), report, minima_known)

    IF(minima_known) THEN ! known minima restart MD
       CALL print(minhop,"Conclusion: old minima")
       minhop%temperature(wid) = minhop%temperature(wid)*1.1
       cmd%positions   = minhop%prev_report(wid)%positions

    ELSE IF (report%Epot - minhop%prev_report(wid)%Epot > minhop%Eaccept(wid)) THEN
       CALL print(minhop, "Conclusion: new minima, but got rejected due to Eaccept")
       minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta3
       minhop%Eaccept(wid)     = minhop%Eaccept(wid) * minhop%alpha2

    ELSE ! Hurra, we accept a new minimum
      CALL print(minhop, "Conclusion: new minima got accepted")
      minhop%temperature(wid) = minhop%temperature(wid) * minhop%beta3
      minhop%Eaccept(wid) = minhop%Eaccept(wid) * minhop%alpha1
      CALL glbopt_history_add(minhop%history(hid), report)
      minhop%prev_report = report
      minhop%n_accepted = minhop%n_accepted + 1
    END IF

    cmd%id              = GLBOPT_CMD_MDGOPT
    cmd%positions       = minhop%prev_report(wid)%positions
    cmd%temperature     = minhop%temperature(wid)

    minhop%iterations = minhop%iterations + 1
    IF(minhop%iw > 0) &
    WRITE(minhop%iw,*) "GLOBAL_OPT| Minima Hopping Interation ", minhop%iterations
    !IF(minhop%iw > 0) &
    !WRITE (minhop%iw,'(A,A,A,I9,A,A)') ' GLOBAL_OPT| ',REPEAT("*", 15),' ',&
    !     minhop%iterations,'. Minima Hopping Iteration ',REPEAT("*", 15)

 END SUBROUTINE glbopt_minhop_steer


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_minhop_init(minhop, root_section, n_walkers, iw, error)
    TYPE(glbopt_minhop_type)                 :: minhop
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: n_walkers
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, n_histories
    TYPE(cp_logger_type), POINTER            :: logger

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

    minhop%temperature(:) = 100/kelvin !TODO: make input parameter
    minhop%share_history = .TRUE.
    minhop%Eaccept(:)     = minhop%Eaccept0
    minhop%n_walkers = n_walkers
    minhop%iw = iw

    !DO i=1, cfg%n_walkers
    !  minhop%walkers(i)%temperature = 100*kelvin
    !  minhop%walkers(i)%Eaccept     = cfg%Eaccept0
    !END DO

 END SUBROUTINE glbopt_minhop_init

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_minhop_finalize(minhop)
    TYPE(glbopt_minhop_type)                 :: minhop
    INTEGER                                  :: i

    DO i=1, SIZE(minhop%history)
      CALL glbopt_history_finalize(minhop%history(i))
    END DO
 END SUBROUTINE glbopt_minhop_finalize

 
! *****************************************************************************
! *****************************************************************************
   SUBROUTINE print(minhop, message)
      TYPE(glbopt_minhop_type), INTENT(INOUT)         :: minhop
      CHARACTER(LEN=*)                         :: message
      IF (minhop%iw>0) &
        WRITE (minhop%iw,*) "GLOBAL_OPT| "//TRIM(message)
   END SUBROUTINE print


!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE glbopt_minhop_register_input(glbopt_section,error)
    TYPE(section_type), POINTER              :: glbopt_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(section_type), POINTER              :: section
    TYPE(keyword_type), POINTER              :: keyword

    NULLIFY(section, keyword)
    CALL section_create(section,name="MINHOP",&
         description="blah TODO.",&
         repeats=.FALSE., required=.TRUE., supported_feature=.TRUE.,error=error)
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="EMIN",&
        description="balba TODO.",&
        usage="EMIN <REAL>", &
        type_of_var=real_t,unit_str="hartree",required=.TRUE.,error=error)
    CALL section_add_keyword(section, keyword, error=error)
    CALL keyword_release(keyword, error=error)
    CALL section_add_subsection(glbopt_section,section,error=error)
    CALL section_release(section,error=error)
  END SUBROUTINE glbopt_minhop_register_input

!! *****************************************************************************
!!> \brief Main driver to perform global optimization
!!> \param error variable to control error logging, stopping,...
!!>        see module cp_error_handling
!! *****************************************************************************
!   SUBROUTINE glbopt_drive_minhop(glbopt_cfg, root_section, para_env, error)
!    TYPE(glbopt_config_type), INTENT(in)     :: glbopt_cfg
!    TYPE(section_vals_type), POINTER         :: root_section
!    TYPE(cp_para_env_type), POINTER          :: para_env
!    TYPE(cp_error_type), INTENT(inout)       :: error
!
!   TYPE(cp_logger_type), POINTER            :: logger
!   INTEGER                                  :: calling_walker, tag, iw
!
!   ! CHARACTER, ALLOCATABLE, DIMENSION(:)     :: key
!   !hop_counter, &
!   !                                             i, key_length, &
!   !                                             minima_counter, msg, &
!   !                                             shutdown_walkers, tag
!   ! LOGICAL                                  :: match, was_set
!   ! REAL(kind=dp), DIMENSION(&
!   !   n_colvars+glbopt_cfg%n_hashes)         :: key_parts  
!   ! INTEGER(kind=int_4), &
!   !   DIMENSION(glbopt_cfg%n_hashes)         :: hashes
!   INTEGER(kind=int_4), DIMENSION(glbopt_cfg%bloom_size/(4*8)) :: shared_bitfield
!   INTEGER(kind=int_4), DIMENSION(glbopt_cfg%bloom_size/(4*8)) :: recv_bitfield
!   REAL(kind=dp)                            :: Epot, Elowest
!
!   Elowest = HUGE(0.0_dp)
!
!    !NULLIFY(logger)
!    !logger => cp_error_get_logger(error)
!    !iw = cp_print_key_unit_nr(logger,root_section,&
!    !    "GLOBAL_OPT%PRINT%MINHOP_RUN_INFO",extension=".minhopLog",error=error)
!
!    iw = 6
!
!
!    IF (iw>0)  WRITE (iw,'(A)') " GLOBAL_OPT Minhop| standing by"
!    DO WHILE(.TRUE.)
!      ! the first message from walkers always contains their bitfield
!      calling_walker = mp_any_source
!      tag = mp_any_tag
!      CALL mp_recv(recv_bitfield, calling_walker, tag, para_env%group)
!      shared_bitfield = IOR(shared_bitfield, recv_bitfield)
!
!      SELECT CASE (tag)
!        CASE (glbopt_cmd_sync)
!          ! walker wants to sync - send shared bitfield
!          CALL mp_send(shared_bitfield, calling_walker, glbopt_cmd_sync, para_env%group)
!          !IF (iw>0)  WRITE (iw,'(A, I10)') " GLOBAL_OPT Minhop| synced with walker ", calling_walker\
!
!        CASE (glbopt_cmd_report)
!          ! walker reports a minimum
!          CALL mp_recv(Epot, calling_walker, tag, para_env%group)
!          Elowest = MIN(Elowest, Epot)
!          IF (iw>0)  WRITE (iw,'(A,I3,A,E15.8,A,E15.8)') ' GLOBAL_OPT Minhop| Walker ',&
!            calling_walker, ' reported minima with Epot=',Epot,'  Elowest=',Elowest
!        CASE DEFAULT
!          STOP "glbopt_drive_minhop: unkown tag."
!      END SELECT
!
!
!    END DO
!
!!      key_length = SIZE(TRANSFER(key_parts, key)) 
!!      ALLOCATE(key(key_length))
!!      bloomfilter = 0
!!      minima_counter = 0
!!      hop_counter = 0
!!      shutdown_walkers = 0
!!
!!      DO WHILE(.TRUE.)
!!        hop_counter = hop_counter + 1
!!        calling_walker = mp_any_source
!!              tag = 42
!!              key_parts = 0.0_dp
!!              hashes = 0
!!              CALL mp_recv(key_parts, calling_walker, tag, para_env%group)
!!              !WRITE(*,*) "Minhop: received new colvars from ", calling_walker
!!        !WRITE(*,*) "parts: ", key_parts
!!              ! following advise from:
!!              ! http://spyced.blogspot.ch/2009/01/all-you-ever-wanted-to-know-about.html
!!              ! Encode the colvar values as a character array and calculate hash
!!              DO i=1, glbopt_cfg%n_hashes
!!          key = TRANSFER(key_parts, key)    
!!          hashes(i) = b3hs_hash_key_jenkins(key, glbopt_cfg%bloom_size) - 1
!!          IF(hashes(i) > glbopt_cfg%bloom_size) STOP "hash value too large"
!!          IF(hashes(i) < 0 ) STOP "hash value too small"
!!          key_parts(n_colvars + i) = hashes(i)
!!        END DO
!!
!!        !write (*,*) "hahshes: ",hashes
!!
!!        ! lockup and register minima at the same time 
!!        match = .TRUE.
!!        DO i=1, glbopt_cfg%n_hashes
!!          was_set = bitfield_set(hashes(i), bloomfilter)
!!          match = match .AND. was_set
!!        END DO                                                                                     
!!        
!!        IF(match) THEN
!!          msg = glbopt_msg_minima_old
!!        ELSE                                                                                                                   
!!          msg = glbopt_msg_minima_new
!!          minima_counter = minima_counter + 1 
!!        END IF
!!        
!!        IF(hop_counter > glbopt_cfg%hop_count_max) THEN
!!          msg = glbopt_msg_shutdown
!!          shutdown_walkers = shutdown_walkers + 1 
!!        END IF
!!        
!!        WRITE (*,*) "minima_counter: ", minima_counter, "hop_counter: ",hop_counter
!!        
!!        CALL mp_send(msg, calling_walker, 42, para_env%group)
!!        
!!        
!!        IF(shutdown_walkers == glbopt_cfg%n_walkers) THEN
!!           WRITE (*,*) "Minhop: shutdown all walkers - exiting main loop."
!!           EXIT
!!        END IF
!!        
!!      END DO
!!      
!!      DEALLOCATE(key)
!!      
!   END SUBROUTINE glbopt_drive_minhop
 
  

  
 
 
END MODULE glbopt_minhop

