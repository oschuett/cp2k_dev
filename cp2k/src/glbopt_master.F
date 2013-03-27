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
                                             glbopt_cmd_shutdown,&
                                             glbopt_cmd_sync,&
                                             glbopt_cmd_report,&
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

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_master'


 PUBLIC :: glbopt_master_type
 PUBLIC :: glbopt_master_init, glbopt_master_finalize
 PUBLIC :: glbopt_master_give_command


 !TYPE walker_state_type
 ! REAL(KIND=dp)                                     :: Eaccept
 ! REAL(KIND=dp)                                     :: temperature
 ! REAL(KIND=dp), DIMENSION(:), ALLOCATABLE          :: prev_positions
 !END TYPE walker_state_type

 TYPE glbopt_master_type
   PRIVATE
   TYPE(glbopt_history_type)                           :: history
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
   REAL(KIND=dp)                                       :: Elowest = Huge(1.0_dp)
   INTEGER                                             :: iterations = 0
   INTEGER                                             :: n_accepted = 0
   INTEGER                                             :: iw
 END TYPE glbopt_master_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_give_command(master, report, cmd)
    TYPE(glbopt_master_type)   :: master
    TYPE(glbopt_report_type)   :: report
    TYPE(glbopt_command_type)  :: cmd
    INTEGER :: wid
    LOGICAL :: minima_known

    master%iterations = master%iterations + 1
    IF(master%iw > 0) &
    WRITE (master%iw,'(A,A,A,I9,A,A)') ' GLOBAL_OPT| ',REPEAT("*", 15),' ',&
         master%iterations,'. Minima Hopping Iteration ',REPEAT("*", 15)

    IF(master%iterations > 30) THEN
       cmd%id = GLBOPT_CMD_SHUTDOWN
       RETURN
    ENDIF

    wid = report%walker_id

    CALL glbopt_history_lookup(master%history, report, minima_known)
    IF(minima_known) THEN ! known minima restart MD
       CALL print(master,"Conclusion: old minima")
       master%temperature(wid) = master%temperature(wid)*1.1
       cmd%positions   = master%prev_report(wid)%positions

    ELSE IF (report%Epot - master%prev_report(wid)%Epot > master%Eaccept(wid)) THEN
       CALL print(master, "Conclusion: new minima, but got rejected due to Eaccept")
       master%temperature(wid) = master%temperature(wid) * master%beta3
       master%Eaccept(wid)     = master%Eaccept(wid) * master%alpha2

    ELSE ! Hurra, we accept a new minimum
      CALL print(master, "Conclusion: new minima got accepted")
      master%temperature(wid) = master%temperature(wid) * master%beta3
      master%Eaccept(wid) = master%Eaccept(wid) * master%alpha1
      CALL glbopt_history_add(master%history, report)
      master%prev_report = report
      master%Elowest = MIN(master%Elowest, report%Epot)
      master%n_accepted = master%n_accepted + 1
    END IF

    cmd%id              = GLBOPT_CMD_MDGOPT
    cmd%positions       = master%prev_report(wid)%positions
    cmd%temperature     = master%temperature(wid)

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

    INTEGER                                  :: i
    TYPE(cp_logger_type), POINTER            :: logger


    CALL glbopt_history_init(master%history)

    ALLOCATE(master%prev_report(cfg%n_walkers))
    ALLOCATE(master%temperature(cfg%n_walkers))
    ALLOCATE(master%Eaccept(cfg%n_walkers))

    master%temperature(:) = 100/kelvin !TODO: make input parameter
    master%Eaccept(:)     = master%Eaccept0

    ! getting an output unit for logging
    logger => cp_error_get_logger(error)
    master%iw = cp_print_key_unit_nr(logger,root_section,&
          "GLOBAL_OPT%PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    !DO i=1, cfg%n_walkers
    !  master%walkers(i)%temperature = 100*kelvin
    !  master%walkers(i)%Eaccept     = cfg%Eaccept0
    !END DO

 END SUBROUTINE glbopt_master_init

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_finalize(master)
    TYPE(glbopt_master_type)                 :: master
    !nothing todo, yet
    CALL glbopt_history_finalize(master%history)
 END SUBROUTINE glbopt_master_finalize

 
! *****************************************************************************
! *****************************************************************************
   SUBROUTINE print(master, message)
      TYPE(glbopt_master_type), INTENT(INOUT)         :: master
      CHARACTER(LEN=*)                         :: message
      IF (master%iw>0) &
        WRITE (master%iw,*) "GLOBAL_OPT| "//TRIM(message)
   END SUBROUTINE print

!! *****************************************************************************
!!> \brief Main driver to perform global optimization
!!> \param error variable to control error logging, stopping,...
!!>        see module cp_error_handling
!! *****************************************************************************
!   SUBROUTINE glbopt_drive_master(glbopt_cfg, root_section, para_env, error)
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
!    !    "GLOBAL_OPT%PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)
!
!    iw = 6
!
!
!    IF (iw>0)  WRITE (iw,'(A)') " GLOBAL_OPT Master| standing by"
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
!          !IF (iw>0)  WRITE (iw,'(A, I10)') " GLOBAL_OPT Master| synced with walker ", calling_walker\
!
!        CASE (glbopt_cmd_report)
!          ! walker reports a minimum
!          CALL mp_recv(Epot, calling_walker, tag, para_env%group)
!          Elowest = MIN(Elowest, Epot)
!          IF (iw>0)  WRITE (iw,'(A,I3,A,E15.8,A,E15.8)') ' GLOBAL_OPT Master| Walker ',&
!            calling_walker, ' reported minima with Epot=',Epot,'  Elowest=',Elowest
!        CASE DEFAULT
!          STOP "glbopt_drive_master: unkown tag."
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
!!              !WRITE(*,*) "Master: received new colvars from ", calling_walker
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
!!           WRITE (*,*) "Master: shutdown all walkers - exiting main loop."
!!           EXIT
!!        END IF
!!        
!!      END DO
!!      
!!      DEALLOCATE(key)
!!      
!   END SUBROUTINE glbopt_drive_master
 
  

  
 
 
END MODULE glbopt_master

