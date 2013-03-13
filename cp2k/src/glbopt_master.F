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
                                             glbopt_msg_shutdown,&
                                             glbopt_msg_sync,&
                                             glbopt_msg_report
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

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_master'

 PUBLIC :: glbopt_drive_master

 CONTAINS


! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE glbopt_drive_master(glbopt_cfg, root_section, para_env, error)
    TYPE(glbopt_config_type), INTENT(in)     :: glbopt_cfg
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

   TYPE(cp_logger_type), POINTER            :: logger
   INTEGER                                  :: calling_walker, tag, iw

   ! CHARACTER, ALLOCATABLE, DIMENSION(:)     :: key
   !hop_counter, &
   !                                             i, key_length, &
   !                                             minima_counter, msg, &
   !                                             shutdown_walkers, tag
   ! LOGICAL                                  :: match, was_set
   ! REAL(kind=dp), DIMENSION(&
   !   n_colvars+glbopt_cfg%n_hashes)         :: key_parts  
   ! INTEGER(kind=int_4), &
   !   DIMENSION(glbopt_cfg%n_hashes)         :: hashes
   INTEGER(kind=int_4), DIMENSION(glbopt_cfg%bloom_size/(4*8)) :: shared_bitfield
   INTEGER(kind=int_4), DIMENSION(glbopt_cfg%bloom_size/(4*8)) :: recv_bitfield
   REAL(kind=dp)                            :: Epot, Elowest

   Elowest = HUGE(0.0_dp)

    !NULLIFY(logger)
    !logger => cp_error_get_logger(error)
    !iw = cp_print_key_unit_nr(logger,root_section,&
    !    "GLOBAL_OPT%PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    iw = 6


    IF (iw>0)  WRITE (iw,'(A)') " GLOBAL_OPT Master| standing by"
    DO WHILE(.TRUE.)
      ! the first message from walkers always contains their bitfield
      calling_walker = mp_any_source
      tag = mp_any_tag
      CALL mp_recv(recv_bitfield, calling_walker, tag, para_env%group)
      shared_bitfield = IOR(shared_bitfield, recv_bitfield)

      SELECT CASE (tag)
        CASE (glbopt_msg_sync)
          ! walker wants to sync - send shared bitfield
          CALL mp_send(shared_bitfield, calling_walker, glbopt_msg_sync, para_env%group)
          !IF (iw>0)  WRITE (iw,'(A, I10)') " GLOBAL_OPT Master| synced with walker ", calling_walker\

        CASE (glbopt_msg_report)
          ! walker reports a minimum
          CALL mp_recv(Epot, calling_walker, tag, para_env%group)
          Elowest = MIN(Elowest, Epot)
          IF (iw>0)  WRITE (iw,'(A,I3,A,E15.8,A,E15.8)') ' GLOBAL_OPT Master| Walker ',&
            calling_walker, ' reported minima with Epot=',Epot,'  Elowest=',Elowest
        CASE DEFAULT
          STOP "glbopt_drive_master: unkown tag."
      END SELECT


    END DO

!      key_length = SIZE(TRANSFER(key_parts, key)) 
!      ALLOCATE(key(key_length))
!      bloomfilter = 0
!      minima_counter = 0
!      hop_counter = 0
!      shutdown_walkers = 0
!
!      DO WHILE(.TRUE.)
!        hop_counter = hop_counter + 1
!        calling_walker = mp_any_source
!              tag = 42
!              key_parts = 0.0_dp
!              hashes = 0
!              CALL mp_recv(key_parts, calling_walker, tag, para_env%group)
!              !WRITE(*,*) "Master: received new colvars from ", calling_walker
!        !WRITE(*,*) "parts: ", key_parts
!              ! following advise from:
!              ! http://spyced.blogspot.ch/2009/01/all-you-ever-wanted-to-know-about.html
!              ! Encode the colvar values as a character array and calculate hash
!              DO i=1, glbopt_cfg%n_hashes
!          key = TRANSFER(key_parts, key)    
!          hashes(i) = b3hs_hash_key_jenkins(key, glbopt_cfg%bloom_size) - 1
!          IF(hashes(i) > glbopt_cfg%bloom_size) STOP "hash value too large"
!          IF(hashes(i) < 0 ) STOP "hash value too small"
!          key_parts(n_colvars + i) = hashes(i)
!        END DO
!
!        !write (*,*) "hahshes: ",hashes
!
!        ! lockup and register minima at the same time 
!        match = .TRUE.
!        DO i=1, glbopt_cfg%n_hashes
!          was_set = bitfield_set(hashes(i), bloomfilter)
!          match = match .AND. was_set
!        END DO
!        
!        IF(match) THEN
!          msg = glbopt_msg_minima_old
!        ELSE
!          msg = glbopt_msg_minima_new
!          minima_counter = minima_counter + 1 
!        END IF
!        
!        IF(hop_counter > glbopt_cfg%hop_count_max) THEN
!          msg = glbopt_msg_shutdown
!          shutdown_walkers = shutdown_walkers + 1 
!        END IF
!        
!        WRITE (*,*) "minima_counter: ", minima_counter, "hop_counter: ",hop_counter
!        
!        CALL mp_send(msg, calling_walker, 42, para_env%group)
!        
!        
!        IF(shutdown_walkers == glbopt_cfg%n_walkers) THEN
!           WRITE (*,*) "Master: shutdown all walkers - exiting main loop."
!           EXIT
!        END IF
!        
!      END DO
!      
!      DEALLOCATE(key)
!      
   END SUBROUTINE glbopt_drive_master
 
  

  
 
 
END MODULE glbopt_master

