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
                                             glbopt_msg_minima_new,&
                                             glbopt_msg_minima_old,&
                                             glbopt_msg_shutdown
  USE hash_functions,                  ONLY: b3hs_hash_key_jenkins  
  USE kinds,                           ONLY: dp,&
                                             int_4
  USE message_passing,                 ONLY: mp_any_source,&
                                             mp_recv,&
                                             mp_send
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
   SUBROUTINE glbopt_drive_master(glbopt_cfg, para_env, error)
    TYPE(glbopt_config_type), INTENT(in)     :: glbopt_cfg
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

   ! CHARACTER, ALLOCATABLE, DIMENSION(:)     :: key
   ! INTEGER                                  :: calling_walker, hop_counter, &
   !                                             i, key_length, &
   !                                             minima_counter, msg, &
   !                                             shutdown_walkers, tag
   ! LOGICAL                                  :: match, was_set
   ! REAL(kind=dp), DIMENSION(&
   !   n_colvars+glbopt_cfg%n_hashes)         :: key_parts  
   ! INTEGER(kind=int_4), &
   !   DIMENSION(glbopt_cfg%n_hashes)         :: hashes
   ! INTEGER(kind=int_4), &
   !   DIMENSION(glbopt_cfg%bloom_size/(4*8)) :: bloomfilter


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

