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
  USE force_env_types,                 ONLY: force_env_type,&
                                             force_env_get
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set,&
                                             section_vals_get
  USE kinds,                           ONLY: dp, dp_size, int_4,&
                                             default_path_length,&
                                             default_string_length

  USE timings,                         ONLY: timeset,&
                                             timestop
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_para_env,                     ONLY: cp_para_env_create, cp_para_env_retain
  USE message_passing,                 ONLY: &
       mp_any_source, mp_bcast, mp_comm_dup, mp_comm_free, mp_comm_split, &
       mp_comm_split_direct, mp_environ, mp_recv, mp_send, mp_sum, mp_sync
 
       
  USE f77_interface,                   ONLY: create_force_env,&
                                             destroy_force_env,&
                                             f77_default_para_env => default_para_env,&
                                             f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE geo_opt,                         ONLY: cp_geo_opt
  USE input_cp2k,                      ONLY: create_cp2k_input_reading,&
                                             empty_initial_variables
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_env_create,&
                                             md_env_release,&
                                             md_environment_type,&
                                             set_md_env
  USE md_run,                          ONLY: qs_mol_dyn
  USE parallel_rng_types,              ONLY: reset_to_next_rng_substream
  USE mdctrl_types,                    ONLY: mdctrl_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type,&
                                             cp_subsys_get,&
                                             pack_subsys_particles,&
                                             unpack_subsys_particles
  USE colvar_methods,                  ONLY: colvar_eval_glob_f
  USE hash_functions,                  ONLY: b3hs_hash_key_jenkins
  USE glbopt_types,                    ONLY: glbopt_config_type, glbopt_mdctrl_data_type,&
       glbopt_msg_minima_new, glbopt_msg_minima_old, glbopt_msg_shutdown
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
   SUBROUTINE glbopt_drive_master(glbopt_cfg, para_env, n_colvars, error)
      TYPE(glbopt_config_type), INTENT(in)     :: glbopt_cfg
      TYPE(cp_para_env_type), POINTER          :: para_env
      INTEGER, INTENT(in)                      :: n_colvars
      TYPE(cp_error_type), INTENT(inout)       :: error
      
      ! ========= local vars ===========  
      REAL(kind=dp), dimension(n_colvars + glbopt_cfg%n_hashes) :: key_parts
      INTEGER(kind=int_4), dimension(glbopt_cfg%n_hashes)  :: hashes                             
      CHARACTER, dimension(:), allocatable     :: key
      INTEGER(kind=int_4), dimension(glbopt_cfg%bloom_size/(4*8))       :: bloomfilter
      INTEGER                                  :: key_length
      INTEGER                                  :: tag, calling_walker, i, foo
      LOGICAL ::  was_set, match
      INTEGER :: msg, minima_counter, hop_counter, shutdown_walkers
      
      key_length = size(transfer(key_parts, key)) 
      ALLOCATE(key(key_length))
      bloomfilter = 0
      minima_counter = 0
      hop_counter = 0
      shutdown_walkers = 0
      
      do while(.TRUE.)
        hop_counter = hop_counter + 1
        calling_walker = mp_any_source
      	tag = 42
      	key_parts = 0.0_dp
      	hashes = 0
      	CALL mp_recv(key_parts, calling_walker, tag, para_env%group)
      	!WRITE(*,*) "Master: received new colvars from ", calling_walker
        !WRITE(*,*) "parts: ", key_parts
      	! following advise from:
      	! http://spyced.blogspot.ch/2009/01/all-you-ever-wanted-to-know-about.html
      	! Encode the colvar values as a character array and calculate hash
      	DO i=1, glbopt_cfg%n_hashes
          key = transfer(key_parts, key)    
          hashes(i) = b3hs_hash_key_jenkins(key, glbopt_cfg%bloom_size) - 1
          IF(hashes(i) > glbopt_cfg%bloom_size) STOP "hash value too large"
          IF(hashes(i) < 0 ) STOP "hash value too small"
          key_parts(n_colvars + i) = hashes(i)
        END DO
        
        !write (*,*) "hahshes: ",hashes
        
        ! lockup and register minima at the same time 
        match = .TRUE.
        DO i=1, glbopt_cfg%n_hashes
          was_set = bitfield_set(hashes(i), bloomfilter)
          match = match .AND. was_set
        END DO
        
        if(match) then
          msg = glbopt_msg_minima_old
        else
          msg = glbopt_msg_minima_new
          minima_counter = minima_counter + 1 
        end if
        
        if(hop_counter > glbopt_cfg%hop_count_max) then
          msg = glbopt_msg_shutdown
          shutdown_walkers = shutdown_walkers + 1 
        end if
        
        write (*,*) "minima_counter: ", minima_counter, "hop_counter: ",hop_counter
        
        CALL mp_send(msg, calling_walker, 42, para_env%group)
        
        
        if(shutdown_walkers == glbopt_cfg%n_walkers) then
           write (*,*) "Master: shutdown all walkers - exiting main loop."
           exit
        end if
        
      END DO
      
      DEALLOCATE(key)
      
   END SUBROUTINE glbopt_drive_master
 
  
! *****************************************************************************
!> \brief 
! *****************************************************************************
  FUNCTION bitfield_set(idx, bitfield) RESULT(was_set)
     INTEGER(kind=int_4), INTENT(in)                     :: idx
     INTEGER(kind=int_4), dimension(:), INTENT(inout)    :: bitfield
     LOGICAL :: was_set
     !local vars
     INTEGER                                 :: offset
     INTEGER                                 :: val, mask
     ! Caution: This code only works with 4 byte integers
     
     offset = ISHFT(idx, -5) + 1
     val = IAND(idx, 31)
     mask = 2**val
     
     IF(offset > size(bitfield)) THEN
     	  write (*,*) "idx: ", idx, "offset: ",offset, "size(bitfield): ",size(bitfield)
     	  STOP "set_bitfield: overflow"
     end if
     
     IF(offset < 1) THEN
     	  write (*,*) "idx: ", idx, "offset: ",offset, "size(bitfield): ",size(bitfield)
     	  STOP "set_bitfield: underflow"
     end if
     
     was_set = IAND(bitfield(offset), mask) == mask
     bitfield(offset) = IOR(bitfield(offset), mask)
  END FUNCTION bitfield_set
  
 
 
END MODULE glbopt_master

