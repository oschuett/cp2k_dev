!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_walker
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

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_walker'

 PUBLIC ::glbopt_drive_walker

  
 CONTAINS 
 
 
! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE glbopt_drive_walker(glbopt_cfg, f_env_id, n_colvars, para_env, walker_group, my_walker_id, error)
    TYPE(glbopt_config_type), INTENT(IN)        :: glbopt_cfg
    INTEGER, INTENT(in)                      :: n_colvars
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: walker_group
    INTEGER, INTENT(IN)                      :: my_walker_id
    TYPE(cp_error_type), INTENT(inout)       :: error
  
    ! ====== Local Variables ======
    CHARACTER(len=*), PARAMETER :: routineN = 'run_walker', &
       routineP = moduleN//':'//routineN
    LOGICAL                                  :: failure
    
    INTEGER                                  :: handle, output_unit, ierr,&
                                                f_env_id, i, n_atoms
    TYPE(mdctrl_type), TARGET                :: mdctrl
    TYPE(mdctrl_type), POINTER               :: mdctrl_p
    TYPE(glbopt_mdctrl_data_type), TARGET    :: mdctrl_data
    
    TYPE(cp_logger_type), POINTER            :: logger
        
  
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    
    
    TYPE(cp_subsys_type), POINTER            :: subsys
    REAL(kind=dp), dimension(n_colvars)      :: colvar_values
    
    REAL(KIND=dp), DIMENSION(:), allocatable :: latest_positions
           
    REAL(KIND=dp), DIMENSION(glbopt_cfg%epot_ndown + glbopt_cfg%epot_nup + 1),&
                                      TARGET :: epot_history
    INTEGER :: tag, msg, master_id
    !LOGICAL :: known_minima
    
    failure = .FALSE.
    
    ! ====== Usual Subroutine Setup ======
    CALL timeset(routineN,handle)
    
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_io_unit(logger)
    
        	
    CALL f_env_add_defaults(f_env_id, f_env, error, failure)
    force_env => f_env%force_env
    
    WRITE(*,*) "allocating epot_history und co"
    CALL force_env_get(force_env, globenv=globenv, subsys=subsys, error=error)
    CALL cp_subsys_get(subsys, natom=n_atoms, error=error)
    ALLOCATE(latest_positions(1:3*n_atoms))
    
    
    
    WRITE(*,*) "setting up random number streams"
    
    ! We want different random-number-streams for each walker
    DO i = 1, my_walker_id
       CALL reset_to_next_rng_substream(globenv%gaussian_rng_stream, error)
    END DO

    
    WRITE(*,*) "setting up md callback"
    mdctrl_data%cfg = glbopt_cfg
    mdctrl_data%epot_history => epot_history
    mdctrl%glbopt => mdctrl_data
    mdctrl_p => mdctrl
    
    WRITE(*,*) "setup done"
    DO while(.TRUE.) !main-loop of minima hopping
      ! store starting configuration
      CALL pack_subsys_particles(subsys, r=latest_positions, error=error)
      mdctrl_data%epot_history = 0.0 
      mdctrl_data%md_mincounter = 0
      !WRITE(*,*) "starting md"
      
      CALL qs_mol_dyn(force_env, globenv, mdctrl=mdctrl_p, error=error)
      !WRITE(*,*) "starting geoopt"
      CALL cp_geo_opt(force_env, globenv, error=error)
      ! TODO: find out if minimazation was successfull
      ! TODO: find out if its a new minima
      
      !WRITE(*,*) "collecting colvars"
      colvar_values = 0.0_dp
      DO i=1, n_colvars
         CALL colvar_eval_glob_f(i,force_env,error=error)
         !WRITE (*,*) "Colvar ", i, "value: ", subsys%colvar_p(i)%colvar%ss
         colvar_values(i) = subsys%colvar_p(i)%colvar%ss
      END DO
      
      master_id = para_env%num_pe-1
      !WRITE(*,*) "Walker: sending to master", my_walker_id
      CALL mp_send(colvar_values,master_id,42,para_env%group)
      !WRITE(*,*) "Walker: waiting for answer", my_walker_id
      tag = 42
      CALL mp_recv(msg, master_id, tag, para_env%group)
      !WRITE(*,*) "Walker: got msg from master: ", msg, my_walker_id
      
      
      SELECT CASE (msg)
       CASE (glbopt_msg_minima_new)
          ! just continue with next loop cylce
       CASE (glbopt_msg_minima_old)
          ! reset to lastest starting positions
          CALL unpack_subsys_particles(subsys, r=latest_positions, error=error)
       CASE (glbopt_msg_shutdown)
          write (*,*) "walker: received shutdown signal ", my_walker_id
      	  exit
       CASE DEFAULT
          STOP "Walker: received unkown message"
      END SELECT
      
      
    END DO
     
    
    
    ! ====== Clean up ======
    
    DEALLOCATE(latest_positions)
    	    
    CALL f_env_rm_defaults(f_env, error,ierr)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    CALL destroy_force_env(f_env_id,ierr)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)


  END SUBROUTINE glbopt_drive_walker
 
 
END MODULE glbopt_walker

