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
  USE colvar_methods,                  ONLY: colvar_eval_glob_f
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type,&
                                             pack_subsys_particles,&
                                             unpack_subsys_particles
  USE f77_interface,                   ONLY: destroy_force_env,&
                                             f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE geo_opt,                         ONLY: cp_geo_opt
  USE glbopt_types,                    ONLY: glbopt_config_type,&
                                             glbopt_mdctrl_data_type,&
                                             glbopt_msg_minima_new,&
                                             glbopt_msg_minima_old,&
                                             glbopt_msg_shutdown
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE md_run,                          ONLY: qs_mol_dyn
  USE mdctrl_types,                    ONLY: mdctrl_type
  USE message_passing,                 ONLY: mp_recv,&
                                             mp_send
  USE parallel_rng_types,              ONLY: reset_to_next_rng_substream
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE input_section_types,             ONLY: section_vals_val_set,section_vals_val_get
  
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
    TYPE(glbopt_config_type), INTENT(IN)     :: glbopt_cfg
    INTEGER                                  :: f_env_id
    INTEGER, INTENT(in)                      :: n_colvars
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(in)                      :: walker_group, my_walker_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'glbopt_drive_walker', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, stat, master_id, &
                                                msg, n_atoms, output_unit, tag
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: latest_positions
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cp_logger_type), POINTER            :: logger
    REAL(kind=dp), DIMENSION(n_colvars)      :: colvar_values
    REAL(KIND=dp), DIMENSION(glbopt_cfg%&
      epot_ndown+glbopt_cfg%epot_nup+1), &
      TARGET                                 :: epot_history
    TYPE(glbopt_mdctrl_data_type), TARGET    :: mdctrl_data
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mdctrl_type), POINTER               :: mdctrl_p
    TYPE(mdctrl_type), TARGET                :: mdctrl
    REAL(kind=dp)                            :: temp
    
! ====== Local Variables ======
!LOGICAL :: known_minima

    failure = .FALSE.
    
    ! ====== Usual Subroutine Setup ======
    CALL timeset(routineN,handle)
    
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_io_unit(logger)
    
                
    CALL f_env_add_defaults(f_env_id, f_env, error, failure)
    if(failure) stop "glbopt_drive_walker: f_env_add_defaults failed" 
    force_env => f_env%force_env
    CALL force_env_get(force_env, globenv=globenv, subsys=subsys, error=error)
    CALL cp_subsys_get(subsys, natom=n_atoms, error=error)
    ALLOCATE(latest_positions(1:3*n_atoms), stat=stat)
    if(stat /= 0) stop "glbopt_drive_walker: allocation failed"
     
    
    ! We want different random-number-streams for each walker
    DO i = 1, my_walker_id
       CALL reset_to_next_rng_substream(globenv%gaussian_rng_stream, error)
    END DO

    ! setup mdctrl callback
    mdctrl_data%cfg = glbopt_cfg
    mdctrl_data%epot_history => epot_history
    mdctrl%glbopt => mdctrl_data
    mdctrl_p => mdctrl
    
    WRITE(*,*) "setup done"
    DO WHILE(.TRUE.) !main-loop of minima hopping
      ! store starting configuration
      ! backup position as possible restart point 
      CALL pack_subsys_particles(subsys, r=latest_positions, error=error)
      
      ! reset md history
      mdctrl_data%epot_history = 0.0 
      mdctrl_data%md_mincounter = 0
      
      ! addjust temperture
      !motion_section => section_vals_get_subs_vals(force_env%root_section,"MOTION",error=error)
      !md_section => section_vals_get_subs_vals(motion_section,"MD",error=error)
      temp = 0.00013300624774681083
      CALL section_vals_val_set(force_env%root_section, "MOTION%MD%TEMPERATURE", r_val=temp, error=error)
      !CALL section_vals_val_get(force_env%root_section, "MOTION%MD%TEMPERATURE", r_val=temp, error=error)
      write (*,*) "got temp: ", temp
      CALL qs_mol_dyn(force_env, globenv, mdctrl=mdctrl_p, error=error)
      CALL cp_geo_opt(force_env, globenv, error=error)
      ! TODO: find out if minimazation was successfull
      ! TODO: find out if its a new minima
      
      
      ! calculate new colvar values and send them to the master
      DO i=1, n_colvars
         CALL colvar_eval_glob_f(i,force_env,error=error)
         colvar_values(i) = subsys%colvar_p(i)%colvar%ss
      END DO
      master_id = para_env%num_pe-1
      CALL mp_send(colvar_values,master_id,42,para_env%group)
      tag = 42
      CALL mp_recv(msg, master_id, tag, para_env%group)
      
      
      ! Follow the instruction just received from the master
      SELECT CASE (msg)
       CASE (glbopt_msg_minima_new)
          ! just continue with next loop cylce
       CASE (glbopt_msg_minima_old)
          ! reset to lastest starting positions
          CALL unpack_subsys_particles(subsys, r=latest_positions, error=error)
       CASE (glbopt_msg_shutdown)
          WRITE (*,*) "walker: received shutdown signal ", my_walker_id
                EXIT
       CASE DEFAULT
          STOP "Walker: received unkown message"
      END SELECT
      
      
    END DO
    
    
    ! ====== Clean up ======
    DEALLOCATE(latest_positions)
    CALL f_env_rm_defaults(f_env, error,stat)
    if(stat /= 0) stop "glbopt_drive_walker: f_env_rm_defaults failed"
    CALL destroy_force_env(f_env_id, stat)
    if(stat /= 0) stop "glbopt_drive_walker: destrey_forve_env failed"
    CALL timestop(handle)
  END SUBROUTINE glbopt_drive_walker
 
 
END MODULE glbopt_walker

