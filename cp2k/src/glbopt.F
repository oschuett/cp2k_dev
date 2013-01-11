!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt
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
  USE glbopt_types,                    ONLY: glbopt_config_type, glbopt_mdctrl_data_type
  
  USE glbopt_walker,                   ONLY: glbopt_drive_walker
  USE glbopt_master,                   ONLY: glbopt_drive_master
  
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt'

 PUBLIC :: run_global_opt, glbopt_config_type
 
  CONTAINS 
 
 

! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE run_global_opt(root_section, para_env, input_path, error)
    TYPE(section_vals_type), POINTER         :: root_section, colvar_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    TYPE(cp_error_type), INTENT(inout)       :: error
  
    ! ====== Local Variables ======
    TYPE(glbopt_config_type)                 :: glbopt_cfg
    
    !TODO: create some kind of warning, that RUN_TYPE musst be NONE
  

    CALL section_vals_val_get(root_section,"GLOBAL_OPT%NUMBER_OF_WALKERS",&
       i_val=glbopt_cfg%n_walkers,error=error)
    !TODO: read from input or make educated guess
    glbopt_cfg%n_hashes = 3
    glbopt_cfg%bloom_size = 2**8
    glbopt_cfg%epot_ndown = 2
    glbopt_cfg%epot_nup = 2
    glbopt_cfg%md_mincount_max = 3
    glbopt_cfg%hop_count_max = 100
 
    CALL run_global_opt_really(glbopt_cfg, root_section, para_env, input_path, error)    	    
    
   
  END SUBROUTINE run_global_opt  
  

   !===========================================================================
   !===========================================================================
   !=========================================================================== 
   SUBROUTINE run_global_opt_really(glbopt_cfg, root_section, para_env, input_path, error)
   TYPE(glbopt_config_type)                  :: glbopt_cfg
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    TYPE(cp_error_type), INTENT(inout)       :: error
  
    ! ====== Local Variables ======
    TYPE(section_vals_type), POINTER         :: colvar_section
    LOGICAL                                  :: failure
    
    INTEGER                                  :: handle, stat, ierr
    
    
    INTEGER                                  :: n_groups_created
    TYPE(cp_para_env_type), POINTER          :: walker_para_env
    INTEGER, DIMENSION(0:para_env%num_pe-2), TARGET  :: group_distribution
    INTEGER, DIMENSION(:), POINTER          :: group_distribution_p
        
    INTEGER :: i, pe_per_walker, n_colvars, f_env_id, my_walker_id
    INTEGER :: subgroup, subgroup_size, subgroup_rank, walker_group
                                
    LOGICAL                                  :: im_the_master
    
    
    ! setup timings
    CHARACTER(len=*), PARAMETER :: routineN = 'run_gobal_opt_really', &
       routineP = moduleN//':'//routineN
    CALL timeset(routineN,handle)
    failure = .FALSE.
    
    
    WRITE(*,*) "Global optimisation reporting for duty :-)"
    
    
    
    ! make sure there are colvars
    colvar_section => section_vals_get_subs_vals(root_section,"FORCE_EVAL%SUBSYS%COLVAR",error=error)
    CALL section_vals_get(colvar_section,n_repetition=n_colvars,error=error)
    IF(n_colvars == 0) THEN
       WRITE (*,*) "No colvars found. Aborting."
       STOP  !TODO find official way of quiting
    END IF
    
        
    ! ====== Setup of MPI-Groups ======
    IF (para_env%num_pe < glbopt_cfg%n_walkers + 1)THEN
       WRITE(*,*) "There are not enough processes for n_walkers + 1. Aborting."
       STOP !TODO find official way of quiting
    END IF
    
    pe_per_walker = (para_env%num_pe-1)/glbopt_cfg%n_walkers
    
    WRITE(*,*) "!!! Number of processors: ",para_env%num_pe
    WRITE(*,*) "!!! Number of walkers: ", glbopt_cfg%n_walkers
    WRITE(*,*) "!!! Number of colvars: ", n_colvars
    
    !IF ((para_env%num_pe > 1) .AND. (MOD(para_env%num_pe-n_masters, n_walkers) == 0)) THEN
    !  WRITE(*,FMT="(A,I0,A,A,I0,A,I0,A,I0,A)") "I have ",para_env%num_pe," processors and I will divide them ",&
    !    "into ",n_masters ," master and ",n_walkers," walkers with ",pe_per_walker," processors each."
    !ELSE
    !    WRITE(*,*) "number of processors-1 is not divisible by n_walkers."
    !ENDIF
    ! the full mpi group is first split in a slave group and a master group, the latter being at most 1 process
    
    
    ! the last task becomes the master. Keeps node-alignment of other tasks.
    im_the_master = (para_env%mepos == para_env%num_pe-1)
    
    ! First split split para_env%group into a master- and a walkers-groups...
    IF (im_the_master) THEN
    	CALL mp_comm_split_direct(para_env%group, subgroup, 1)
    	CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
    	WRITE (*,*) "Hello, this is a Master - there are ",subgroup_size, " of us."
    	CPPostcondition(subgroup_size==1, cp_failure_level,routineP,error,failure)
    ELSE
    	CALL mp_comm_split_direct(para_env%group, subgroup, 2)
    	CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
    	WRITE (*,*) "Hello, this is a Walker - there are ",subgroup_size, " of us."
    	CPPostcondition(subgroup_size==para_env%num_pe-1, cp_failure_level,routineP,error,failure)   	
    ENDIF

    
    IF (im_the_master) THEN
    	CALL glbopt_drive_master(glbopt_cfg, para_env, n_colvars, error)
    ELSE
        ! ...then split walkers-group into n_walkers groups - one for each walker.
      	!ALLOCATE(group_distribution(0:subgroup_size-1),STAT=stat)
      	group_distribution_p => group_distribution
        CALL mp_comm_split(subgroup, walker_group, n_groups_created, group_distribution_p, n_subgroups=glbopt_cfg%n_walkers)
        my_walker_id = group_distribution(subgroup_rank)
        ! Did we get the desired number of walker-groups?
        CPPostcondition(n_groups_created==glbopt_cfg%n_walkers, cp_failure_level,routineP,error,failure)
 	!DEALLOCATE(group_distribution)
 	f_env_id = walker_create_f_env(root_section, input_path, walker_group, my_walker_id, error)
        CALL glbopt_drive_walker(glbopt_cfg, f_env_id, n_colvars, para_env, walker_group, my_walker_id, error)
        CALL mp_comm_free(walker_group)
    ENDIF
    
    CALL mp_comm_free(subgroup)
    CALL timestop(handle)

  END SUBROUTINE run_global_opt_really
 

  ! *****************************************************************************
  
  FUNCTION walker_create_f_env(root_section, input_path, mpi_comm, walker_id, error) RESULT(f_env_id)
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: mpi_comm
    INTEGER, INTENT(IN)                      :: walker_id
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: f_env_id
    
    
    !local vars
    CHARACTER(len=default_string_length)     :: project_name, new_project_name, walker_name
    LOGICAL                                  :: failure
    INTEGER                                  :: ierr
    CHARACTER(len=default_path_length)       :: output_path
    CHARACTER(len=*), PARAMETER :: routineN = 'walker_create_f_env', &
       routineP = moduleN//':'//routineN
    failure = .FALSE.
    
    ! change the PROJECT_NAME to avoid clashes during backup of RESTART-files
    CALL section_vals_val_get(root_section, "GLOBAL%PROJECT_NAME", c_val=project_name,error=error)
    walker_name = "walker_"//ADJUSTL(cp_to_string(walker_id))
    new_project_name = TRIM(project_name)//"-"//walker_name
    CALL section_vals_val_set(root_section, "GLOBAL%PROJECT_NAME", c_val=new_project_name,error=error)
    
    output_path = TRIM(new_project_name)//".out"
    WRITE (*,*) output_path
    
    ! ====== Setting up the force_env ======
    CALL create_force_env(f_env_id, &
                          input_path=input_path, &
                          input=root_section, &
                          output_path=output_path, &
                          mpi_comm=mpi_comm,&
                          ierr=ierr)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    
  END FUNCTION walker_create_f_env
  
  
END MODULE glbopt

