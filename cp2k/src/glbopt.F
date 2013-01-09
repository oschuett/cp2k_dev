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
                                             section_vals_val_set
  USE kinds,                           ONLY: dp,&
                                             default_path_length,&
                                             default_string_length

  USE timings,                         ONLY: timeset,&
                                             timestop
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_para_env,                     ONLY: cp_para_env_create
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
  USE glbopt_env_types,                ONLY: glbopt_env_type
  USE mdctrl_types,                    ONLY: mdctrl_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type,&
                                             cp_subsys_get,&
                                             pack_subsys_particles,&
                                             unpack_subsys_particles
  
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt'

 PUBLIC :: run_global_opt
  
 CONTAINS 
 
  
! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE run_global_opt(root_section, para_env, input_path, error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    TYPE(cp_error_type), INTENT(inout)       :: error
  
    ! ====== Local Variables ======
    CHARACTER(len=*), PARAMETER :: routineN = 'run_gobal_opt', &
       routineP = moduleN//':'//routineN
    LOGICAL                                  :: failure = .FALSE.
    
    INTEGER                                  :: handle, stat, step_start_val, ierr, walker_f_env_id
    
    
    
    
!    TYPE(gopt_f_type), POINTER               :: gopt_env
!    TYPE(gopt_param_type), POINTER           :: gopt_param
    INTEGER                                  :: n_walkers, n_masters, n_groups_created
    TYPE(cp_para_env_type), POINTER          :: walker_para_env
    TYPE(glbopt_env_type), POINTER           :: glbopt_env
        
    INTEGER :: i,  primus_slave, slave_group, slave_rank, num_slaves, pe_per_walker, walker_num
    INTEGER :: subgroup, subgroup_size, subgroup_rank, walker_group
    
     !  INTEGER :: dest, handle, i, i_job_to_restart, ierr, ijob, ijob_current, &
  ! ijob_end, ijob_start, iunit, n_jobs_to_run, new_group, new_output_unit, &
 !  new_rank, new_size, ngroups, num_slaves, output_unit, primus_slave, &
 ! , source, stat, tag, todo
    INTEGER, DIMENSION(:), POINTER           :: processors_partition, &
                                             group_distribution, &
                                             master_slave_partition, &
                                             slave_distribution, &
                                             slave_status
                                             
    LOGICAL                                  :: found, master, &
                                             run_OK, slave, im_a_master
    CHARACTER(len=default_path_length)       :: walker_output_path
    TYPE(f_env_type), POINTER                :: walker_f_env
    TYPE(force_env_type), POINTER            :: walker_force_env
    TYPE(cp_error_type)			      :: walker_error
    TYPE(global_environment_type), POINTER   :: walker_globenv
    CHARACTER(len=default_string_length)     :: project_name, walker_project_name
    
    
    ! ====== Usual Subroutine Setup ======
    CALL timeset(routineN,handle)
    
    
    WRITE(*,*) "!!! Juchu jetzt in eigenem module :-)"
    !TODO: create some kind of warning, that RUN_TYPE musst be NONE
    
    ! ====== Reading Input Section ======
    CALL section_vals_val_get(root_section,"GLOBAL_OPT%NUMBER_OF_WALKERS",i_val=n_walkers,error=error)
    
    
    
    !WRITE(unit_nr,'(T15,A)') "Success !!! now we could start global optimization :-)"
            
    n_masters = 1 ! fixed for now - might change in the future
    
        
    ! ====== Setup of MPI-Groups ======
    
    pe_per_walker = (para_env%num_pe-1)/n_walkers
    
    WRITE(*,*) "!!! Number of processors: ", para_env%num_pe
    WRITE(*,*) "!!! Number of masters: ", n_masters
    WRITE(*,*) "!!! Number of walkers: ", n_walkers
        
    IF (para_env%num_pe < n_masters + n_walkers)THEN
       WRITE(*,*) "There are not enough processes for n_masters+ n_walkers."
    ENDIF
    
    
    !IF ((para_env%num_pe > 1) .AND. (MOD(para_env%num_pe-n_masters, n_walkers) == 0)) THEN
    !  WRITE(*,FMT="(A,I0,A,A,I0,A,I0,A,I0,A)") "I have ",para_env%num_pe," processors and I will divide them ",&
    !    "into ",n_masters ," master and ",n_walkers," walkers with ",pe_per_walker," processors each."
    !ELSE
    !    WRITE(*,*) "number of processors-1 is not divisible by n_walkers."
    !ENDIF
    ! the full mpi group is first split in a slave group and a master group, the latter being at most 1 process
    
    
    im_a_master = (para_env%mepos < n_masters) ! task with lower ranks become masters
    
    ! First split split para_env%group into a master- and a walkers-groups...
    IF (im_a_master) THEN
    	CALL mp_comm_split_direct(para_env%group, subgroup, 1)
    	CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
    	WRITE (*,*) "Hello, this is a Master - there are ",subgroup_size, " of us."
    	CPPostcondition(subgroup_size==n_masters, cp_failure_level,routineP,error,failure)
    ELSE
    	CALL mp_comm_split_direct(para_env%group, subgroup, 2)
    	CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
    	WRITE (*,*) "Hello, this is a Walker - there are ",subgroup_size, " of us."
    	CPPostcondition(subgroup_size==para_env%num_pe-n_masters, cp_failure_level,routineP,error,failure)   	
    ENDIF
          
    
    IF (im_a_master) THEN
    	!CALL glbopt_run_master()
    ELSE
       ! ...then split walkers-group into n_walkers groups - one for each walker.
      	ALLOCATE(group_distribution(0:subgroup_size-1),STAT=stat)
        CALL mp_comm_split(subgroup, walker_group, n_groups_created, group_distribution, n_subgroups=n_walkers)
        walker_num = group_distribution(subgroup_rank)
        ! Did we get the desired number of walker-groups?
        CPPostcondition(n_groups_created==n_walkers, cp_failure_level,routineP,error,failure)
 	DEALLOCATE(group_distribution)
 	
        CALL run_walker(root_section, walker_group, input_path, walker_num, error)
         
        CALL mp_comm_free(walker_group)
    ENDIF
    
    CALL mp_comm_free(subgroup)
    CALL timestop(handle)

  END SUBROUTINE run_global_opt
 
  
  ! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE run_walker(root_section, mpi_comm, input_path, walker_num, error)
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: mpi_comm
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: walker_num
    TYPE(cp_error_type), INTENT(inout)       :: error
  
    ! ====== Local Variables ======
    CHARACTER(len=*), PARAMETER :: routineN = 'run_walker', &
       routineP = moduleN//':'//routineN
    LOGICAL                                  :: failure = .FALSE.
    
    INTEGER                                  :: handle, stat, output_unit, ierr, f_env_id, i, n_atoms
    
    TYPE(cp_logger_type), POINTER            :: logger
        
    CHARACTER(len=default_path_length)       :: output_path
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    CHARACTER(len=default_string_length)     :: project_name, new_project_name, walker_name
    TYPE(glbopt_env_type), POINTER           :: glbopt_env
    TYPE(mdctrl_type), POINTER               :: mdctrl
    TYPE(cp_subsys_type), POINTER            :: subsys
 
     
    ! ====== Usual Subroutine Setup ======
    CALL timeset(routineN,handle)
    
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_io_unit(logger)
    
    
    
    
    ! change the PROJECT_NAME to avoid clashes during backup of RESTART-files
    CALL section_vals_val_get(root_section, "GLOBAL%PROJECT_NAME", c_val=project_name,error=error)
    walker_name = "walker_"//ADJUSTL(cp_to_string(walker_num))
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
        	
    CALL f_env_add_defaults(f_env_id, f_env, error, failure)
    force_env => f_env%force_env
    CALL force_env_get(force_env, globenv=globenv, subsys=subsys, error=error)
    CALL cp_subsys_get(subsys, natom=n_atoms, error=error)
    
    
    ! ====== creating glbopt ======
    ALLOCATE(glbopt_env)
    glbopt_env%walker_num = walker_num
    glbopt_env%mincount_max = 3
    glbopt_env%epot_ndown = 2
    glbopt_env%epot_nup = 2
    ALLOCATE(glbopt_env%epot_history(glbopt_env%epot_ndown + glbopt_env%epot_nup + 1))
    ALLOCATE(glbopt_env%latest_positions(1:3*n_atoms))
    
    
    ! We want different random-number-streams for each walker
    DO i = 1, walker_num
       CALL reset_to_next_rng_substream(globenv%gaussian_rng_stream, error)
    END DO
        
    
    ALLOCATE(mdctrl)
    mdctrl%glbopt_env => glbopt_env
    
    DO i = 1, 3 !main-loop of minima hopping
    
      ! store starting configuration
      CALL pack_subsys_particles(subsys, r=glbopt_env%latest_positions, error=error)
      glbopt_env%epot_history = 0.0 
      glbopt_env%mincount = 0
      CALL qs_mol_dyn(force_env, globenv, mdctrl=mdctrl, error=error)
      CALL cp_geo_opt(force_env, globenv, error=error)
      ! find out if minimazation was successfull
      ! find out if its a new minima
      CALL unpack_subsys_particles(subsys, r=glbopt_env%latest_positions, error=error)
   END DO
    
    
    
    ! ====== Clean up ======
    CALL f_env_rm_defaults(f_env, error,ierr)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    CALL destroy_force_env(f_env_id,ierr)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)


  END SUBROUTINE run_walker
 
 
END MODULE glbopt

