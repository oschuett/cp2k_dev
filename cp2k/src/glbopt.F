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
  USE glbopt_env_types,                ONLY: glbopt_env_type
  USE mdctrl_types,                    ONLY: mdctrl_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type,&
                                             cp_subsys_get,&
                                             pack_subsys_particles,&
                                             unpack_subsys_particles
  USE colvar_methods,                  ONLY: colvar_eval_glob_f
  USE hash_functions,                  ONLY: b3hs_hash_key_jenkins
  
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt'

 PUBLIC :: run_global_opt
 
 
 INTEGER, PARAMETER  :: msg_minima_new = 111
 INTEGER, PARAMETER  :: msg_minima_old = 222
 INTEGER, PARAMETER  :: msg_shutdown   = 999
 
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
    CHARACTER(len=*), PARAMETER :: routineN = 'run_gobal_opt', &
       routineP = moduleN//':'//routineN
    LOGICAL                                  :: failure = .FALSE.
    
    INTEGER                                  :: handle, stat, ierr
    
    TYPE(glbopt_env_type), POINTER           :: glbopt_env
    
    INTEGER                                  :: n_groups_created
    TYPE(cp_para_env_type), POINTER          :: walker_para_env
    INTEGER, DIMENSION(:), POINTER           :: group_distribution               
        
    INTEGER :: i, pe_per_walker
    INTEGER :: subgroup, subgroup_size, subgroup_rank, walker_group
                                
    LOGICAL                                  :: im_a_master
    
    
    ! ====== Usual Subroutine Setup ======
    CALL timeset(routineN,handle)
    
    
    WRITE(*,*) "Global optimisation reporting for duty :-)"
    !TODO: create some kind of warning, that RUN_TYPE musst be NONE
    
    
    
    !WRITE(unit_nr,'(T15,A)') "Success !!! now we could start global optimization :-)"
            
    
    ! ====== Creating glbopt_env and reading input section =====
    ALLOCATE(glbopt_env)
    glbopt_env%n_masters = 1 ! fixed for now - might change in the future
    CALL section_vals_val_get(root_section,"GLOBAL_OPT%NUMBER_OF_WALKERS",i_val=glbopt_env%n_walkers,error=error)
    glbopt_env%md_mincount_max = 3
    glbopt_env%epot_ndown = 2
    glbopt_env%epot_nup = 2
    glbopt_env%hop_count_max = 100
    glbopt_env%full_para_env => para_env
    CALL cp_para_env_retain(para_env, error)
    
    !TODO: read from input or make educated guess
    glbopt_env%bloom_size = 2**8
    glbopt_env%n_hashes = 3
    
    
    ! make sure there are colvars
    colvar_section => section_vals_get_subs_vals(root_section,"FORCE_EVAL%SUBSYS%COLVAR",error=error)
    CALL section_vals_get(colvar_section,n_repetition=glbopt_env%n_colvars,error=error)
    IF(glbopt_env%n_colvars == 0) THEN
       WRITE (*,*) "No colvars found. Aborting."
       STOP  !TODO find official way of quiting
    END IF
    
        
    ! ====== Setup of MPI-Groups ======
    IF (para_env%num_pe < glbopt_env%n_masters + glbopt_env%n_walkers)THEN
       WRITE(*,*) "There are not enough processes for n_masters + n_walkers. Aborting."
       STOP !TODO find official way of quiting
    END IF
    
    pe_per_walker = (para_env%num_pe-1)/glbopt_env%n_walkers
    
    WRITE(*,*) "!!! Number of processors: ", para_env%num_pe
    WRITE(*,*) "!!! Number of masters: ", glbopt_env%n_masters
    WRITE(*,*) "!!! Number of walkers: ", glbopt_env%n_walkers
    WRITE(*,*) "!!! Number of colvars: ", glbopt_env%n_colvars
   
    
   
    
    
    
    !IF ((para_env%num_pe > 1) .AND. (MOD(para_env%num_pe-n_masters, n_walkers) == 0)) THEN
    !  WRITE(*,FMT="(A,I0,A,A,I0,A,I0,A,I0,A)") "I have ",para_env%num_pe," processors and I will divide them ",&
    !    "into ",n_masters ," master and ",n_walkers," walkers with ",pe_per_walker," processors each."
    !ELSE
    !    WRITE(*,*) "number of processors-1 is not divisible by n_walkers."
    !ENDIF
    ! the full mpi group is first split in a slave group and a master group, the latter being at most 1 process
    
    
    im_a_master = (para_env%mepos >= para_env%num_pe-glbopt_env%n_masters) ! task with highest ranks become masters
    
    ! First split split para_env%group into a master- and a walkers-groups...
    IF (im_a_master) THEN
    	CALL mp_comm_split_direct(para_env%group, subgroup, 1)
    	CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
    	WRITE (*,*) "Hello, this is a Master - there are ",subgroup_size, " of us."
    	CPPostcondition(subgroup_size==glbopt_env%n_masters, cp_failure_level,routineP,error,failure)
    ELSE
    	CALL mp_comm_split_direct(para_env%group, subgroup, 2)
    	CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
    	WRITE (*,*) "Hello, this is a Walker - there are ",subgroup_size, " of us."
    	CPPostcondition(subgroup_size==para_env%num_pe-glbopt_env%n_masters, cp_failure_level,routineP,error,failure)   	
    ENDIF

    
    IF (im_a_master) THEN
    	CALL run_master(glbopt_env, error)
    ELSE
        glbopt_env%my_master_id = para_env%num_pe-1 ! assign all walkers to same master 
        ! ...then split walkers-group into n_walkers groups - one for each walker.
      	ALLOCATE(group_distribution(0:subgroup_size-1),STAT=stat)
        CALL mp_comm_split(subgroup, walker_group, n_groups_created, group_distribution, n_subgroups=glbopt_env%n_walkers)
        glbopt_env%my_walker_id = group_distribution(subgroup_rank)
        ! Did we get the desired number of walker-groups?
        CPPostcondition(n_groups_created==glbopt_env%n_walkers, cp_failure_level,routineP,error,failure)
 	DEALLOCATE(group_distribution)
        CALL run_walker(glbopt_env, root_section, walker_group, input_path, error)
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
   SUBROUTINE run_master(glbopt_env, error)
      TYPE(glbopt_env_type), POINTER           :: glbopt_env
      TYPE(cp_error_type), INTENT(inout)       :: error
      
      ! ========= local vars ===========  
      REAL(kind=dp), dimension(glbopt_env%n_colvars + glbopt_env%n_hashes) :: key_parts
      INTEGER(kind=int_4), dimension(glbopt_env%n_hashes)  :: hashes                             
      CHARACTER, dimension(:), allocatable     :: key
      INTEGER(kind=int_4), dimension(glbopt_env%bloom_size/(4*8))       :: bloomfilter
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
      	CALL mp_recv(key_parts, calling_walker, tag, glbopt_env%full_para_env%group)
      	!WRITE(*,*) "Master: received new colvars from ", calling_walker
        !WRITE(*,*) "parts: ", key_parts
      	! following advise from:
      	! http://spyced.blogspot.ch/2009/01/all-you-ever-wanted-to-know-about.html
      	! Encode the colvar values as a character array and calculate hash
      	DO i=1, glbopt_env%n_hashes
          key = transfer(key_parts, key)    
          hashes(i) = b3hs_hash_key_jenkins(key, glbopt_env%bloom_size) - 1
          IF(hashes(i) > glbopt_env%bloom_size) STOP "hash value too large"
          IF(hashes(i) < 0 ) STOP "hash value too small"
          key_parts(glbopt_env%n_colvars + i) = hashes(i)
        END DO
        
        !write (*,*) "hahshes: ",hashes
        
        ! lockup and register minima at the same time 
        match = .TRUE.
        DO i=1, glbopt_env%n_hashes
          was_set = bitfield_set(hashes(i), bloomfilter)
          match = match .AND. was_set
        END DO
        
        if(match) then
          msg = msg_minima_old
        else
          msg = msg_minima_new
          minima_counter = minima_counter + 1 
        end if
        
        if(hop_counter > glbopt_env%hop_count_max) then
          msg = msg_shutdown
          shutdown_walkers = shutdown_walkers + 1 
        end if
        
        write (*,*) "minima_counter: ", minima_counter, "hop_counter: ",hop_counter
        
        CALL mp_send(msg, calling_walker, 42, glbopt_env%full_para_env%group)
        
        
        if(shutdown_walkers == glbopt_env%n_walkers) then
           write (*,*) "Master: shutdown all walkers - exiting main loop."
           exit
        end if
        
      END DO
      
      DEALLOCATE(key)
      
   END SUBROUTINE
   
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
  
! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE run_walker(glbopt_env, root_section, mpi_comm, input_path, error)
    TYPE(glbopt_env_type), POINTER           :: glbopt_env
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: mpi_comm
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    TYPE(cp_error_type), INTENT(inout)       :: error
  
    ! ====== Local Variables ======
    CHARACTER(len=*), PARAMETER :: routineN = 'run_walker', &
       routineP = moduleN//':'//routineN
    LOGICAL                                  :: failure = .FALSE.
    
    INTEGER                                  :: handle, stat, output_unit, ierr,&
                                                f_env_id, i, n_atoms
    
    TYPE(cp_logger_type), POINTER            :: logger
        
    CHARACTER(len=default_path_length)       :: output_path
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    CHARACTER(len=default_string_length)     :: project_name, new_project_name, walker_name
    
    TYPE(mdctrl_type), POINTER               :: mdctrl
    TYPE(cp_subsys_type), POINTER            :: subsys
    REAL(kind=dp), dimension(glbopt_env%n_colvars) :: colvar_values
    
    REAL(KIND=dp), DIMENSION(:), allocatable :: latest_positions
                                           
    INTEGER :: src, tag, msg
    !LOGICAL :: known_minima
    
    ! ====== Usual Subroutine Setup ======
    CALL timeset(routineN,handle)
    
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_io_unit(logger)
    
    
    ! change the PROJECT_NAME to avoid clashes during backup of RESTART-files
    CALL section_vals_val_get(root_section, "GLOBAL%PROJECT_NAME", c_val=project_name,error=error)
    walker_name = "walker_"//ADJUSTL(cp_to_string(glbopt_env%my_walker_id))
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
    
    WRITE(*,*) "allocating epot_history und co"
    CALL force_env_get(force_env, globenv=globenv, subsys=subsys, error=error)
    CALL cp_subsys_get(subsys, natom=n_atoms, error=error)
    ALLOCATE(latest_positions(1:3*n_atoms))
    ALLOCATE(glbopt_env%epot_history(glbopt_env%epot_ndown + glbopt_env%epot_nup + 1))
    
    
    WRITE(*,*) "setting up random number streams"
    
    ! We want different random-number-streams for each walker
    DO i = 1, glbopt_env%my_walker_id
       CALL reset_to_next_rng_substream(globenv%gaussian_rng_stream, error)
    END DO
    
    WRITE(*,*) "setting up md callback"
    
    ALLOCATE(mdctrl)
    mdctrl%glbopt_env => glbopt_env
    
    
    WRITE(*,*) "setup done"
    DO while(.TRUE.) !main-loop of minima hopping
      ! store starting configuration
      CALL pack_subsys_particles(subsys, r=latest_positions, error=error)
      glbopt_env%epot_history = 0.0 
      glbopt_env%md_mincounter = 0
      !WRITE(*,*) "starting md"
      CALL qs_mol_dyn(force_env, globenv, mdctrl=mdctrl, error=error)
      !WRITE(*,*) "starting geoopt"
      CALL cp_geo_opt(force_env, globenv, error=error)
      ! TODO: find out if minimazation was successfull
      ! TODO: find out if its a new minima
      
      !WRITE(*,*) "collecting colvars"
      colvar_values = 0.0_dp
      DO i=1, glbopt_env%n_colvars
         CALL colvar_eval_glob_f(i,force_env,error=error)
         !WRITE (*,*) "Colvar ", i, "value: ", subsys%colvar_p(i)%colvar%ss
         colvar_values(i) = subsys%colvar_p(i)%colvar%ss
      END DO
      !WRITE(*,*) "Walker: sending to master", glbopt_env%my_walker_id
      CALL mp_send(colvar_values,glbopt_env%my_master_id,42,glbopt_env%full_para_env%group)
      !WRITE(*,*) "Walker: waiting for answer", glbopt_env%my_walker_id
      src = glbopt_env%my_master_id
      tag = 42
      CALL mp_recv(msg, src, tag, glbopt_env%full_para_env%group)
      !WRITE(*,*) "Walker: got msg from master: ", msg, glbopt_env%my_walker_id
      
      
      SELECT CASE (msg)
       CASE (msg_minima_new)
          ! just continue with next loop cylce
       CASE (msg_minima_old)
          ! reset to lastest starting positions
          CALL unpack_subsys_particles(subsys, r=latest_positions, error=error)
       CASE (msg_shutdown)
          write (*,*) "walker: received shutdown signal ", glbopt_env%my_walker_id
      	  exit
       CASE DEFAULT
          STOP "Walker: received unkown message"
      END SELECT
      
      
    END DO
     
    
    
    ! ====== Clean up ======
    DEALLOCATE(glbopt_env%epot_history)
    DEALLOCATE(latest_positions)
    	    
    CALL f_env_rm_defaults(f_env, error,ierr)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    CALL destroy_force_env(f_env_id,ierr)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)


  END SUBROUTINE run_walker
 
 
END MODULE glbopt

