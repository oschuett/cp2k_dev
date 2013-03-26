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
                                             glbopt_cmd_sync,&
                                             glbopt_cmd_report,&
                                             glbopt_cmd_shutdown,&
                                             glbopt_command_type,&
                                             glbopt_report_type,&
                                             glbopt_walker_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp, sp, int_4, int_8
  USE md_run,                          ONLY: qs_mol_dyn
  USE mdctrl_types,                    ONLY: mdctrl_type
  USE message_passing,                 ONLY: mp_recv,&
                                             mp_send
  USE parallel_rng_types,              ONLY: reset_to_next_rng_substream
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE input_section_types,             ONLY: section_vals_val_set,section_vals_val_get
  USE hash_functions,                  ONLY: b3hs_hash_key_jenkins
  USE physcon,                         ONLY: kelvin, femtoseconds
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length
  USE f77_interface,                   ONLY: create_force_env
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE particle_list_types,             ONLY: particle_list_type
  USE kinds,                           ONLY: dp
  USE glbopt_utils,                    ONLY: particles_connected

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_walker'

 PUBLIC :: glbopt_walker_init, glbopt_walker_finalize
 PUBLIC :: glbopt_walker_execute_command, glbopt_walker_assemble_report


 CONTAINS




! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_walker_init(walker, glbopt_cfg, root_section, input_path,&
                output_unit, para_env, walker_group, walker_id, error)
       TYPE(glbopt_walker_type), INTENT(INOUT)         :: walker
       TYPE(glbopt_config_type), INTENT(IN)     :: glbopt_cfg
       TYPE(section_vals_type), POINTER         :: root_section
       CHARACTER(LEN=*), INTENT(IN)             :: input_path
       INTEGER, INTENT(IN)                      :: output_unit
       TYPE(cp_para_env_type), POINTER          :: para_env
       INTEGER, INTENT(in)                      :: walker_group, walker_id
       TYPE(cp_error_type), INTENT(inout)       :: error

       TYPE(global_environment_type), POINTER   :: globenv
       CHARACTER(len=default_path_length)       :: output_path
       CHARACTER(len=default_string_length)     :: new_project_name, &
                                                   project_name, walker_name
       TYPE(force_env_type), POINTER            :: force_env
       TYPE(f_env_type), POINTER                :: f_env
       TYPE(cp_subsys_type), POINTER            :: subsys
       TYPE(cp_logger_type), POINTER            :: logger
       INTEGER                                  :: f_env_id
       INTEGER                                  :: i
       LOGICAL                                  :: failure

       failure = .FALSE.

        ! change the PROJECT_NAME to avoid clashes during backup of RESTART-files
       output_path = "__STD_OUT__"
       IF(glbopt_cfg%n_walkers > 1) THEN
          IF(walker_id > 999) STOP "walker_create_f_env: Did not expect so many walkers."
          CALL section_vals_val_get(root_section, "GLOBAL%PROJECT_NAME", c_val=project_name,error=error)
          WRITE(walker_name,"('WALKER',I3.3)"), walker_id
          new_project_name = TRIM(project_name)//"-"//walker_name
          CALL section_vals_val_set(root_section, "GLOBAL%PROJECT_NAME", c_val=new_project_name,error=error)
          output_path = TRIM(new_project_name)//".out"
       END IF

      ! ======= Create f_env =======
      CALL create_force_env(f_env_id, &
                            input_path=input_path, &
                            input=root_section, &
                            output_path=output_path, &
                            output_unit=output_unit,&
                            mpi_comm=walker_group)


      ! ======= More setup stuff =======
      CALL f_env_add_defaults(f_env_id, f_env, error, failure)
      if(failure) stop "glbopt_drive_walker: f_env_add_defaults failed" 
      force_env => f_env%force_env
      CALL force_env_get(force_env, globenv=globenv, subsys=subsys, error=error)
      !CALL cp_subsys_get(subsys, natom=n_atoms, error=error)
      !ALLOCATE(prev_positions(1:3*n_atoms))
      !prev_positions = 0.0_dp

      ! We want different random-number-streams for each walker
      DO i = 1, walker_id
         CALL reset_to_next_rng_substream(globenv%gaussian_rng_stream, error)
      END DO

      ! getting an output unit for logging
      logger => cp_error_get_logger(error)
      walker%iw = cp_print_key_unit_nr(logger,root_section,&
          "GLOBAL_OPT%PRINT%WALKER_RUN_INFO",extension=".walkerLog",error=error)

!    ! setup mdctrl callback
!    walker_env%mdctrl_data%cfg = glbopt_cfg
!    mdctrl_data%epot_history => epot_history
!    mdctrl_data%output_unit = iw
!    mdctrl%glbopt => mdctrl_data
!    mdctrl_p => mdctrl
  END SUBROUTINE glbopt_walker_init



! *****************************************************************************
! *****************************************************************************
   SUBROUTINE print(walker, message)
      TYPE(glbopt_walker_type), INTENT(INOUT)         :: walker
      CHARACTER(LEN=*)                         :: message
      IF (walker%iw>0) &
        WRITE (walker%iw,*) "GLOBAL_OPT| "//TRIM(message)
   END SUBROUTINE print


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_walker_execute_command(walker, cmd, should_stop)
       TYPE(glbopt_walker_type), INTENT(INOUT)         :: walker
       TYPE(glbopt_command_type), INTENT(IN)           :: cmd
       LOGICAL, INTENT(INOUT)                   :: should_stop
       SELECT CASE (cmd%id)
         CASE (glbopt_cmd_shutdown)
           call print(walker, "Received shutdown command, quitting.")
           should_stop = .TRUE.
         CASE DEFAULT
           STOP "Walker: received unkown command"
       END SELECT
   END SUBROUTINE glbopt_walker_execute_command

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_walker_assemble_report(walker, report)
       TYPE(glbopt_walker_type), INTENT(INOUT)     :: walker
       TYPE(glbopt_report_type), INTENT(INOUT)     :: report
   END SUBROUTINE glbopt_walker_assemble_report

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_walker_finalize(walker)
       TYPE(glbopt_walker_type), INTENT(INOUT)     :: walker
   END SUBROUTINE glbopt_walker_finalize

!! *****************************************************************************
!!> \brief Main driver to perform global optimization
!!> \param error variable to control error logging, stopping,...
!!>        see module cp_error_handling
!! *****************************************************************************
!   SUBROUTINE glbopt_drive_walker(glbopt_cfg, root_section, input_path, output_unit, para_env, walker_group, walker_id, error)
!    TYPE(glbopt_config_type), INTENT(IN)     :: glbopt_cfg
!    TYPE(section_vals_type), POINTER         :: root_section
!    CHARACTER(LEN=*), INTENT(IN)             :: input_path
!    INTEGER, INTENT(IN)                      :: output_unit
!    TYPE(cp_para_env_type), POINTER          :: para_env
!    INTEGER, INTENT(in)                      :: walker_group, walker_id
!    TYPE(cp_error_type), INTENT(inout)       :: error
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'glbopt_drive_walker', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: handle, i, stat, master_id, &
!                                                msg, n_atoms, tag,&
!                                                i_mainloop
!    LOGICAL                                  :: failure
!    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: prev_positions
!    TYPE(force_env_type), POINTER            :: force_env
!    TYPE(f_env_type), POINTER                :: f_env
!    TYPE(cp_subsys_type), POINTER            :: subsys
!    TYPE(cp_logger_type), POINTER            :: logger
!    INTEGER                                  :: f_env_id
!
!
!    REAL(KIND=dp), DIMENSION(glbopt_cfg%&
!      epot_ndown+glbopt_cfg%epot_nup+1), &
!      TARGET                                 :: epot_history
!    TYPE(glbopt_mdctrl_data_type), TARGET    :: mdctrl_data
!    TYPE(global_environment_type), POINTER   :: globenv
!    TYPE(mdctrl_type), POINTER               :: mdctrl_p
!    TYPE(mdctrl_type), TARGET                :: mdctrl
!    REAL(kind=dp)                            :: temperature, prev_epot, Eaccept, epot
!    REAL(kind=dp)                            :: Elowest
!    REAL(kind=dp)                            :: timestep
!
!    INTEGER(kind=int_4), DIMENSION(glbopt_cfg%bloom_size/(4*8)) :: bitfield
!    TYPE(INTEGER(kind=int_4)), DIMENSION(glbopt_cfg%n_hashes) :: bloomidx
!    LOGICAL                                  :: minima_known
!    INTEGER                                  :: iw
!    CHARACTER(len=default_path_length)       :: output_path
!    CHARACTER(len=default_string_length)     :: new_project_name, &
!                                                project_name, walker_name
!    INTEGER                                  :: prev_iframe, iframe
!    INTEGER                                  :: gopt_max_iter
!    INTEGER                                  :: n_accepted
!    failure = .FALSE.
!    TYPE(walker_env_type)                    :: walker_env
!    ! ====== Usual Subroutine Setup ======
!    CALL timeset(routineN,handle)
!   
!    CALL init(walker_env, glbopt_cfg, root_section, input_path, output_unit, para_env, walker_group, walker_id, error)
!
!    !output_unit = cp_logger_get_default_io_unit(logger)
!
!    mainloop: DO i_mainloop = 1, HUGE(1)  !infinite loop
!       CALL request_order(walker_env, order)
!       CALL execute_order(walker_env, order, should_stop)
!       IF(should_stop) EXIT
!       CALL report_status(walker_env)
!    END DO mainloop
!
!    CALL finalize(walker_env)
!
!    Eaccept = glbopt_cfg%Eaccept0
!    prev_epot = Huge(1.0_dp)
!    epot = Huge(1.0_dp)
!    Elowest = Huge(1.0_dp)
!    n_accepted = 0
!    ! store starting configuration,  backup position as possible restart point
!    CALL pack_subsys_particles(subsys, r=prev_positions, error=error)
!    CALL section_vals_val_get(root_section, "MOTION%MD%TEMPERATURE", r_val=temperature, error=error)
!    ! fetch original value from input
!    CALL section_vals_val_get(root_section, "MOTION%GEO_OPT%MAX_ITER", i_val=gopt_max_iter, error=error)
!
!    ! -----------main-loop of minima hopping---------------
!    iframe = 0 !indicates the NEXT frame to be written
!    mainloop: DO i_mainloop = 1, HUGE(1)  !infinite loop
!      ! quit when known minima is reached
!      IF(glbopt_cfg%Emin > epot) THEN
!         IF (iw>0)  WRITE (iw,*), 'GLOBAL_OPT| Reached Emin > Epot. Quitting.'
!         EXIT mainloop
!      END IF
!
!      timestep = SQRT(MIN(1.0, 300.0/(temperature*kelvin))) / femtoseconds
!      IF (iw>0) THEN
!        WRITE (iw,'(A,A,A,I9,A,A)') ' GLOBAL_OPT| ',REPEAT("*", 15),' ',&
!          i_mainloop,'. Minima Hopping Iteration ',REPEAT("*", 15)
!        WRITE (iw,'(A,11X,E20.10)') ' GLOBAL_OPT| Lowest minima found so far [Hartree]',Elowest
!        WRITE (iw,'(A,32X,I10)')    ' GLOBAL_OPT| Number of accepted minima', n_accepted
!        WRITE (iw,'(A,20X,E20.10)') ' GLOBAL_OPT| Acceptance Energy [Hartree]',Eaccept
!        WRITE (iw,'(A,29X,F20.3)')  ' GLOBAL_OPT| MD temperature [K]',temperature*kelvin
!        WRITE (iw,'(A,31X,F20.3)')  ' GLOBAL_OPT| MD timestep [fs]',timestep*femtoseconds
!      END IF
!
!      ! quit when system explodes
!      !IF(temperature*kelvin > 1000) THEN
!      !   IF (iw>0)  WRITE (iw,*) " GLOBAL_OPT| Reached T>1000K, quiting minima hopping."
!      !   EXIT mainloop
!      !END IF
!
!      IF (iw>0) WRITE (iw,'(A,25X,I10)') " GLOBAL_OPT| Starting MD at trajectory frame ", iframe
!      mdctrl_data%epot_history = 0.0
!      mdctrl_data%md_mincounter = 0
!      prev_iframe = iframe
!      IF(iframe==0) iframe = 1 ! qs_mol_dyn behaves differently for STEP_START_VAL=0
!      CALL section_vals_val_set(root_section, "MOTION%MD%STEP_START_VAL", i_val=iframe-1, error=error)
!      CALL section_vals_val_set(root_section, "MOTION%MD%TEMPERATURE", r_val=temperature, error=error)
!      CALL section_vals_val_set(root_section, "MOTION%MD%TIMESTEP", r_val=timestep, error=error)
!      CALL qs_mol_dyn(force_env, globenv, mdctrl=mdctrl_p, error=error)
!      iframe = mdctrl_data%itimes + 1
!      IF (iw>0) WRITE (iw,'(A,I4,A)') " GLOBAL_OPT| md ended after ", iframe-prev_iframe, " steps."
!
!      IF (iw>0) WRITE (iw,'(A,9X,I10)') " GLOBAL_OPT| Starting local optimisation at trajectory frame ", iframe
!      CALL section_vals_val_set(root_section, "MOTION%GEO_OPT%STEP_START_VAL", i_val=iframe-1, error=error)
!      CALL section_vals_val_set(root_section, "MOTION%GEO_OPT%MAX_ITER", i_val=iframe+gopt_max_iter, error=error)
!      CALL cp_geo_opt(force_env, globenv, rm_restart_info=.FALSE., error=error)
!      prev_iframe = iframe
!      CALL section_vals_val_get(root_section, "MOTION%GEO_OPT%STEP_START_VAL", i_val=iframe, error=error)
!      iframe = iframe + 2 ! Compensates for different START_VAL interpretation.
!      IF (iw>0) WRITE (iw,'(A,I4,A)') " GLOBAL_OPT| gopt ended after ", iframe-prev_iframe-1, " steps."
!      ! Stop if cluster starts to splits into multiple fragments.
!      IF(.NOT. particles_connected(force_env, glbopt_cfg%gopt_max_dist, error)) THEN
!         IF(iw > 0) WRITE (iw,"(A)") " GLOBAL_OPT| Conclusion: rejecting minima because of GOPT_MAX_DIST."
!         temperature = temperature*glbopt_cfg%beta3 !decreasing temperature
!         CALL unpack_subsys_particles(subsys, r=prev_positions, error=error)
!         CYCLE mainloop
!      END IF
!
!      CALL force_env_get(force_env, potential_energy=epot, error=error)
!      IF (iw>0) WRITE (iw,'(A,21X,E20.10)')' GLOBAL_OPT| Potential Energy [Hartree]',Epot
!
!      CALL calc_bloom_indices(glbopt_cfg, force_env, subsys, bloomidx, error)
!      !IF (iw>0) WRITE (iw,*) "GLOBAL_OPT| got bloom indices: ", bloomidx
!
!      minima_known = bloomfilter_query(bitfield, bloomidx)
!
!      IF(.NOT. minima_known .AND. glbopt_cfg%n_walkers > 1) THEN
!         ! Sync with master and check again.
!         CALL master_sync_bloomfilter(para_env, bitfield)
!         minima_known = bloomfilter_query(bitfield, bloomidx)
!      END IF
!
!      IF(minima_known) THEN
!          IF (iw>0) WRITE (iw,*) "GLOBAL_OPT| Conclusion: old minima"
!          ! known minima restart MD
!          temperature = temperature*glbopt_cfg%beta1
!          CALL unpack_subsys_particles(subsys, r=prev_positions, error=error)
!          CYCLE mainloop
!      END IF
!
!      ! found a new minima
!
!      temperature = temperature*glbopt_cfg%beta3
!
!      IF (epot - prev_epot > Eaccept) THEN
!         ! rejected new minimum due to energy.
!         IF (iw>0) WRITE (iw,*) "GLOBAL_OPT| Conclusion: new minima, but got rejected due to Eaccept"
!         CALL unpack_subsys_particles(subsys, r=prev_positions, error=error)
!         Eaccept = Eaccept * glbopt_cfg%alpha2
!         CYCLE mainloop
!      END IF
!
!     ! Hurra, we accept a new minimum
!     IF (iw>0) WRITE (iw,*) "GLOBAL_OPT| Conclusion: new minima got accepted"
!     CALL bloomfilter_add(bitfield, bloomidx)
!     Eaccept = Eaccept * glbopt_cfg%alpha1
!     prev_epot = epot
!     Elowest = MIN(Elowest, epot)
!     n_accepted = n_accepted + 1
!     CALL pack_subsys_particles(subsys, r=prev_positions, error=error)
!     IF(glbopt_cfg%n_walkers > 1) &
!        CALL master_report_minima(para_env, bitfield, Epot)
!
!! ========== OLD CRAP ==========================================================
!!      ! ok the mimima would be acceptable energy-wise, let's see if its new
!!      ! calculate new colvar values and send them to the master
!!      DO i=1, n_colvars
!!         CALL colvar_eval_glob_f(i,force_env,error=error)
!!         colvar_values(i) = subsys%colvar_p(i)%colvar%ss
!!      END DO
!!      master_id = para_env%num_pe-1
!!      CALL mp_send(colvar_values,master_id,42,para_env%group)
!!      tag = 42
!!      CALL mp_recv(msg, master_id, tag, para_env%group)
!!
!
!      ! Follow the instruction just received from the master
!!      SELECT CASE (msg)
!!       CASE (glbopt_msg_minima_new)
!!           ! found a new minimum - and accepted it
!!           Eaccept = Eaccept * glbopt_cfg%alpha1
!!       CASE (glbopt_msg_minima_old)
!!          ! Found previously found minimum - reset to lastest starting positions
!!          temperature = temperature * glbopt_cfg%beta2
!!          CALL unpack_subsys_particles(subsys, r=prev_positions, error=error)
!!
!!       CASE (glbopt_msg_shutdown)
!!          WRITE (*,*) "walker: received shutdown signal ", walker_id
!!                EXIT
!!       CASE DEFAULT
!!          STOP "Walker: received unkown message"
!!      END SELECT
!!
!
!!          if not unique:
!!170	            self._temperature *= self._beta2
!!171	            self._log('msg', 'Found previously found minimum.')
!!172	            self._log('par')
!!173	            if self._previous_optimum:
!!174	                self._log('msg', 'Restoring last minimum.')
!!175	                self._atoms.positions = self._previous_optimum.positions
!!176	            return
!!177	        # Must have found a unique minimum.
!!178	        self._temperature *= self._beta3
!!179	        self._log('msg', 'Found a new minimum.')
!!180	        self._log('par')
!!181	        if (self._atoms.get_potential_energy() <
!!182	            self._previous_energy + self._Eaccept):
!!183	                self._log('msg', 'Accepted new minimum.')
!!184	                self._Eaccept *= self._alpha1
!!185	                self._log('par')
!!186	                self._record_minimum()
!!187	        else:
!!188	            self._log('msg', 'Rejected new minimum due to energy. '
!!189	                             'Restoring last minimum.')
!!190	            self._atoms.positions = self._previous_optimum.positions
!!191	            self._Eaccept *= self._alpha2
!!192	            self._log('par')
!
!
!
!    END DO mainloop
!
!
!    ! ====== Clean up ======
!    DEALLOCATE(prev_positions)
!    CALL f_env_rm_defaults(f_env, error,stat)
!    if(stat /= 0) stop "glbopt_drive_walker: f_env_rm_defaults failed"
!    CALL destroy_force_env(f_env_id, stat)
!    if(stat /= 0) stop "glbopt_drive_walker: destrey_forve_env failed"
!    CALL timestop(handle)
!  END SUBROUTINE glbopt_drive_walker


!
!! *****************************************************************************
!!> \brief
!! *****************************************************************************
!  SUBROUTINE master_report_minima(para_env, bitfield, Epot)
!     TYPE(cp_para_env_type), POINTER          :: para_env
!     INTEGER(kind=int_4), DIMENSION(:), INTENT(in)    :: bitfield
!     REAL(kind=dp)                            :: Epot
!     INTEGER :: master_id, tag
!
!     master_id = para_env%num_pe-1
!     tag = glbopt_msg_report
!     CALL mp_send(bitfield, master_id, tag, para_env%group)
!     CALL mp_send(Epot,     master_id, tag, para_env%group)
!     !CALL mp_recv(bitfield, master_id, tag, para_env%group)
!  END SUBROUTINE master_report_minima
!
!  
!! *****************************************************************************
!!> \brief
!! *****************************************************************************
!  SUBROUTINE master_sync_bloomfilter(para_env, bitfield)
!     TYPE(cp_para_env_type), POINTER          :: para_env
!     INTEGER(kind=int_4), DIMENSION(:), INTENT(inout)    :: bitfield
!     INTEGER :: master_id, tag
!
!     master_id = para_env%num_pe-1
!     tag = glbopt_msg_sync
!     CALL mp_send(bitfield, master_id, tag, para_env%group)
!     CALL mp_recv(bitfield, master_id, tag, para_env%group)
!  END SUBROUTINE master_sync_bloomfilter
!
!! *****************************************************************************
!!> \brief
!! *****************************************************************************
!  SUBROUTINE bloomfilter_add(bitfield, bloomidx)
!     INTEGER(kind=int_4), DIMENSION(:), INTENT(inout)    :: bitfield
!     TYPE(INTEGER(kind=int_4)), DIMENSION(:), INTENT(IN) :: bloomidx
!     INTEGER                 :: i
!     DO i=1, SIZE(bloomidx)
!        CALL bitfield_set(bloomidx(i), bitfield)
!     END DO
!     !TODO count occupancy, warn if too high
!  END SUBROUTINE bloomfilter_add
!
!
!! *****************************************************************************
!!> \brief
!! *****************************************************************************
!  FUNCTION bloomfilter_query(bitfield, bloomidx) RESULT(is_member)
!     INTEGER(kind=int_4), DIMENSION(:), INTENT(inout)    :: bitfield
!     TYPE(INTEGER(kind=int_4)), DIMENSION(:), INTENT(IN) :: bloomidx
!     LOGICAL :: is_member
!     INTEGER                 :: i
!     is_member = .TRUE.
!     DO i=1, SIZE(bloomidx)
!        is_member = is_member .AND. bitfield_get(bloomidx(i), bitfield)
!     END DO
!  END FUNCTION bloomfilter_query
!
!
!! *****************************************************************************
!!> \brief
!! *****************************************************************************
!  SUBROUTINE calc_bloom_indices(glbopt_cfg, force_env, subsys, bloomidx, error)
!     TYPE(glbopt_config_type), INTENT(IN)                        :: glbopt_cfg
!     TYPE(force_env_type), POINTER                               :: force_env
!     TYPE(cp_subsys_type), POINTER                               :: subsys
!     TYPE(INTEGER(kind=int_4)), DIMENSION(glbopt_cfg%n_hashes), INTENT(OUT) :: bloomidx
!     TYPE(cp_error_type), INTENT(inout)                          :: error
!
!
!     REAL(kind=dp), DIMENSION(glbopt_cfg%n_colvars+glbopt_cfg%n_hashes) :: key_parts
!     CHARACTER, ALLOCATABLE, DIMENSION(:)                        :: key
!     INTEGER                                                     :: key_length, i, handle
!     INTEGER(kind=int_4), DIMENSION(glbopt_cfg%n_hashes)         :: hashes
!
!
!     CHARACTER(len=*), PARAMETER :: routineN = 'calc_bloom_indices', &
!      routineP = moduleN//':'//routineN
!     CALL timeset(routineN,handle)
!
!     key_parts = 0.0_dp
!     DO i=1, glbopt_cfg%n_colvars
!        CALL colvar_eval_glob_f(i,force_env,error=error)
!        key_parts(i) = subsys%colvar_p(i)%colvar%ss
!     END DO
!
!     !WRITE (*,*) "colvars: ", key_parts(1:glbopt_cfg%n_colvars)
!     key_parts = REAL(FLOOR(key_parts/glbopt_cfg%colvar_precision, int_8), dp)*glbopt_cfg%colvar_precision
!     !WRITE (*,*) "colvars: ", key_parts(1:glbopt_cfg%n_colvars)
!
!     !WRITE (*,*) "colvars: ", key_parts(1:glbopt_cfg%n_colvars)
!     !WRITE (*,*) "blub: ", key_parts(1)*10**8, REAL(FLOOR(key_parts(1)*10**8, int_8), dp)/10**8, REAL(key_parts(1), sp)
!     !key_parts = FLOOR(key_parts*10**8)/10**8
!     !WRITE (*,*) "colvars: ", key_parts(1:glbopt_cfg%n_colvars)
!
!     ! following advise from:
!     ! http://spyced.blogspot.ch/2009/01/all-you-ever-wanted-to-know-about.html
!     ! Encode the colvar values as a character array and calculate hash
!     key_length = SIZE(TRANSFER(key_parts, key))
!     ALLOCATE(key(key_length))
!     DO i=1, glbopt_cfg%n_hashes
!       key = TRANSFER(key_parts, key)
!       bloomidx(i) = b3hs_hash_key_jenkins(key, glbopt_cfg%bloom_size) - 1
!       IF(bloomidx(i) > glbopt_cfg%bloom_size) STOP "hash value too large"
!       IF(bloomidx(i) < 0 ) STOP "hash value too small"
!       key_parts(glbopt_cfg%n_colvars + i) = bloomidx(i)
!     END DO
!
!     CALL timestop(handle)
!  END SUBROUTINE
!
!! *****************************************************************************
!!> \brief
!! *****************************************************************************
!  SUBROUTINE bitfield_set(idx, bitfield)
!    INTEGER(kind=int_4), INTENT(in)          :: idx
!    INTEGER(kind=int_4), DIMENSION(:), &
!      INTENT(inout)                          :: bitfield
!
!    INTEGER                                  :: mask, offset, val
!
!! Caution: This code only works with 4 byte integers
!
!     offset = ISHFT(idx, -5) + 1
!     val = IAND(idx, 31)
!     mask = 2**val
!
!     IF(offset > SIZE(bitfield)) THEN
!       WRITE (*,*) "idx: ", idx, "offset: ",offset, "size(bitfield): ",SIZE(bitfield)
!       STOP "set_bitfield: overflow"
!     END IF
!
!     IF(offset < 1) THEN
!       WRITE (*,*) "idx: ", idx, "offset: ",offset, "size(bitfield): ",SIZE(bitfield)
!       STOP "set_bitfield: underflow"
!     END IF
!
!     bitfield(offset) = IOR(bitfield(offset), mask)
!  END SUBROUTINE bitfield_set
!
!! *****************************************************************************
!!> \brief
!! *****************************************************************************
!  FUNCTION bitfield_get(idx, bitfield) RESULT(is_set)
!    INTEGER(kind=int_4), INTENT(in)          :: idx
!    INTEGER(kind=int_4), DIMENSION(:), &
!      INTENT(in)                             :: bitfield
!    LOGICAL                                  :: is_set
!
!    INTEGER                                  :: mask, offset, val
!
!! Caution: This code only works with 4 byte integers
!
!     offset = ISHFT(idx, -5) + 1
!     val = IAND(idx, 31)
!     mask = 2**val
!
!     IF(offset > SIZE(bitfield)) THEN
!       WRITE (*,*) "idx: ", idx, "offset: ",offset, "size(bitfield): ",SIZE(bitfield)
!       STOP "set_bitfield: overflow"
!     END IF
!     IF(offset < 1) THEN
!       WRITE (*,*) "idx: ", idx, "offset: ",offset, "size(bitfield): ",SIZE(bitfield)
!       STOP "set_bitfield: underflow"
!     END IF
!
!     is_set = IAND(bitfield(offset), mask) == mask
!  END FUNCTION bitfield_get
!
END MODULE glbopt_walker

