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
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type,&
                                             pack_subsys_particles
  USE f77_interface,                   ONLY: create_force_env,&
                                             destroy_force_env,&
                                             f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE geo_opt,                         ONLY: cp_geo_opt
  USE glbopt_types,                    ONLY: glbopt_cmd_mdgopt,&
                                             glbopt_cmd_setiframe,&
                                             glbopt_cmd_shutdown,&
                                             glbopt_command_type,&
                                             glbopt_mdctrl_data_type,&
                                             glbopt_report_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE md_run,                          ONLY: qs_mol_dyn
  USE mdctrl_types,                    ONLY: mdctrl_type
  USE parallel_rng_types,              ONLY: reset_to_next_rng_substream
  USE physcon,                         ONLY: angstrom,&
                                             femtoseconds,&
                                             kelvin
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_walker'

 PUBLIC :: glbopt_walker_init, glbopt_walker_finalize
 PUBLIC :: glbopt_walker_execute_command, glbopt_walker_assemble_report
 PUBLIC :: glbopt_walker_type

  TYPE glbopt_walker_type
   PRIVATE
   INTEGER                                  :: id
   INTEGER                                  :: iw
   TYPE(cp_error_type)                      :: error
   INTEGER                                  :: f_env_id
   TYPE(f_env_type), POINTER                :: f_env
   TYPE(force_env_type), POINTER            :: force_env
   TYPE(cp_subsys_type), POINTER            :: subsys
   INTEGER                                  :: iframe !indicates the NEXT frame to be written
   TYPE(section_vals_type), POINTER         :: root_section
   TYPE(global_environment_type), POINTER   :: globenv
   INTEGER                                  :: gopt_max_iter
   INTEGER                                  :: epot_ndown
   INTEGER                                  :: epot_nup
   INTEGER                                  :: md_mincount_max
   REAL(KIND=dp)                            :: md_max_dist
 END TYPE glbopt_walker_type

 CONTAINS




! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_walker_init(walker, para_env, root_section,&
                 input_path, walker_id, error)
    TYPE(glbopt_walker_type), INTENT(INOUT)  :: walker
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(in)                      :: walker_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=default_path_length)       :: output_path
    CHARACTER(len=default_string_length)     :: new_project_name, &
                                                project_name, walker_name
    INTEGER                                  :: i
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

       failure = .FALSE.

        walker%error = error
        walker%root_section => root_section
        walker%iframe = 0
        walker%id = walker_id

       ! change the PROJECT_NAME to avoid clashes during backup of RESTART-files
       !output_path = "__STD_OUT__"
       IF(walker_id > 999) STOP "walker_create_f_env: Did not expect so many walkers."
       CALL section_vals_val_get(root_section, "GLOBAL%PROJECT_NAME", c_val=project_name,error=error)
       WRITE(walker_name,"('WALKER',I3.3)"), walker_id
       new_project_name = TRIM(project_name)//"-"//walker_name
       CALL section_vals_val_set(root_section, "GLOBAL%PROJECT_NAME", c_val=new_project_name,error=error)
       output_path = TRIM(new_project_name)//".out"

      ! ======= Create f_env =======
      CALL create_force_env(walker%f_env_id, &
                            input_path=input_path, &
                            input=root_section, &
                            output_path=output_path, &
                            mpi_comm=para_env%group)


      ! ======= More setup stuff =======
      CALL f_env_add_defaults(walker%f_env_id, walker%f_env, walker%error, failure)
      IF(failure) STOP "glbopt_drive_walker: f_env_add_defaults failed"
      walker%force_env => walker%f_env%force_env
      CALL force_env_get(walker%force_env, globenv=walker%globenv, subsys=walker%subsys, error=walker%error)

      ! We want different random-number-streams for each walker
      DO i = 1, walker_id
         CALL reset_to_next_rng_substream(walker%globenv%gaussian_rng_stream, walker%error)
      END DO

      ! getting an output unit for logging
      logger => cp_error_get_logger(walker%error)
      walker%iw = cp_print_key_unit_nr(logger,root_section,&
          "GLOBAL_OPT%PRINT%WALKER_RUN_INFO",extension=".walkerLog", error=walker%error)


      ! fetch original value from input
      CALL section_vals_val_get(root_section, "MOTION%GEO_OPT%MAX_ITER", i_val=walker%gopt_max_iter, error=walker%error)

      walker%epot_nup = 2
      walker%epot_ndown = 2
      walker%md_mincount_max = 3
      walker%md_max_dist = 3.0/angstrom

  END SUBROUTINE glbopt_walker_init



! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_walker_finalize(walker)
    TYPE(glbopt_walker_type), INTENT(INOUT)  :: walker

    INTEGER                                  :: ierr

     CALL f_env_rm_defaults(walker%f_env, walker%error)
     CALL destroy_force_env(walker%f_env_id, ierr)
     IF(ierr /= 0) STOP "glbopt_walker_finalize: destroy_force_env failed"
   END SUBROUTINE glbopt_walker_finalize


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_walker_execute_command(walker, cmd, should_stop)
    TYPE(glbopt_walker_type), INTENT(INOUT)  :: walker
    TYPE(glbopt_command_type), INTENT(IN)    :: cmd
    LOGICAL, INTENT(INOUT)                   :: should_stop

       SELECT CASE (cmd%cmd_id)
         CASE (GLBOPT_CMD_SETIFRAME)
           walker%iframe = cmd%iframe
         CASE (GLBOPT_CMD_MDGOPT)
           CALL run_mdgopt(walker, cmd%temperature)

         CASE (GLBOPT_CMD_SHUTDOWN)
           IF(walker%iw>0) WRITE(walker%iw,*) "GLOBAL_OPT| Received shutdown command, quitting."
           should_stop = .TRUE.
         CASE DEFAULT
           STOP "Walker: received unkown command"
       END SELECT

   END SUBROUTINE glbopt_walker_execute_command

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_walker_assemble_report(walker, report)
    TYPE(glbopt_walker_type), INTENT(INOUT)  :: walker
    TYPE(glbopt_report_type), INTENT(OUT)    :: report

    INTEGER                                  :: n_atoms

       report%walker_id = walker%id
       report%iframe = walker%iframe
       CALL force_env_get(walker%force_env, potential_energy=report%Epot, error=walker%error)
       IF (walker%iw>0) WRITE (walker%iw,'(A,21X,E20.10)')' GLOBAL_OPT| Potential Energy [Hartree]',report%Epot
       CALL cp_subsys_get(walker%subsys, natom=n_atoms, error=walker%error)
       ALLOCATE(report%positions(3*n_atoms))
       CALL pack_subsys_particles(walker%subsys, r=report%positions, error=walker%error)
   END SUBROUTINE glbopt_walker_assemble_report



! *****************************************************************************
! *****************************************************************************
   SUBROUTINE run_mdgopt(walker, temperature)
    TYPE(glbopt_walker_type), INTENT(INOUT)  :: walker
    REAL(kind=dp), INTENT(IN)                :: temperature

    INTEGER                                  :: prev_iframe
    REAL(kind=dp)                            :: timestep
    TYPE(glbopt_mdctrl_data_type), TARGET    :: mdctrl_data
    TYPE(mdctrl_type), POINTER               :: mdctrl_p
    TYPE(mdctrl_type), TARGET                :: mdctrl

! setup mdctrl callback

      ALLOCATE(mdctrl_data%epot_history(walker%epot_ndown+walker%epot_nup+1))
      mdctrl_data%output_unit = walker%iw
      mdctrl_data%epot_history = 0.0
      mdctrl_data%md_mincounter = 0
      mdctrl_data%epot_nup = walker%epot_nup
      mdctrl_data%epot_ndown = walker%epot_ndown
      mdctrl_data%md_mincount_max = walker%md_mincount_max
      mdctrl_data%md_max_dist = walker%md_max_dist
      mdctrl%glbopt => mdctrl_data
      mdctrl_p => mdctrl


      !TODO: 300K is hard encoded
      timestep = SQRT(MIN(1.0, 300.0/(temperature*kelvin))) / femtoseconds
      IF (walker%iw>0) THEN
        WRITE (walker%iw,'(A,29X,F20.3)')  ' GLOBAL_OPT| MD temperature [K]',temperature*kelvin
        WRITE (walker%iw,'(A,31X,F20.3)')  ' GLOBAL_OPT| MD timestep [fs]',timestep*femtoseconds
        WRITE (walker%iw,'(A,25X,I10)') " GLOBAL_OPT| Starting MD at trajectory frame ", walker%iframe
      END IF

      mdctrl_data%epot_history = 0.0
      mdctrl_data%md_mincounter = 0
      prev_iframe = walker%iframe
      IF(walker%iframe==0) walker%iframe = 1 ! qs_mol_dyn behaves differently for STEP_START_VAL=0
      CALL section_vals_val_set(walker%root_section, "MOTION%MD%STEP_START_VAL", i_val=walker%iframe-1, error=walker%error)
      CALL section_vals_val_set(walker%root_section, "MOTION%MD%TEMPERATURE", r_val=temperature, error=walker%error)
      CALL section_vals_val_set(walker%root_section, "MOTION%MD%TIMESTEP", r_val=timestep, error=walker%error)

      CALL qs_mol_dyn(walker%force_env, walker%globenv, mdctrl=mdctrl_p, error=walker%error)

      walker%iframe = mdctrl_data%itimes + 1
      IF (walker%iw>0) WRITE (walker%iw,'(A,I4,A)') " GLOBAL_OPT| md ended after ", walker%iframe-prev_iframe, " steps."
      IF (walker%iw>0) WRITE (walker%iw,'(A,9X,I10)') " GLOBAL_OPT| Starting local optimisation at trajectory frame ", walker%iframe
      CALL section_vals_val_set(walker%root_section, "MOTION%GEO_OPT%STEP_START_VAL", i_val=walker%iframe-1, error=walker%error)
      CALL section_vals_val_set(walker%root_section, "MOTION%GEO_OPT%MAX_ITER",&
                                 i_val=walker%iframe+walker%gopt_max_iter, error=walker%error)

      CALL cp_geo_opt(walker%force_env, walker%globenv, rm_restart_info=.FALSE., error=walker%error)

      prev_iframe = walker%iframe
      CALL section_vals_val_get(walker%root_section, "MOTION%GEO_OPT%STEP_START_VAL", i_val=walker%iframe, error=walker%error)
      walker%iframe = walker%iframe + 2 ! Compensates for different START_VAL interpretation.
      IF (walker%iw>0) WRITE (walker%iw,'(A,I4,A)') " GLOBAL_OPT| gopt ended after ", walker%iframe-prev_iframe-1, " steps."

   END SUBROUTINE run_mdgopt


END MODULE glbopt_walker

