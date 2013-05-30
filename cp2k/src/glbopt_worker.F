!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_worker
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
                                             glbopt_mdctrl_data_type
  USE swarm_types,                     ONLY: swarm_cmd_shutdown,&
                                             swarm_command_type,&
                                             swarm_report_type
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
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_worker'

 PUBLIC :: glbopt_worker_init, glbopt_worker_finalize
 PUBLIC :: glbopt_worker_execute
 PUBLIC :: glbopt_worker_type

  TYPE glbopt_worker_type
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
 END TYPE glbopt_worker_type


 CONTAINS




! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_worker_init(worker, para_env, root_section,&
                 input_path, worker_id, error)
    TYPE(glbopt_worker_type), INTENT(INOUT)  :: worker
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(in)                      :: worker_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=default_path_length)       :: output_path
    CHARACTER(len=default_string_length)     :: new_project_name, &
                                                project_name, worker_name
    INTEGER                                  :: i
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: worker_section
    REAL(kind=dp)                            :: dist_in_angstrom
       failure = .FALSE.


        worker%error = error
        worker%root_section => root_section
        worker%iframe = 0
        worker%id = worker_id

       ! change the PROJECT_NAME to avoid clashes during backup of RESTART-files
       !output_path = "__STD_OUT__"
       IF(worker_id > 999) STOP "worker_create_f_env: Did not expect so many workers."
       CALL section_vals_val_get(root_section, "GLOBAL%PROJECT_NAME", c_val=project_name,error=error)
       WRITE(worker_name,"('WORKER',I3.3)"), worker_id
       new_project_name = TRIM(project_name)//"-"//worker_name
       CALL section_vals_val_set(root_section, "GLOBAL%PROJECT_NAME", c_val=new_project_name,error=error)
       output_path = TRIM(new_project_name)//".out"

      ! ======= Create f_env =======
      CALL create_force_env(worker%f_env_id, &
                            input_path=input_path, &
                            input=root_section, &
                            output_path=output_path, &
                            mpi_comm=para_env%group)


      ! ======= More setup stuff =======
      CALL f_env_add_defaults(worker%f_env_id, worker%f_env, worker%error, failure)
      IF(failure) STOP "glbopt_drive_worker: f_env_add_defaults failed"
      worker%force_env => worker%f_env%force_env
      CALL force_env_get(worker%force_env, globenv=worker%globenv, subsys=worker%subsys, error=worker%error)

      ! We want different random-number-streams for each worker
      DO i = 1, worker_id
         CALL reset_to_next_rng_substream(worker%globenv%gaussian_rng_stream, worker%error)
      END DO

      ! getting an output unit for logging
      logger => cp_error_get_logger(worker%error)
      worker%iw = cp_print_key_unit_nr(logger,root_section,&
          "GLBOPT%PRINT%WORKER_RUN_INFO",extension=".workerLog", error=worker%error)


      ! fetch original value from input
      CALL section_vals_val_get(root_section, "MOTION%GEO_OPT%MAX_ITER", i_val=worker%gopt_max_iter, error=worker%error)
      worker_section => section_vals_get_subs_vals(root_section, "GLBOPT%WORKER", error=error)

      CALL section_vals_val_get(worker_section,"MD_EPOT_N_UP", i_val=worker%epot_nup, error=error)
      CALL section_vals_val_get(worker_section,"MD_EPOT_N_DOWN", i_val=worker%epot_ndown, error=error)
      CALL section_vals_val_get(worker_section,"MD_MINCOUNT_MAX", i_val=worker%md_mincount_max, error=error)
      CALL section_vals_val_get(worker_section,"MD_MAX_DIST", r_val=dist_in_angstrom, error=error)
      worker%md_max_dist = dist_in_angstrom/angstrom
  END SUBROUTINE glbopt_worker_init



! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_worker_finalize(worker)
    TYPE(glbopt_worker_type), INTENT(INOUT)  :: worker

    INTEGER                                  :: ierr

     CALL f_env_rm_defaults(worker%f_env, worker%error)
     CALL destroy_force_env(worker%f_env_id, ierr)
     IF(ierr /= 0) STOP "glbopt_worker_finalize: destroy_force_env failed"
   END SUBROUTINE glbopt_worker_finalize


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE glbopt_worker_execute(worker, cmd, report, should_stop)
    TYPE(glbopt_worker_type), INTENT(INOUT)  :: worker
    TYPE(swarm_command_type), INTENT(IN)     :: cmd
    TYPE(swarm_report_type), INTENT(OUT)     :: report
    LOGICAL, INTENT(INOUT)                   :: should_stop

    INTEGER                                  :: n_atoms

       SELECT CASE (cmd%cmd_id)
         CASE (GLBOPT_CMD_SETIFRAME)
           worker%iframe = cmd%iframe
         CASE (GLBOPT_CMD_MDGOPT)
           CALL run_mdgopt(worker, cmd%temperature)
         CASE DEFAULT
           STOP "Worker: received unkown command"
       END SELECT

       ! assemble report
       report%worker_id = worker%id
       report%iframe = worker%iframe
       CALL force_env_get(worker%force_env, potential_energy=report%Epot, error=worker%error)
       IF (worker%iw>0) WRITE (worker%iw,'(A,21X,E20.10)')' GLBOPT| Potential Energy [Hartree]',report%Epot
       CALL cp_subsys_get(worker%subsys, natom=n_atoms, error=worker%error)
       ALLOCATE(report%positions(3*n_atoms))
       CALL pack_subsys_particles(worker%subsys, r=report%positions, error=worker%error)
   END SUBROUTINE glbopt_worker_execute


! *****************************************************************************
! *****************************************************************************
   SUBROUTINE run_mdgopt(worker, temperature)
    TYPE(glbopt_worker_type), INTENT(INOUT)  :: worker
    REAL(kind=dp), INTENT(IN)                :: temperature

    INTEGER                                  :: prev_iframe
    REAL(kind=dp)                            :: timestep
    TYPE(glbopt_mdctrl_data_type), TARGET    :: mdctrl_data
    TYPE(mdctrl_type), POINTER               :: mdctrl_p
    TYPE(mdctrl_type), TARGET                :: mdctrl

! setup mdctrl callback

      ALLOCATE(mdctrl_data%epot_history(worker%epot_ndown+worker%epot_nup+1))
      mdctrl_data%output_unit = worker%iw
      mdctrl_data%epot_history = 0.0
      mdctrl_data%md_mincounter = 0
      mdctrl_data%epot_nup = worker%epot_nup
      mdctrl_data%epot_ndown = worker%epot_ndown
      mdctrl_data%md_mincount_max = worker%md_mincount_max
      mdctrl_data%md_max_dist = worker%md_max_dist
      mdctrl%glbopt => mdctrl_data
      mdctrl_p => mdctrl


      !TODO: 300K is hard encoded
      timestep = SQRT(MIN(1.0, 300.0/(temperature*kelvin))) / femtoseconds
      IF (worker%iw>0) THEN
        WRITE (worker%iw,'(A,29X,F20.3)')  ' GLBOPT| MD temperature [K]',temperature*kelvin
        WRITE (worker%iw,'(A,31X,F20.3)')  ' GLBOPT| MD timestep [fs]',timestep*femtoseconds
        WRITE (worker%iw,'(A,25X,I10)') " GLBOPT| Starting MD at trajectory frame ", worker%iframe
      END IF

      mdctrl_data%epot_history = 0.0
      mdctrl_data%md_mincounter = 0
      prev_iframe = worker%iframe
      IF(worker%iframe==0) worker%iframe = 1 ! qs_mol_dyn behaves differently for STEP_START_VAL=0
      CALL section_vals_val_set(worker%root_section, "MOTION%MD%STEP_START_VAL", i_val=worker%iframe-1, error=worker%error)
      CALL section_vals_val_set(worker%root_section, "MOTION%MD%TEMPERATURE", r_val=temperature, error=worker%error)
      CALL section_vals_val_set(worker%root_section, "MOTION%MD%TIMESTEP", r_val=timestep, error=worker%error)

      CALL qs_mol_dyn(worker%force_env, worker%globenv, mdctrl=mdctrl_p, error=worker%error)

      worker%iframe = mdctrl_data%itimes + 1
      IF (worker%iw>0) WRITE (worker%iw,'(A,I4,A)') " GLBOPT| md ended after ", worker%iframe-prev_iframe, " steps."
      IF (worker%iw>0) WRITE (worker%iw,'(A,9X,I10)') " GLBOPT| Starting local optimisation at trajectory frame ", worker%iframe
      CALL section_vals_val_set(worker%root_section, "MOTION%GEO_OPT%STEP_START_VAL", i_val=worker%iframe-1, error=worker%error)
      CALL section_vals_val_set(worker%root_section, "MOTION%GEO_OPT%MAX_ITER",&
                                 i_val=worker%iframe+worker%gopt_max_iter, error=worker%error)

      CALL cp_geo_opt(worker%force_env, worker%globenv, rm_restart_info=.FALSE., error=worker%error)

      prev_iframe = worker%iframe
      CALL section_vals_val_get(worker%root_section, "MOTION%GEO_OPT%STEP_START_VAL", i_val=worker%iframe, error=worker%error)
      worker%iframe = worker%iframe + 2 ! Compensates for different START_VAL interpretation.
      IF (worker%iw>0) WRITE (worker%iw,'(A,I4,A)') " GLBOPT| gopt ended after ", worker%iframe-prev_iframe-1, " steps."

   END SUBROUTINE run_mdgopt


END MODULE glbopt_worker

