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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_interface,                   ONLY: create_force_env
  USE glbopt_master,                   ONLY: glbopt_drive_master
  USE glbopt_types,                    ONLY: glbopt_config_type
  USE glbopt_walker,                   ONLY: glbopt_drive_walker
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length
  USE message_passing,                 ONLY: mp_comm_free,&
                                             mp_comm_split,&
                                             mp_comm_split_direct,&
                                             mp_environ
  USE timings,                         ONLY: timeset,&
                                             timestop

  USE physcon,                         ONLY: angstrom



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
   SUBROUTINE run_global_opt(root_section, para_env, input_path, output_unit, error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(glbopt_config_type)                 :: glbopt_cfg

    TYPE(section_vals_type), POINTER         :: colvar_section

    CALL section_vals_val_get(root_section,"GLOBAL_OPT%NUMBER_OF_WALKERS",&
       i_val=glbopt_cfg%n_walkers,error=error)
    CALL section_vals_val_get(root_section,"GLOBAL_OPT%EMIN",&
       r_val=glbopt_cfg%Emin,error=error)
    !TODO: read from input or make educated guess
    ! TODO: more hierarchic naming scheme eg. bloom_* md_* ...
    glbopt_cfg%n_hashes = 10
    glbopt_cfg%bloom_size = 2**10
    glbopt_cfg%epot_ndown = 2
    glbopt_cfg%epot_nup = 2
    glbopt_cfg%md_mincount_max = 3
    glbopt_cfg%hop_count_max = 100
    glbopt_cfg%beta1  = 1.1
    ! glbopt_cfg%beta2  = 1.1
    glbopt_cfg%beta3  = 1.0 / 1.1
    glbopt_cfg%Eaccept0 = 0.005
    glbopt_cfg%alpha1 = 0.98
    glbopt_cfg%alpha2 = 1 / 0.98
    glbopt_cfg%md_max_dist   = 3.0/angstrom
    glbopt_cfg%gopt_max_dist = 3.0/angstrom

    ! make sure there are colvars
    colvar_section => section_vals_get_subs_vals(root_section,"FORCE_EVAL%SUBSYS%COLVAR",error=error)
    CALL section_vals_get(colvar_section,n_repetition=glbopt_cfg%n_colvars,error=error)
    IF(glbopt_cfg%n_colvars == 0) STOP "No colvars found. Aborting."

    CALL glbopt_driver(glbopt_cfg, root_section, para_env, input_path, output_unit, error)

  END SUBROUTINE run_global_opt

! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
   SUBROUTINE glbopt_driver(glbopt_cfg, root_section, para_env, input_path, output_unit, error)
    TYPE(glbopt_config_type)                 :: glbopt_cfg
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER :: f_env_id, handle
    CHARACTER(len=*), PARAMETER :: routineN = 'glbopt_driver', &
      routineP = moduleN//':'//routineN
    CALL timeset(routineN,handle)

    IF(output_unit > 0) WRITE(output_unit,"(A)") " GLOBAL_OPT| Reporting for duty :-)"

    IF(glbopt_cfg%n_walkers == 1) THEN
       IF(output_unit > 0) WRITE(output_unit,"(A)") " GLOBAL_OPT| Running in single walker mode."
       CALL glbopt_drive_walker(glbopt_cfg, root_section, input_path, output_unit, para_env, para_env%group, 0, error)
    ELSE
       IF(output_unit > 0) WRITE(output_unit,"(A)") " GLOBAL_OPT| Running in master / walkers mode."
       CALL driver_master_walkers(glbopt_cfg, root_section, para_env, input_path, output_unit, error)
    ENDIF

    CALL timestop(handle)
   END SUBROUTINE glbopt_driver

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE driver_master_walkers(glbopt_cfg, root_section, para_env, input_path, output_unit, error)
    TYPE(glbopt_config_type)                 :: glbopt_cfg
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER :: f_env_id, my_walker_id, n_groups_created, &
      pe_per_walker, subgroup, subgroup_rank, subgroup_size, walker_group
    LOGICAL                                  :: im_the_master
    INTEGER, DIMENSION(:), POINTER           :: group_distribution_p
    INTEGER, &
      DIMENSION(0:para_env%num_pe-2), TARGET :: group_distribution


    ! ====== Setup of MPI-Groups ======
    IF (MOD(para_env%num_pe-1, glbopt_cfg%n_walkers) /= 0) &
       STOP "number of processors-1 is not divisible by n_walkers."
    IF (para_env%num_pe < glbopt_cfg%n_walkers + 1) &
       STOP "There are not enough processes for n_walkers + 1. Aborting."

    pe_per_walker = (para_env%num_pe-1)/glbopt_cfg%n_walkers

    WRITE(*,*) "!!! Number of mpi ranks: ",para_env%num_pe
    WRITE(*,*) "!!! Number of walkers: ", glbopt_cfg%n_walkers
    WRITE(*,*) "!!! Number of colvars: ", glbopt_cfg%n_colvars

    ! the last task becomes the master. Keeps node-alignment of other tasks.
    im_the_master = (para_env%mepos == para_env%num_pe-1)

    ! First split split para_env%group into a master- and a walkers-groups...
    IF (im_the_master) THEN
       CALL mp_comm_split_direct(para_env%group, subgroup, 1)
       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
       !WRITE (*,*) "Hello, this is a Master - there are ",subgroup_size, " of us."
       IF(subgroup_size/=1) STOP "glbopt: mp_comm_split_direct failed (master)"
    ELSE
       CALL mp_comm_split_direct(para_env%group, subgroup, 2)
       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
       !WRITE (*,*) "Hello, this is a Walker - there are ",subgroup_size, " of us."
       IF(subgroup_size/=para_env%num_pe-1) STOP "glbopt: mp_comm_split_direct failed (walker)"
    ENDIF


    IF (im_the_master) THEN
       CALL glbopt_drive_master(glbopt_cfg, root_section, para_env, error)
    ELSE
       ! ...then split walkers-group into n_walkers groups - one for each walker.
       group_distribution_p => group_distribution
       CALL mp_comm_split(subgroup, walker_group, n_groups_created, group_distribution_p, n_subgroups=glbopt_cfg%n_walkers)
       my_walker_id = group_distribution(subgroup_rank)
       IF(n_groups_created/=glbopt_cfg%n_walkers) STOP "glbopt: mp_comm_split failed."
       ! f_env_id = walker_create_f_env(root_section, input_path, walker_group, my_walker_id, error)
       CALL glbopt_drive_walker(glbopt_cfg, root_section, input_path, output_unit, para_env, walker_group, my_walker_id, error)
       CALL mp_comm_free(walker_group)
    ENDIF

    CALL mp_comm_free(subgroup)

  END SUBROUTINE driver_master_walkers


!! *****************************************************************************
!  FUNCTION walker_create_f_env(root_section, input_path, mpi_comm, walker_id, error) RESULT(f_env_id)
!    TYPE(section_vals_type), POINTER         :: root_section
!    CHARACTER(LEN=*), INTENT(IN)             :: input_path
!    INTEGER, INTENT(IN)                      :: mpi_comm, walker_id
!    TYPE(cp_error_type), INTENT(inout)       :: error
!    INTEGER                                  :: ierr, f_env_id
!
!    CHARACTER(len=default_path_length)       :: output_path
!    CHARACTER(len=default_string_length)     :: new_project_name, &
!                                                project_name, walker_name
!
!! change the PROJECT_NAME to avoid clashes during backup of RESTART-files
!    IF(walker_id > 999) STOP "walker_create_f_env: Did not expect so many walkers."
!    CALL section_vals_val_get(root_section, "GLOBAL%PROJECT_NAME", c_val=project_name,error=error)
!    WRITE(walker_name,"('WALKER',I3.3)"), walker_id
!    !walker_name = "walker"//ADJUSTL(cp_to_string(walker_id))
!    new_project_name = TRIM(project_name)//"-"//walker_name
!    CALL section_vals_val_set(root_section, "GLOBAL%PROJECT_NAME", c_val=new_project_name,error=error)
!
!    output_path = TRIM(new_project_name)//".out"
!    WRITE (*,*) output_path
!
!    CALL create_force_env(f_env_id, &
!                          input_path=input_path, &
!                          input=root_section, &
!                          output_path=output_path, &
!                          mpi_comm=mpi_comm,&
!                          ierr=ierr)
!
!    IF(ierr/=0) STOP "glbopt: failed to create force_env"
!
!  END FUNCTION walker_create_f_env
!
!
END MODULE glbopt

