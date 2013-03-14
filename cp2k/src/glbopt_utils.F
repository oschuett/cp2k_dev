!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_utils
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
                                             glbopt_msg_sync,&
                                             glbopt_msg_report,&
                                             glbopt_msg_shutdown
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
  USE physcon,                         ONLY: kelvin
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

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_utils'

 PUBLIC :: particles_connected


 CONTAINS



! *****************************************************************************
! *****************************************************************************
  FUNCTION particles_connected(force_env, max_distance, error) RESULT(all_connected)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp)                            :: max_distance
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: all_connected

    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    REAL(KIND=dp)                            :: dr(3), dr_norm
    INTEGER                                  :: iparticle, jparticle
    LOGICAL, DIMENSION(:), ALLOCATABLE       :: marked
    INTEGER, DIMENSION(:), ALLOCATABLE       :: stack
    INTEGER                                  :: stack_size

    CALL force_env_get(force_env, subsys=subsys, error=error)
    CALL cp_subsys_get(subsys,particles=particles,error=error)

    ALLOCATE(stack(SIZE(particles%els)), marked(SIZE(particles%els)))
    marked = .FALSE.; stack_size = 0

    ! First particle taken as root of flooding, mark it and push to stack
    marked(1) = .TRUE.; stack(1) = 1; stack_size = 1

    !WRITE (*,*) "max dist: ",max_distance
    DO WHILE (stack_size > 0)
       iparticle = stack(stack_size); stack_size=stack_size-1  !pop
       !WRITE (*,*) "working on particle: ",iparticle
       DO jparticle = 1, SIZE(particles%els)
         dr = particles%els(iparticle)%r(:) - particles%els(jparticle)%r(:)
         dr_norm = SQRT(DOT_PRODUCT(dr,dr))
         !WRITE (*,*) "dist ", iparticle, jparticle, dr_norm, dr_norm < max_distance
         IF(dr_norm < max_distance) THEN ! they are close = they are connected
            IF(.NOT. marked(jparticle)) THEN
                marked(jparticle) = .TRUE.
                stack(stack_size+1) = jparticle; stack_size=stack_size+1; !push
            END IF
         END IF
       END DO
    END DO

    all_connected = ALL(marked) !did we visit every particle?

    !IF(.NOT. all_connected) THEN
    !   WRITE (*,*) "GLOBAL_OPT: unconnected particles"
    !   DO iparticle = 1, SIZE(particles%els)
    !     IF(.NOT. marked(iparticle)) WRITE (*,*) "particle ", iparticle
    !   END DO
    !   STOP "particles not connected"
    !END IF

  END FUNCTION particles_connected

END MODULE glbopt_utils

