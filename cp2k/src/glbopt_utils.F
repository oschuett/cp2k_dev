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
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dp
  USE particle_list_types,             ONLY: particle_list_type
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

    INTEGER                                  :: iparticle, jparticle, &
                                                stack_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: stack
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: marked
    REAL(KIND=dp)                            :: dr(3), dr_norm
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

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

