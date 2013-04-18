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
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_utils'

 PUBLIC :: particles_connected


 CONTAINS



! *****************************************************************************
! *****************************************************************************
  FUNCTION particles_connected(positions, max_distance) RESULT(all_connected)
    REAL(KIND=dp), DIMENSION(:)              :: positions
    REAL(KIND=dp)                            :: max_distance
    LOGICAL                                  :: all_connected

    INTEGER                                  :: iparticle, jparticle, &
                                                n_particles, stack_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: stack
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: marked
    REAL(KIND=dp)                            :: dr(3), dr_norm

    n_particles = SIZE(positions) / 3
    ALLOCATE(stack(n_particles), marked(n_particles))
    marked = .FALSE.; stack_size = 0

    ! First particle taken as root of flooding, mark it and push to stack
    marked(1) = .TRUE.; stack(1) = 1; stack_size = 1


    !WRITE (*,*) "max dist: ",max_distance
    DO WHILE (stack_size > 0)
       iparticle = stack(stack_size); stack_size=stack_size-1  !pop
       !WRITE (*,*) "working on particle: ",iparticle
       DO jparticle = 1, n_particles
         dr = positions(3*iparticle-2:3*iparticle) - positions(3*jparticle-2:3*jparticle)
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

