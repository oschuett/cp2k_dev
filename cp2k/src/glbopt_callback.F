!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_callback
  USE glbopt_types,                    ONLY: glbopt_mdctrl_data_type
  USE md_ener_types,                   ONLY: md_ener_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type

  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE particle_list_types,             ONLY: particle_list_type
  USE kinds,                           ONLY: dp



#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_callback'

 PUBLIC :: glbopt_md_callback

  
 CONTAINS 
 
 
  
! *****************************************************************************
!> \brief Main driver to perform global optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE glbopt_md_callback(mdctrl_data, md_env, should_stop, error)
    TYPE(glbopt_mdctrl_data_type), POINTER   :: mdctrl_data
    TYPE(md_environment_type), POINTER       :: md_env
    LOGICAL, INTENT(inout)                   :: should_stop
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'glbopt_md_callback', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure = .FALSE., &
                                                passed_minimum
    INTEGER, POINTER                         :: itimes
    TYPE(md_ener_type), POINTER              :: md_ener

! ====== Local Variables ======

    CPPrecondition(ASSOCIATED(mdctrl_data), cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(md_env), cp_failure_level,routineP,error,failure)



    ! add new potential energy value to history
    NULLIFY(md_ener, itimes)
    CALL get_md_env(md_env=md_env, md_ener=md_ener, itimes=itimes, error=error)
    !WRITE (*,*) 'potential energy OLE',  glbopt_env%walker_num, md_ener%epot
    mdctrl_data%epot_history = EOSHIFT(mdctrl_data%epot_history, shift=-1)
    mdctrl_data%epot_history(1) = md_ener%epot
    !DO i=1,SIZE(glbopt_env%epot_history,1)
    !  WRITE (*,*), glbopt_env%epot_history(i)
    !END DO

    !WRITE (*,*) "hello word - glbopt_md_callback epot: ", md_ener%epot

    ! check if we passed a minimum
    passed_minimum = .TRUE.
    DO i=1, mdctrl_data%cfg%epot_nup
      IF(mdctrl_data%epot_history(i) <= mdctrl_data%epot_history(i+1)) &
         passed_minimum = .FALSE.
    END DO

    DO i=mdctrl_data%cfg%epot_nup+1, mdctrl_data%cfg%epot_nup+mdctrl_data%cfg%epot_ndown
      IF(mdctrl_data%epot_history(i) >= mdctrl_data%epot_history(i+1)) &
         passed_minimum = .FALSE.
    END DO

    !WRITE (*,*) "passed minima : ", mdctrl_data%md_mincounter

    ! count the passed minima and stop md_run when mincount_max is reached.
    IF(passed_minimum) &
       mdctrl_data%md_mincounter = mdctrl_data%md_mincounter + 1

    IF(mdctrl_data%md_mincounter >= mdctrl_data%cfg%md_mincount_max) THEN
       should_stop = .TRUE.
       IF(mdctrl_data%output_unit > 0) WRITE (mdctrl_data%output_unit,&
         "(A,I4,A)") " GLOBAL_OPT| Stoping MD because of MD_MINCOUNTER after ",itimes," steps."
    END IF

   ! IF(itimes > 5) THEN
   !    ! Check for particles that wander off (= with very weak forces)
   !    CALL get_md_env(md_env=md_env, force_env=force_env, error=error)
   !    CALL force_env_get(force_env, subsys=subsys, error=error)
   !    CALL cp_subsys_get(subsys,particles=particles,error=error)
   !    DO iparticle = 1, SIZE(particles%els)
   !       force = particles%els(iparticle)%f(:)
   !       force_norm = SQRT(DOT_PRODUCT(force,force)) 
   !       IF(force_norm < mdctrl_data%cfg%md_max_force) THEN
   !          IF(mdctrl_data%output_unit > 0) &
   !             WRITE (mdctrl_data%output_unit,"(A,I4,A)")&
   !               " GLOBAL_OPT| Stoping MD because of MD_MAX_FORCE after ", itimes," steps."
   !          should_stop = .TRUE.
   !       END IF
   !    END DO
   ! END IF


    CALL check_connectivity(md_env, error)

    ! stop when certain number of steps are reached.
    IF(itimes >= 30) THEN
       should_stop = .TRUE.
       IF(mdctrl_data%output_unit > 0) &
          WRITE (mdctrl_data%output_unit,"(A,I4,A)")&
            " GLOBAL_OPT| stoping MD because of STEPS after ",itimes," steps."
    END IF

  END SUBROUTINE glbopt_md_callback

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE check_connectivity(md_env, error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    REAL(KIND=dp)                            :: dr(3), dr_norm
    INTEGER                                  :: iparticle, jparticle

    CALL get_md_env(md_env=md_env, force_env=force_env, error=error)
    CALL force_env_get(force_env, subsys=subsys, error=error)
    CALL cp_subsys_get(subsys,particles=particles,error=error)
    DO iparticle = 1, SIZE(particles%els)
       DO jparticle = 1, SIZE(particles%els)
       dr = particles%els(iparticle)%r(:) - particles%els(jparticle)%r(:)
       dr_norm = SQRT(DOT_PRODUCT(dr,dr))
       WRITE (*,*) "dist ", iparticle, jparticle
       END DO
    END DO
  END SUBROUTINE check_connectivity

END MODULE glbopt_callback

