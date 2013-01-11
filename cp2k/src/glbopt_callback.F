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
    TYPE(md_ener_type), POINTER              :: md_ener

! ====== Local Variables ======

    CPPrecondition(ASSOCIATED(mdctrl_data), cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(md_env), cp_failure_level,routineP,error,failure)
        
    !WRITE (*,*) "hello word - glbopt_md_callback"
    
    ! add new potential energy value to history
    NULLIFY(md_ener)
    CALL get_md_env(md_env=md_env, md_ener=md_ener, error=error)
    !WRITE (*,*) 'potential energy OLE',  glbopt_env%walker_num, md_ener%epot
    mdctrl_data%epot_history = EOSHIFT(mdctrl_data%epot_history, shift=-1)
    mdctrl_data%epot_history(1) = md_ener%epot
    !DO i=1,SIZE(glbopt_env%epot_history,1)
    !  WRITE (*,*), glbopt_env%epot_history(i)
    !END DO
    
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
    
    !WRITE (*,*) "passed minima : ", glbopt_env%md_mincounter
    
    ! count the passed minima and stop md_run when mincount_max is reached.
    IF(passed_minimum) &
       mdctrl_data%md_mincounter = mdctrl_data%md_mincounter + 1
    
    IF(mdctrl_data%md_mincounter >= mdctrl_data%cfg%md_mincount_max) &
       should_stop = .TRUE.
    
    
    
  END SUBROUTINE glbopt_md_callback
 
END MODULE glbopt_callback

