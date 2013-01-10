!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_env_types

#include "cp_common_uses.h"
 USE kinds,                           ONLY: dp
 USE cp_para_types,                   ONLY: cp_para_env_type
  
 IMPLICIT NONE
 PRIVATE

 
 

  TYPE glbopt_env_type
     !PRIVATE  !TODO: make this private with getters & setters
     INTEGER                                    :: ref_count
     ! config parameters
     INTEGER                                    :: n_masters
     INTEGER                                    :: n_walkers
     INTEGER                                    :: n_hashes
     INTEGER                                    :: n_colvars
     INTEGER                                    :: bloom_size
     INTEGER                                    :: epot_ndown
     INTEGER                                    :: epot_nup
     INTEGER                                    :: md_mincount_max
     INTEGER                                    :: hop_count_max
     ! stuff passed on for master/slave communication
     TYPE(cp_para_env_type), POINTER            :: full_para_env
     INTEGER                                    :: my_walker_id
     INTEGER                                    :: my_master_id
     ! stuff needed for md_callback
     INTEGER                                    :: md_mincounter
     REAL(KIND=dp), DIMENSION(:), POINTER       :: epot_history
     !mixed
     
  END TYPE glbopt_env_type

  
 PUBLIC :: glbopt_env_type, retain_glbopt_env, release_glbopt_env
  
 CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_env_types' 
CONTAINS
	
! *****************************************************************************
!> \brief retains the given global opt env
!> \param glbopt_env the global opt environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      see doc/ReferenceCounting.html
!> \par History
!>      11.2011 created [ole]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE retain_glbopt_env(glbopt_env, error)
    TYPE(glbopt_env_type), POINTER           :: glbopt_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'glbopt_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(glbopt_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(glbopt_env%ref_count>0,cp_failure_level,routineP,error,failure)
       glbopt_env%ref_count = glbopt_env%ref_count + 1
    END IF
  END SUBROUTINE retain_glbopt_env
  

! *****************************************************************************
!> \brief releases the given global opt env
!> \param glbopt_env the global opt environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      see doc/ReferenceCounting.html
!> \par History
!>      11.2012 created [ole]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE release_glbopt_env(glbopt_env, error)
    TYPE(glbopt_env_type), POINTER           :: glbopt_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'glbopt_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(glbopt_env)) THEN
       CPPrecondition(glbopt_env%ref_count>0,cp_failure_level,routineP,error,failure)
       glbopt_env%ref_count = glbopt_env%ref_count - 1
       IF (glbopt_env%ref_count==0) THEN
          !NULLIFY(glbopt_env%foobar)
          DEALLOCATE(glbopt_env,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE release_glbopt_env
  
 
END MODULE glbopt_env_types

