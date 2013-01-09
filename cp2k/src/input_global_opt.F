!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the input structure for global optimization
!> \par History
!>      11.2012 created [Ole]
!> \author Ole
! *****************************************************************************
MODULE input_global_opt
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE f77_blas
  USE input_constants,                 ONLY: low_print_level,&
                                             opt_force_matching
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: char_t,&
                                             real_t,&
                                             integer_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_global_opt'
  PUBLIC :: create_global_opt_section

CONTAINS

! *****************************************************************************
!> \brief creates the optimize_input section
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_global_opt_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_global_opt_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: sub_section, subsubsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="GLOBAL_OPT",&
         description="blah TODO.",&
         repeats=.FALSE., required=.TRUE., supported_feature=.TRUE.,error=error)
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="NUMBER_OF_WALKERS",&
        description="balba TODO.",&
        usage="NUMBER_OF_WALKERS <INTEGER>", &
        type_of_var=integer_t,required=.TRUE.,error=error)
    CALL section_add_keyword(section, keyword, error=error)
    CALL keyword_release(keyword, error=error)
    
  END SUBROUTINE create_global_opt_section

END MODULE input_global_opt

