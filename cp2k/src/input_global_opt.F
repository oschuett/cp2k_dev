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
  USE f77_blas
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: integer_t, real_t
  USE input_constants
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_global_opt'
  PUBLIC :: create_global_opt_section

CONTAINS

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE create_global_opt_section(section,error)
    TYPE(section_type), POINTER              :: section, subsection, print_key
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_global_opt_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

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

    CALL keyword_create(keyword, name="EMIN",&
        description="balba TODO.",&
        usage="EMIN <REAL>", &
        type_of_var=real_t,unit_str="hartree",required=.TRUE.,error=error)
    CALL section_add_keyword(section, keyword, error=error)
    CALL keyword_release(keyword, error=error)


    NULLIFY(subsection)
    CALL section_create(subsection,name="PRINT",&
               description="Controls the printing properties during a global optimization run",&
               n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
               error=error)
    NULLIFY(print_key)
    CALL cp_print_key_section_create(print_key,"WALKER_RUN_INFO",&
               description="Controls the printing of the walker's basic information during the global optimization", &
               print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
               error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"MASTER_RUN_INFO",&
               description="Controls the printing of the masters's basic information during the global optimization", &
               print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
               error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_global_opt_section

END MODULE input_global_opt

