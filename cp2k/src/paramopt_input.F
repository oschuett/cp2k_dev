!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE paramopt_input
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_constants,                 ONLY: add_last_numeric,&
                                             low_print_level,&
                                             swarm_do_glbopt,&
                                             swarm_do_paramopt
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: integer_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'paramopt_input'

 PUBLIC :: paramopt_declare_input

  CONTAINS



!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE paramopt_declare_input(swarm_section,error)
    TYPE(section_type), POINTER              :: swarm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: paramopt_section

    NULLIFY(paramopt_section, keyword)

    CALL section_create(paramopt_section,name="PARAM_OPT",&
         description="blah TODO.",&
         repeats=.FALSE., supported_feature=.TRUE.,error=error)

    CALL section_add_subsection(swarm_section, paramopt_section, error=error)
    CALL section_release(paramopt_section,error=error)

  END SUBROUTINE paramopt_declare_input
END MODULE paramopt_input

