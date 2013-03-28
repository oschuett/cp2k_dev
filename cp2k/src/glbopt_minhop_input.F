!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_minhop_input
  USE kinds,                           ONLY: dp,&
                                             int_4
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_type
  USE physcon,                         ONLY: kelvin, femtoseconds
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
  !USE input_constants


#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_minhop_input'


 PUBLIC :: glbopt_create_minhop_section

 CONTAINS

!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE glbopt_create_minhop_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(keyword_type), POINTER              :: keyword

    NULLIFY(section, keyword)
    CALL section_create(section,name="MINIMA_HOPPING",&
         description="blah TODO.",&
         repeats=.FALSE., required=.TRUE., supported_feature=.TRUE.,error=error)
    !NULLIFY(keyword)
    !CALL keyword_create(keyword, name="EMIN",&
    !    description="balba TODO.",&
    !    usage="EMIN <REAL>", &
    !    type_of_var=real_t,unit_str="hartree",required=.TRUE.,error=error)
    !CALL section_add_keyword(section, keyword, error=error)
    !CALL keyword_release(keyword, error=error)
  END SUBROUTINE glbopt_create_minhop_section

 
END MODULE glbopt_minhop_input

