!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_input
  USE glbopt_minhop_input,             ONLY: glbopt_create_minhop_section
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length
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
  USE input_constants,                 ONLY: low_print_level, add_last_numeric


#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_input'

 PUBLIC :: create_global_opt_section
 
  CONTAINS 
 
 

!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE create_global_opt_section(glbopt_section,error)
    TYPE(section_type), POINTER              :: glbopt_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(section_type), POINTER              :: subsection, printkey
    TYPE(keyword_type), POINTER              :: keyword

    NULLIFY(glbopt_section, subsection, printkey, keyword)

    CALL section_create(glbopt_section,name="GLOBAL_OPT",&
         description="blah TODO.",&
         repeats=.FALSE., required=.TRUE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword, name="NUMBER_OF_WALKERS",&
        description="balba TODO.",&
        usage="NUMBER_OF_WALKERS <INTEGER>", &
        type_of_var=integer_t,required=.TRUE.,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="EMIN",&
        description="balba TODO.",&
        usage="EMIN <REAL>", &
        type_of_var=real_t,unit_str="hartree",required=.TRUE.,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL section_create(subsection,name="PRINT",&
               description="Controls the printing properties during a global optimization run",&
               n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
               error=error)

    CALL cp_print_key_section_create(printkey,"WALKER_RUN_INFO",&
               description="Controls the printing of the walker's basic information during the global optimization", &
               print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
               error=error)
    CALL section_add_subsection(subsection,printkey,error=error)
    CALL section_release(printkey,error=error)

    CALL cp_print_key_section_create(printkey,"MASTER_RUN_INFO",&
               description="Controls the printing of the masters's basic information during the global optimization", &
               print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
               error=error)
    CALL section_add_subsection(subsection,printkey,error=error)
    CALL section_release(printkey,error=error)

    CALL section_add_subsection(glbopt_section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL glbopt_create_minhop_section(subsection, error)
    CALL section_add_subsection(glbopt_section,subsection,error=error)
    CALL section_release(subsection,error=error)

    !CALL glbopt_walker_register_input(glbopt_section, error)
  END SUBROUTINE create_global_opt_section

END MODULE glbopt_input

