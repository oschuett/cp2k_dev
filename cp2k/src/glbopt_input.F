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
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length, dp
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


!TODO: maybe not put default values

!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE create_global_opt_section(glbopt_section,error)
    TYPE(section_type), POINTER              :: glbopt_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(section_type), POINTER              :: print_section, printkey
    TYPE(keyword_type), POINTER              :: keyword

    NULLIFY(glbopt_section, print_section, printkey, keyword)

    CALL section_create(glbopt_section,name="GLOBAL_OPT",&
         description="blah TODO.",&
         repeats=.FALSE., required=.TRUE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword, name="NUMBER_OF_WALKERS",&
        description="balba TODO.",&
        type_of_var=integer_t,required=.TRUE.,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="E_MIN",&
        description="balba TODO.",&
        default_r_val=-1*HUGE(1.0_dp),type_of_var=real_t,unit_str="hartree",error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL section_create(print_section,name="PRINT",&
         description="Controls the printing properties during a global optimization run",&
         n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,error=error)

    CALL cp_print_key_section_create(printkey,"WALKER_RUN_INFO",&
               description="Controls the printing of the walker's basic information during the global optimization", &
               print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
               error=error)
    CALL section_add_subsection(print_section,printkey,error=error)
    CALL section_release(printkey,error=error)

    CALL cp_print_key_section_create(printkey,"MASTER_RUN_INFO",&
               description="Controls the printing of the masters's basic information during the global optimization", &
               print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
               error=error)
    CALL section_add_subsection(print_section,printkey,error=error)
    CALL section_release(printkey,error=error)

    CALL section_add_subsection(glbopt_section,print_section,error=error)
    CALL section_release(print_section,error=error)

    CALL glbopt_master_declare_input(glbopt_section, error)
  END SUBROUTINE create_global_opt_section

!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE glbopt_master_declare_input(glbopt_section,error)
    TYPE(section_type), POINTER              :: glbopt_section
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(keyword_type), POINTER              :: keyword
    NULLIFY(keyword)

    CALL keyword_create(keyword, name="MAX_ITER",&
        description="The maximum number iterations the master should perform",&
        type_of_var=integer_t,default_i_val=HUGE(1),error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL glbopt_minhop_declare_input(glbopt_section, error)
  END SUBROUTINE glbopt_master_declare_input


!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE glbopt_minhop_declare_input(glbopt_section,error)
    TYPE(section_type), POINTER              :: glbopt_section
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(section_type), POINTER              :: minhop_section
    TYPE(keyword_type), POINTER              :: keyword

    NULLIFY(minhop_section, keyword)

    CALL section_create(minhop_section,name="MINIMA_HOPPING",&
         description="blah TODO.",&
         repeats=.FALSE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword, name="BETA_1",&
         description="balba TODO.",&
         default_r_val=1.1_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="BETA_2",&
         description="balba TODO.", &
         default_r_val=1.1_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="BETA_3",&
         description="balba TODO.", &
         default_r_val=1.0/1.1_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="ALPHA_1",&
         description="balba TODO.", &
         default_r_val=0.98_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="ALPHA_2",&
         description="balba TODO.", &
         default_r_val=1.0/0.98_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="E_ACCEPT_INIT",&
         description="Initial value of acceptance Energy",  &
         default_r_val=0.005_dp, type_of_var=real_t,unit_str="hartree",error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="TEMPERATURE_INIT",&
         description="Initially used temperature in Kelvin",  &
         default_r_val=100.0_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="TEMPERATURE_THROTTLE",&
         description="The maximum temperature used in Kelvin",  &
         default_r_val=100000.0_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="SHARE_HISTORY",&
         description="blab TODO",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(minhop_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(glbopt_section,minhop_section,error=error)
    CALL section_release(minhop_section,error=error)
  END SUBROUTINE glbopt_minhop_declare_input

END MODULE glbopt_input

