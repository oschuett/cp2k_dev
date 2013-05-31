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
  USE input_cp2k_subsys,               ONLY: create_coord_section

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_input'

 PUBLIC :: glbopt_declare_input

  CONTAINS



!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE glbopt_declare_input(swarm_section,error)
    TYPE(section_type), POINTER              :: swarm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: glbopt_section

    NULLIFY(glbopt_section, keyword)

    CALL section_create(glbopt_section,name="GLOBAL_OPT",&
         description="blah TODO.",&
         repeats=.FALSE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword, name="E_MIN",&
        description="balba TODO.",&
        default_r_val=-1*HUGE(1.0_dp),type_of_var=real_t,unit_str="hartree",error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)

    CALL keyword_release(keyword, error=error)
    CALL keyword_create(keyword, name="MD_EPOT_N_UP",&
        description="balba TODO.",&
        type_of_var=integer_t,default_i_val=2,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="MD_EPOT_N_DOWN",&
        description="balba TODO.",&
        type_of_var=integer_t,default_i_val=2,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="MD_MINCOUNT_MAX",&
        description="balba TODO.",&
        type_of_var=integer_t,default_i_val=3,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="MD_MAX_DIST",&
         description="balba TODO.", &
         default_r_val=2.0_dp, unit_str="angstrom",type_of_var=real_t,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    !CALL create_coord_section(glbopt_section, error=error)

    CALL minhop_declare_input(glbopt_section, error=error)

    CALL section_add_subsection(swarm_section, glbopt_section,error=error)
    CALL section_release(glbopt_section,error=error)
  END SUBROUTINE glbopt_declare_input

  
!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE minhop_declare_input(glbopt_section,error)
    TYPE(section_type), POINTER              :: glbopt_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: minhop_section

    NULLIFY(minhop_section, keyword)

    CALL section_create(minhop_section,name="MINIMA_HOPPING",&
         description="blah TODO.",&
         repeats=.FALSE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword, name="BETA_2",&
         description="balba TODO.We do not desctinguis between beta_1 and beta_2 as in the orig paper", &
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

    CALL keyword_create(keyword, name="TEMPERATURE_MAX",&
         description="The maximum temperature used in Kelvin",  &
         default_r_val=100000.0_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="TEMPERATURE_MIN",&
         description="The maximum temperature used in Kelvin",  &
         default_r_val=1.0_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="SHARE_HISTORY",&
         description="blab TODO",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(minhop_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_DIST",&
         description="balba TODO.", &
         default_r_val=2.0_dp, unit_str="angstrom",type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL section_add_subsection(glbopt_section,minhop_section,error=error)
    CALL section_release(minhop_section,error=error)
  END SUBROUTINE minhop_declare_input
END MODULE glbopt_input

