!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE swarm_input
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
  USE glbopt_input,                    ONLY: glbopt_declare_input
  USE paramopt_input,                  ONLY: paramopt_declare_input
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_input'

 PUBLIC :: create_swarm_section

  CONTAINS


!TODO: maybe not put default values

!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE create_swarm_section(swarm_section,error)
    TYPE(section_type), POINTER              :: swarm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_section, printkey

    NULLIFY(swarm_section, print_section, printkey, keyword)

    CALL section_create(swarm_section,name="SWARM",&
         description="blah TODO.",&
         repeats=.FALSE., required=.TRUE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword, name="BEHAVIOR",&
         description="Which behaviour should control the swarm.",&
         usage="BEHAVIOR <STRING>",&
         required=.TRUE.,&
         default_i_val=swarm_do_glbopt,&
         enum_c_vals=s2a("GLOBAL_OPT","PARAM_OPT"),&
         enum_desc=s2a("Runs global geometry optimisation",&
                       "Runs global parameter optimisation"),&
         enum_i_vals=(/swarm_do_glbopt, swarm_do_paramopt/),&
         error=error)
    CALL section_add_keyword(swarm_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NUMBER_OF_WORKERS",&
        description="balba TODO.",&
        type_of_var=integer_t,required=.TRUE.,error=error)
    CALL section_add_keyword(swarm_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="REPLAY_COMMUNICATION_LOG",&
           description="TODO",repeats=.FALSE.,&
           usage="REPLAY_COMMUNICATION_LOG <CHARACTER>", default_lc_val="swarm_translog_replay.xyz", &
           error=error)
    CALL section_add_keyword(swarm_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_ITER",&
        description="The maximum number iterations the master should perform",&
        type_of_var=integer_t,default_i_val=HUGE(1),error=error)
    CALL section_add_keyword(swarm_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL section_create(print_section,name="PRINT",&
         description="Controls the printing properties during a global optimization run",&
         n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,error=error)

    CALL cp_print_key_section_create(printkey,"WORKER_RUN_INFO",&
               description="Controls the printing of the worker's basic information during the global optimization", &
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

    CALL cp_print_key_section_create(printkey,"COMMUNICATION_LOG",&
            description="Log all the communication between workers and master. Needed for restart.",&
            print_level=low_print_level, common_iter_levels=1,&
            filename="",unit_str="angstrom",error=error)
    CALL section_add_subsection(print_section,printkey,error=error)
    CALL section_release(printkey,error=error)

    CALL section_add_subsection(swarm_section,print_section,error=error)
    CALL section_release(print_section,error=error)


    CALL glbopt_declare_input(swarm_section, error)
    CALL paramopt_declare_input(swarm_section, error)

    !TODO: needed?
    !CALL section_release(swarm_section,error=error)
  END SUBROUTINE create_swarm_section


END MODULE swarm_input

