!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input section for MP2
!> \par History
!>      05.2011 created
!> \author MDB
! *****************************************************************************
MODULE input_cp2k_mp2
  USE bibliography,                    ONLY: DelBen2012
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE f77_blas
  USE input_constants
  USE input_cp2k_hfx,                  ONLY: create_hfx_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mp2'

  PUBLIC :: create_mp2_section

CONTAINS

! *****************************************************************************
!> \brief creates the input section for the mp2 part
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author MDB
! *****************************************************************************
  SUBROUTINE create_mp2_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mp2_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"WF_CORRELATION",&
            description="Sets up the Wavefunction-based Correlation parameters if requested ",&
            n_keywords=8, n_subsections=6, repeats=.TRUE., required=.TRUE.,&
            citations=(/DelBen2012/),&
            error=error)

       NULLIFY(keyword, print_key, subsection)

       CALL keyword_create(&
            keyword=keyword,&
            name="METHOD",&
            citations=(/DelBen2012/),&
            description="Which method should be used to compute the MP2 energy",&
            usage="METHOD MP2_GPW",&
            enum_c_vals=s2a("NONE","DIRECT_CANONICAL","MP2_GPW","RI_MP2_GPW","RI_RPA_GPW","RI_SOS_LAPLACE"),&
            enum_i_vals=(/mp2_method_none,mp2_method_direct,mp2_method_gpw,ri_mp2_method_gpw,ri_rpa_method_gpw,ri_mp2_laplace/),&
            enum_desc=s2a("Skip MP2 calculation",&
                          "Use the direct mp2 canonical approach",&
                          "Use the GPW approach to MP2",&
                          "Use the GPW approach to RI-MP2",&
                          "Use the GPW approach to RI-RPA",&
                          "Use the GPW approach to RI-SOS-Laplace-MP2"),&
            default_i_val=mp2_method_direct, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="MEMORY",&
            description="Maximum allowed total memory usage during MP2 methods [Mb].",&
            usage="MEMORY 1500 ",&
            default_r_val=1.024E+3_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="SCALE_S",&
            description="Scaling factor of the singlet energy component (opposite spin, OS). ",&
            usage="SCALE_S  1.0",&
            default_r_val=1.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
  
       CALL keyword_create(&
            keyword=keyword,&
            name="SCALE_T",&
            description="Scaling factor of the triplet energy component (same spin, SS).",&
            usage="SCALE_T  1.0",&
            default_r_val=1.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="GROUP_SIZE",&
            variants=(/"NUMBER_PROC"/),&
            description="Group size used in the computation of the integrals. Default is to use all processors (GROUP_SIZE=-1)."//&
             "A smaller group size (for example the node size), might a better choice if the actual MP2 time is large "//&
             "compared to integral computation time. This is usually the case if the total number of processors is not too large.",&
            usage="GROUP_SIZE 32",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="ROW_BLOCK",&
            variants=(/"ROW_BLOCK_SIZE"/),&
            description="Size of the row block used in the SCALAPACK block cyclic data distribution."//&
             "Default is (ROW_BLOCK=-1) is automatic. "//&
             "A proper choice can speedup the parallel matrix multiplication in the case of RI-RPA and RI-SOS-MP2-Laplace.",&
            usage="ROW_BLOCK 512",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="COL_BLOCK",&
            variants=(/"COL_BLOCK_SIZE"/),&
            description="Size of the column block used in the SCALAPACK block cyclic data distribution."//&
             "Default is (COL_BLOCK=-1) is automatic. "//&
             "A proper choice can speedup the parallel matrix multiplication in the case of RI-RPA and RI-SOS-MP2-Laplace.",&
            usage="COL_BLOCK 512",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"MP2_INFO",&
            description="Controls the printing basic info about MP2 method", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL create_mp2_direct(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_wfc_gpw(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_ri_mp2(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_ri_rpa(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_ri_laplace(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_mp2_potential(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_mp2_section

  SUBROUTINE create_mp2_direct(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mp2_direct', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"DIRECT_CANONICAL",&
            description="Parameters influencing the direct canonical method",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(&
            keyword=keyword,&
            name="BIG_SEND",&
            description="Send big messages between processes (useful for >48 processors).",&
            usage="BIG_SEND",&
            default_l_val=.TRUE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_mp2_direct

  SUBROUTINE create_ri_mp2(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ri_mp2', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"RI_MP2",&
            description="Parameters influencing the RI MP2 method",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="BLOCK_SIZE",&
            variants=(/"MESSAGE_SIZE"/),&
            description="Determines the blocking used for communication in RI-MP2. Larger BLOCK_SIZE "//&
                        "reduces communication but requires more memory. The default (-1) is automatic.",&
            usage="BLOCK_SIZE 2",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_ri_mp2

  SUBROUTINE create_ri_laplace(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ri_laplace', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"RI_LAPLACE",&
            description="Parameters influencing the RI-SOS-MP2-Laplace method",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="QUADRATURE_POINTS",&
            variants=(/"LAPLACE_NUM_QUAD_POINTS"/),&
            description="Number of quadrature points for the numerical integration in the RI-SOS-MP2-Laplace method.",&
            usage="QUADRATURE_POINTS 6",&
            default_i_val=5,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SIZE_INTEG_GROUP",&
            variants=(/"LAPLACE_GROUP_SIZE"/),&
            description="Group size for the integration in the Laplace method, that is the number of processes involved in "//&
                        "the computation of each integration point. SIZE_INTEG_GROUP has to be a multiple "//&
                        "of GROUP_SIZE in the WF_CORRELATION section. The default (-1) "//&
                        "is automatic.",&
            usage="SIZE_INTEG_GROUP 16",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_ri_laplace

  SUBROUTINE create_ri_rpa(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ri_rpa', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"RI_RPA",&
            description="Parameters influencing the RI RPA method",&
            n_keywords=2, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword,subsection)
  
       CALL keyword_create(keyword, name="QUADRATURE_POINTS",&
            variants=(/"RPA_NUM_QUAD_POINTS"/),&
            description="Number of quadrature points for the numerical integration in the RI-RPA method.",&
            usage="QUADRATURE_POINTS 60",&
            default_i_val=40,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SIZE_FREQ_INTEG_GROUP",&
            variants=(/"RPA_GROUP_SIZE"/),&
            description="Group size for frequency integration, that is the number of processes involved in "//&
                        "the computation of each integration point. SIZE_FREQ_INTEG_GROUP has to be a multiple "//&
                        "of GROUP_SIZE in the WF_CORRELATION section. The default (-1) "//&
                        "is automatic.",&
            usage="SIZE_FREQ_INTEG_GROUP 16",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! here we generate a hfx subsection to use in the case EXX has to be computed after RPA
       CALL create_hfx_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_ri_rpa

  SUBROUTINE create_wfc_gpw(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_wfc_gpw', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"WFC_GPW",&
            description="Parameters for the GPW approach in Wavefunction-based Correlation methods",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="EPS_GRID",&
            description="Determines a threshold for the GPW based integration",&
            usage="EPS_GRID 1.0E-5 ",required=.FALSE.,type_of_var=real_t,&
            default_r_val=1.0E-6_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_FILTER",&
            description="Determines a threshold for the DBCSR based multiply (usually 10 times smaller than EPS_GRID).",&
            usage="EPS_FILTER 1.0E-6 ",required=.FALSE.,type_of_var=real_t,&
            default_r_val=1.0E-7_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CUTOFF",&
            description="The cutoff of the finest grid level in the MP2 gpw integration.",&
            usage="CUTOFF 300",required=.FALSE.,type_of_var=real_t,&
            default_r_val=300.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REL_CUTOFF",&
            variants=(/"RELATIVE_CUTOFF"/),&
            description="Determines the grid at which a Gaussian is mapped.",&
            usage="REL_CUTOFF 50",required=.FALSE.,type_of_var=real_t,&
            default_r_val=50.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRINT_LEVEL",&
            variants=(/"IOLEVEL"/),&
            description="How much output is written by the individual groups.",&
            usage="PRINT_LEVEL HIGH",&
            default_i_val=silent_print_level,enum_c_vals=&
             s2a("SILENT","LOW","MEDIUM","HIGH","DEBUG"),&
            enum_desc=s2a( "Almost no output",&
                           "Little output", "Quite some output", "Lots of output",&
                           "Everything is written out, useful for debugging purposes only"),&
            enum_i_vals=(/silent_print_level,low_print_level,medium_print_level,&
            high_print_level,debug_print_level/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_wfc_gpw


  SUBROUTINE create_mp2_potential(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mp2_potential', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"INTERACTION_POTENTIAL",&
            description="Parameters the interaction potential in computing the biel integrals",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="POTENTIAL_TYPE",&
            description="Which interaction potential should be used "//&
                         "(Coulomb, TShPSC operator).",&
            usage="POTENTIAL_TYPE TSHPSC",&
            enum_c_vals=s2a("COULOMB","TShPSC"),&
            enum_i_vals=(/do_hfx_potential_coulomb,&
                          do_mp2_potential_TShPSC/),&
            enum_desc=s2a("Coulomb potential: 1/r", &
                          "TShPSC: 1/x - s/Rc  ->  for x<=Rc "//&
                          "(1-s)/Rc - "//&
                          "(x-Rc)/Rc^2 + "//&
                          "(x-Rc)^2/Rc^3 - "//&
                          "(2*n^2-7*n+9-4*s)*(x-Rc)^3/(Rc^4*(n^2-2*n+1)*(n-1)) + "//&
                          "(6-3*s-4*n+n^2)*(x-Rc)^4/(Rc^5*(n^4-4*n^3+6*n^2-4*n+1))   -> for Rc<x<=n*Rc (4th order polynomial) "//&
                          " 0  -> for x>n*Rc"),&
            default_i_val=do_hfx_potential_coulomb, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TRUNCATION_RADIUS",&
            description="Determines truncation radius for the truncated TShPSC potential. "//&
                        "Only valid when doing truncated calculation",&
            usage="TRUNCATION_RADIUS 10.0",required=.FALSE.,type_of_var=real_t,&
            default_r_val=10.0_dp,&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="TShPSC_DATA",&
            description="Location of the file TShPSC.dat that contains the data for the "//&
                        "evaluation of the TShPSC G0 ",&
            usage="TShPSC_DATA /data/t_sh_p_s_c.dat",&
            default_c_val="../../t_sh_p_s_c.dat",&  ! this is tests/QS i.e. regtests are possible
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_mp2_potential

END MODULE input_cp2k_mp2

