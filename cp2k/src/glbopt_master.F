!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_master
  USE swarm_types,                     ONLY: swarm_command_type,&
                                             swarm_report_type
  USE glbopt_minhop,                   ONLY: minhop_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: angstrom,&
                                             kelvin
  USE glbopt_minhop,                   ONLY: minhop_type,minhop_steer,&
                                             minhop_init, minhop_finalize

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt'


 PUBLIC :: glbopt_master_type
 PUBLIC :: glbopt_master_init, glbopt_master_finalize
 PUBLIC :: glbopt_master_steer


 TYPE glbopt_master_type
   PRIVATE
   TYPE(minhop_type) :: minhop
 END TYPE glbopt_master_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_init(glbopt, section, n_walkers, iw, error)
    TYPE(glbopt_master_type)                 :: glbopt
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: n_walkers, iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL minhop_init(glbopt%minhop, section, n_walkers, iw, error)
 END SUBROUTINE glbopt_master_init

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_steer(glbopt, report, cmd)
    TYPE(glbopt_master_type)                 :: glbopt
    TYPE(swarm_report_type)                 :: report
    TYPE(swarm_command_type)                :: cmd

    CALL minhop_steer(glbopt%minhop, report, cmd)
 END SUBROUTINE glbopt_master_steer


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_finalize(glbopt)
    TYPE(glbopt_master_type)                 :: glbopt

    CALL minhop_finalize(glbopt%minhop)
 END SUBROUTINE glbopt_master_finalize


END MODULE glbopt_master

