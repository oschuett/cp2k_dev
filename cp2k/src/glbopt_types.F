!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_types

  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE


 TYPE glbopt_mdctrl_data_type
   INTEGER                                    :: md_mincounter
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE   :: epot_history
   INTEGER                                    :: output_unit
   INTEGER                                    :: itimes
   INTEGER                                    :: epot_nup
   INTEGER                                    :: epot_ndown
   INTEGER                                    :: md_mincount_max
   REAL(KIND=dp)                              :: md_max_dist
 END TYPE glbopt_mdctrl_data_type


 TYPE glbopt_command_type
   INTEGER :: walker_id = -1
   INTEGER :: cmd_id = -1
   INTEGER :: iframe = -1
   REAL(KIND=dp) :: temperature = 0.0
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: positions
 END TYPE glbopt_command_type


 TYPE glbopt_report_type
   INTEGER                                   :: walker_id = -1
   INTEGER                                   :: iframe = -1
   REAL(KIND=dp)                             :: Epot = 0.0
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: positions
 END TYPE glbopt_report_type


 INTEGER, PARAMETER  :: GLBOPT_CMD_SHUTDOWN      = 999
 INTEGER, PARAMETER  :: GLBOPT_CMD_MDGOPT        = 888
 INTEGER, PARAMETER  :: GLBOPT_CMD_SETIFRAME     = 777


 PUBLIC :: GLBOPT_CMD_SHUTDOWN, GLBOPT_CMD_MDGOPT, GLBOPT_CMD_SETIFRAME
 PUBLIC :: glbopt_report_type, glbopt_command_type, glbopt_mdctrl_data_type

END MODULE glbopt_types


