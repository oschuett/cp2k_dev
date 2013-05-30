!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE swarm_types

  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

 TYPE swarm_command_type
   INTEGER :: worker_id = -1
   INTEGER :: cmd_id = -1
   INTEGER :: iframe = -1
   REAL(KIND=dp) :: temperature = 0.0
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: positions
 END TYPE swarm_command_type


 TYPE swarm_report_type
   INTEGER                                   :: worker_id = -1
   INTEGER                                   :: iframe = -1
   REAL(KIND=dp)                             :: Epot = 0.0
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: positions
 END TYPE swarm_report_type


 INTEGER, PARAMETER  :: SWARM_CMD_SHUTDOWN      = 999


 PUBLIC :: SWARM_CMD_SHUTDOWN
 PUBLIC :: swarm_report_type, swarm_command_type

END MODULE swarm_types


