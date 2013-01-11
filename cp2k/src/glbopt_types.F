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

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

                                    
 TYPE glbopt_config_type
   INTEGER                                    :: n_walkers
   INTEGER                                    :: n_hashes
   INTEGER                                    :: bloom_size
   INTEGER                                    :: epot_ndown
   INTEGER                                    :: epot_nup
   INTEGER                                    :: md_mincount_max
   INTEGER                                    :: hop_count_max
 END TYPE glbopt_config_type

 TYPE glbopt_mdctrl_data_type
   INTEGER                                    :: md_mincounter
   REAL(KIND=dp), DIMENSION(:), POINTER       :: epot_history
   TYPE(glbopt_config_type)                   :: cfg
 END TYPE glbopt_mdctrl_data_type

 INTEGER, PARAMETER  :: glbopt_msg_minima_new = 111
 INTEGER, PARAMETER  :: glbopt_msg_minima_old = 222
 INTEGER, PARAMETER  :: glbopt_msg_shutdown   = 999


!  TYPE glbopt_env_type
!     !PRIVATE  !TODO: make this private with getters & setters
!     INTEGER                                    :: ref_count
!     ! config parameters
!     INTEGER                                    :: n_masters
!     INTEGER                                    :: n_walkers
!     INTEGER                                    :: n_hashes
!     INTEGER                                    :: n_colvars
!     INTEGER                                    :: bloom_size
!     INTEGER                                    :: epot_ndown
!     INTEGER                                    :: epot_nup
!     INTEGER                                    :: md_mincount_max
!     INTEGER                                    :: hop_count_max
!     ! stuff passed on for master/slave communication
!     TYPE(cp_para_env_type), POINTER            :: full_para_env
!     INTEGER                                    :: my_walker_id
!     INTEGER                                    :: my_master_id
!     ! stuff needed for md_callback
!     INTEGER                                    :: md_mincounter
!     REAL(KIND=dp), DIMENSION(:), POINTER       :: epot_history
!     !mixed
!     
!  END TYPE glbopt_env_type
!
 PUBLIC :: glbopt_config_type, glbopt_mdctrl_data_type
 PUBLIC :: glbopt_msg_minima_new, glbopt_msg_minima_old, glbopt_msg_shutdown
  
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_types' 
CONTAINS
        
 
END MODULE glbopt_types

