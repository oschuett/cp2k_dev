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
   REAL                                       :: beta1
   REAL                                       :: beta2
   REAL                                       :: beta3
   REAL                                       :: Ediff0
   REAL                                       :: alpha1
   REAL                                       :: alpha2     
   
  
!19	    _default_settings = {
!20	        'T0': 1000.,  # K, initial MD 'temperature'
!21	        'beta1': 1.1,  # temperature adjustment parameter
!22	        'beta2': 1.1,  # temperature adjustment parameter
!23	        'beta3': 1. / 1.1,  # temperature adjustment parameter
!24	        'Ediff0': 0.5,  # eV, initial energy acceptance threshold
!25	        'alpha1': 0.98,  # energy threshold adjustment parameter
!26	        'alpha2': 1. / 0.98,  # energy threshold adjustment parameter
!27	        'mdmin': 2,  # criteria to stop MD simulation (no. of minima)
!28	        'logfile': 'hop.log',  # text log
!29	        'minima_threshold': 0.5,  # A, threshold for identical configs
!30	        'timestep': 1.0,  # fs, timestep for MD simulations
!31	        'optimizer': QuasiNewton,  # local optimizer to use
!32	        'minima_traj': 'minima.traj',  # storage file for minima list
!33	        'fmax' : 0.05, # eV/A, max force for optimizations
!34	                          }
   
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

