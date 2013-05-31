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
  USE swarm_message,                   ONLY: swarm_message_type

  USE glbopt_minhop,                   ONLY: minhop_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: angstrom,&
                                             kelvin
  USE glbopt_minhop,                   ONLY: minhop_type,minhop_steer,&
                                             minhop_init, minhop_finalize
  USE swarm_message,                   ONLY: swarm_message_type,&
                                             swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_haskey,&
                                             swarm_message_length
  USE atoms_input,                     ONLY: read_atoms_input
  USE f77_interface,                   ONLY: create_force_env,&
                                             destroy_force_env,&
                                             f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type,&
                                             pack_subsys_particles
  USE f77_interface,                   ONLY: create_force_env,&
                                             destroy_force_env,&
                                             f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE geo_opt,                         ONLY: cp_geo_opt
  USE glbopt_types,                    ONLY: glbopt_mdctrl_data_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE md_run,                          ONLY: qs_mol_dyn
  USE mdctrl_types,                    ONLY: mdctrl_type
  USE parallel_rng_types,              ONLY: reset_to_next_rng_substream
  USE physcon,                         ONLY: angstrom,&
                                             femtoseconds,&
                                             kelvin
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE swarm_message,                   ONLY: swarm_message_type,&
                                             swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_length

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt'


 PUBLIC :: glbopt_master_type
 PUBLIC :: glbopt_master_init, glbopt_master_finalize
 PUBLIC :: glbopt_master_steer


 TYPE glbopt_master_type
   PRIVATE
   REAL(KIND=dp)                                       :: Elowest = HUGE(1.0_dp)
   REAL(KIND=dp)                                       :: Emin = TINY(1.0_dp)
   INTEGER                                             :: iw = 0
   !REAL(KIND=dp), DIMENSION(:), POINTER                :: start_positions => Null()
   !INTEGER                                             :: n_workers = -1
   !TYPE(cp_error_type)                                 :: error
   TYPE(minhop_type) :: minhop

 END TYPE glbopt_master_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_init(this, para_env, swarm_section, input_path, n_walkers, iw, error)
    TYPE(glbopt_master_type)                 :: this
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: swarm_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: n_walkers, iw
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(section_vals_type), POINTER         :: glbopt_section


    this%iw = iw

    glbopt_section => section_vals_get_subs_vals(swarm_section, "GLOBAL_OPT", error=error)
    CALL section_vals_val_get(glbopt_section,"E_MIN", r_val=this%Emin, error=error)

    CALL minhop_init(this%minhop, glbopt_section, n_walkers, iw, error)
 END SUBROUTINE glbopt_master_init


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_steer(this, report, cmd, should_stop)
    TYPE(glbopt_master_type)                 :: this
    TYPE(swarm_message_type)                 :: report
    TYPE(swarm_message_type)                 :: cmd
    LOGICAL, INTENT(INOUT)                   :: should_stop

    REAL(KIND=dp) :: report_Epot
    INTEGER :: report_worker_id

    CALL swarm_message_get(report, "worker_id", report_worker_id)

    IF(.NOT. swarm_message_haskey(report, "initial_hello")) THEN
       CALL swarm_message_get(report, "Epot", report_Epot)
       this%Elowest = MIN(this%Elowest, report_Epot)
       IF(this%iw>0) THEN
          WRITE (this%iw,'(A,46X,I8)')&
                 " GLBOPT| Reporting worker ", report_worker_id
          WRITE (this%iw,'(A,20X,E15.8)')&
                 " GLBOPT| Reported potential Energy [Hartree] ", report_Epot
          WRITE (this%iw,'(A,13X,E15.8)')&
           " GLBOPT| Lowest reported potential Energy [Hartree] ", this%Elowest
       END IF
    END IF

    IF(this%Elowest<this%Emin) THEN
       IF(this%iw>0)  WRITE (this%iw, *)  "GLBOPT| Reached Emin > Epot. Quitting."
       should_stop = .TRUE.
    ELSE
       CALL minhop_steer(this%minhop, report, cmd)
    ENDIF
 END SUBROUTINE glbopt_master_steer


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_finalize(this)
    TYPE(glbopt_master_type)                 :: this

    CALL minhop_finalize(this%minhop)
 END SUBROUTINE glbopt_master_finalize


END MODULE glbopt_master

