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
   !INTEGER                                             :: n_workers = -1
   !TYPE(cp_error_type)                                 :: error
   TYPE(minhop_type) :: minhop
 END TYPE glbopt_master_type

 CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_init(glbopt, swarm_section, n_walkers, iw, error)
    TYPE(glbopt_master_type)                 :: glbopt
    TYPE(section_vals_type), POINTER         :: swarm_section
    INTEGER, INTENT(IN)                      :: n_walkers, iw
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(section_vals_type), POINTER         :: glbopt_section

    glbopt%iw = iw

    glbopt_section => section_vals_get_subs_vals(swarm_section, "GLOBAL_OPT", error=error)
    CALL section_vals_val_get(glbopt_section,"E_MIN", r_val=glbopt%Emin, error=error)

    CALL minhop_init(glbopt%minhop, glbopt_section, n_walkers, iw, error)
 END SUBROUTINE glbopt_master_init

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_steer(glbopt, report, cmd, should_stop)
    TYPE(glbopt_master_type)                 :: glbopt
    TYPE(swarm_message_type)                 :: report
    TYPE(swarm_message_type)                 :: cmd
    LOGICAL, INTENT(INOUT)                   :: should_stop

    REAL(KIND=dp) :: report_Epot
    INTEGER :: report_worker_id

    CALL swarm_message_get(report, "worker_id", report_worker_id)
    CALL swarm_message_get(report, "Epot", report_Epot)

    glbopt%Elowest = MIN(glbopt%Elowest, report_Epot)
    IF(glbopt%iw>0) THEN
       WRITE (glbopt%iw,'(A,46X,I8)')&
              " GLBOPT| Reporting worker ", report_worker_id
       WRITE (glbopt%iw,'(A,20X,E15.8)')&
              " GLBOPT| Reported potential Energy [Hartree] ", report_Epot
       WRITE (glbopt%iw,'(A,13X,E15.8)')&
        " GLBOPT| Lowest reported potential Energy [Hartree] ",glbopt%Elowest
    END IF

    !IF(glbopt%Elowest<glbopt%Emin .AND. .NOT. master%should_stop) THEN
    IF(glbopt%Elowest<glbopt%Emin) THEN
       IF(glbopt%iw>0)  WRITE (glbopt%iw, *)  "GLBOPT| Reached Emin > Epot. Quitting."
       should_stop = .TRUE.
    ELSE
       CALL minhop_steer(glbopt%minhop, report, cmd)
    ENDIF
 END SUBROUTINE glbopt_master_steer


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_finalize(glbopt)
    TYPE(glbopt_master_type)                 :: glbopt

    CALL minhop_finalize(glbopt%minhop)
 END SUBROUTINE glbopt_master_finalize


END MODULE glbopt_master

