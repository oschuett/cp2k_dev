!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_history
  USE glbopt_types,                    ONLY: glbopt_report_type
  USE kinds,                           ONLY: dp, int_8

 IMPLICIT NONE
 PRIVATE

  TYPE glbopt_history_type
    PRIVATE
    TYPE(history_entry_type), POINTER :: root      => Null()
    REAL(KIND=dp)                     :: PRECISION =  1.0e-6
    INTEGER                           :: iw = -1
    !.. config params
  END TYPE glbopt_history_type
  
  ! a linked list
  TYPE history_entry_type
    TYPE(glbopt_report_type) :: report
    TYPE(history_entry_type), POINTER :: next => Null()
  END TYPE history_entry_type


  PUBLIC :: glbopt_history_type
  PUBLIC :: glbopt_history_init, glbopt_history_finalize
  PUBLIC :: glbopt_history_add, glbopt_history_lookup

  CONTAINS

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_history_init(history, iw)
    TYPE(glbopt_history_type)                :: history
    INTEGER                                  :: iw
    history%iw = iw

!TYPE(section_vals_type), POINTER         :: history_section
!TYPE(cp_error_type), INTENT(inout)       :: error
!nothing todo, yet

 END SUBROUTINE glbopt_history_init

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_history_finalize(history)
    TYPE(glbopt_history_type)                :: history

!nothing todo, yet
!TODO: deallocate 

 END SUBROUTINE glbopt_history_finalize


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_history_add(history, report)
    TYPE(glbopt_history_type)                :: history
    TYPE(glbopt_report_type)                 :: report

    TYPE(history_entry_type), POINTER        :: last_entry, new_entry

    ALLOCATE(new_entry)
    new_entry%report = report

    IF(.NOT. ASSOCIATED(history%root)) THEN
       history%root => new_entry
    ELSE
      last_entry => history%root
      DO WHILE(ASSOCIATED(last_entry%next))
        last_entry => last_entry%next
      END DO
      last_entry%next => new_entry
    END IF
 END SUBROUTINE glbopt_history_add


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_history_lookup(history, report, found)
    TYPE(glbopt_history_type), INTENT(IN)    :: history
    TYPE(glbopt_report_type), INTENT(IN)     :: report
    LOGICAL, INTENT(OUT)                     :: found

    TYPE(history_entry_type), POINTER        :: ENTRY
    INTEGER(kind=int_8)                       :: a, b

    ! discretize accoring to history%precision
    a = INT(report%Epot / history%precision, kind=int_8)

    found = .FALSE.
    ENTRY => history%root
    DO WHILE(ASSOCIATED(ENTRY))
      !TODO: could store discretized energies instead of reports in a growing array
      b = INT(entry%report%Epot / history%precision, kind=int_8)
      IF(a == b) THEN
         found = .TRUE.
         RETURN
      END IF
      ENTRY => entry%next
    END DO

 END SUBROUTINE glbopt_history_lookup


END MODULE glbopt_history
