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
  USE swarm_message,                   ONLY: swarm_message_type
  USE kinds,                           ONLY: dp, int_8

 IMPLICIT NONE
 PRIVATE

  TYPE history_type
    PRIVATE
    TYPE(history_entry_type), POINTER :: root      => Null()
    REAL(KIND=dp)                     :: PRECISION =  1.0e-6
    INTEGER                           :: iw = -1
    !.. config params
  END TYPE history_type

  ! a linked list
  TYPE history_entry_type
    REAL(KIND=dp)                     :: Epot = HUGE(1.0_dp)
    TYPE(history_entry_type), POINTER :: next => Null()
  END TYPE history_entry_type


  PUBLIC :: history_type
  PUBLIC :: history_init, history_finalize
  PUBLIC :: history_add, history_lookup

  CONTAINS

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE history_init(history, iw)
    TYPE(history_type)                :: history
    INTEGER                                  :: iw
    history%iw = iw

!TYPE(section_vals_type), POINTER         :: history_section
!TYPE(cp_error_type), INTENT(inout)       :: error
!nothing todo, yet

 END SUBROUTINE history_init

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE history_finalize(history)
    TYPE(history_type)                :: history

!nothing todo, yet
!TODO: deallocate 

 END SUBROUTINE history_finalize


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE history_add(history, Epot)
    TYPE(history_type)                :: history
    REAL(KIND=dp)                     :: Epot


    TYPE(history_entry_type), POINTER        :: last_entry, new_entry

    ALLOCATE(new_entry)
    new_entry%Epot = Epot

    IF(.NOT. ASSOCIATED(history%root)) THEN
       history%root => new_entry
    ELSE
      last_entry => history%root
      DO WHILE(ASSOCIATED(last_entry%next))
        last_entry => last_entry%next
      END DO
      last_entry%next => new_entry
    END IF
 END SUBROUTINE history_add


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE history_lookup(history, Epot, found)
    TYPE(history_type), INTENT(IN)    :: history
    REAL(KIND=dp)                     :: Epot
    LOGICAL, INTENT(OUT)              :: found

    TYPE(history_entry_type), POINTER        :: entry
    INTEGER(kind=int_8)                      :: a, b
    REAL(kind=dp)                            :: closest

    closest = HUGE(1.0)
    ! discretize accoring to history%precision
    a = FLOOR(Epot / history%precision, kind=int_8)

    found = .FALSE.
    ENTRY => history%root
    DO WHILE(ASSOCIATED(ENTRY))
     closest = MIN(closest, ABS(entry%Epot - Epot))
      !TODO: could store discretized energies instead of reports in a growing array
      b = FLOOR(entry%Epot / history%precision, kind=int_8)
      IF(a == b) THEN
         found = .TRUE.
         !IF(history%iw>0) &
         !  WRITE(history%iw, *) "history: Found: ",report%Epot, "dist: ", closest
         RETURN
      END IF
      ENTRY => entry%next
    END DO

    !IF(history%iw>0) &
    !  WRITE(history%iw, *) "history: Did not find: ",report%Epot, "closest dist: ", closest
 END SUBROUTINE history_lookup


END MODULE glbopt_history
