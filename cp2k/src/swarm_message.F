!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE swarm_message

  USE kinds,                           ONLY: int_4, int_8, real_4, real_8, default_string_length
  USE message_passing,                 ONLY: &
       mp_abort, mp_any_source, mp_bcast, mp_comm_free, mp_comm_split, &
       mp_comm_split_direct, mp_environ, mp_recv, mp_send, mp_sum, mp_sync

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

 TYPE swarm_message_type
   PRIVATE
   TYPE(message_entry_type), POINTER :: root => Null()
 END TYPE swarm_message_type

 INTEGER, PARAMETER  :: key_length      = 20
 
 TYPE message_entry_type
   CHARACTER(LEN=key_length)                       :: key
   TYPE(message_entry_type),            POINTER    :: next        => Null()
   CHARACTER(LEN=default_string_length),POINTER    :: value_str   => Null()
   INTEGER(KIND=int_4),                 POINTER    :: value_i4    => Null()
   INTEGER(KIND=int_8),                 POINTER    :: value_i8    => Null()
   REAL(KIND=real_4),                   POINTER    :: value_r4    => Null()
   REAL(KIND=real_8),                   POINTER    :: value_r8    => Null()
   INTEGER(KIND=int_4), DIMENSION(:),   POINTER    :: value_i4_1d => Null()
   INTEGER(KIND=int_8), DIMENSION(:),   POINTER    :: value_i8_1d => Null()
   REAL(KIND=real_4),   DIMENSION(:),   POINTER    :: value_r4_1d => Null()
   REAL(KIND=real_8),   DIMENSION(:),   POINTER    :: value_r8_1d => Null()
 END TYPE message_entry_type

  INTERFACE swarm_message_add
     MODULE PROCEDURE swarm_message_add_str
     MODULE PROCEDURE swarm_message_add_i4,    swarm_message_add_i8
     MODULE PROCEDURE swarm_message_add_r4,    swarm_message_add_r8
     MODULE PROCEDURE swarm_message_add_i4_1d, swarm_message_add_i8_1d
     MODULE PROCEDURE swarm_message_add_r4_1d, swarm_message_add_r8_1d
  END INTERFACE swarm_message_add

  INTERFACE swarm_message_get
     MODULE PROCEDURE swarm_message_get_str
     MODULE PROCEDURE swarm_message_get_i4,    swarm_message_get_i8
     MODULE PROCEDURE swarm_message_get_r4,    swarm_message_get_r8
     MODULE PROCEDURE swarm_message_get_i4_1d, swarm_message_get_i8_1d
     MODULE PROCEDURE swarm_message_get_r4_1d, swarm_message_get_r8_1d
  END INTERFACE swarm_message_get

 PUBLIC :: swarm_message_type, swarm_message_add, swarm_message_get
 PUBLIC :: swarm_message_mpi_send, swarm_message_mpi_recv, swarm_message_mpi_bcast
 PUBLIC :: swarm_message_length, swarm_message_haskey

 CONTAINS


 ! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_bcast(entry, src, group, mepos)
    TYPE(message_entry_type),INTENT(INOUT) :: entry
    INTEGER, INTENT(IN)                    :: src, group, mepos

    INTEGER, DIMENSION(key_length) :: key_arr
    INTEGER                        :: datatype, datasize
    INTEGER, DIMENSION(:), ALLOCATABLE :: value_str_arr
    if(src==mepos) key_arr = TRANSFER(entry%key, key_arr)

    CALL mp_bcast(key_arr, src, group)
    if(src/=mepos) entry%key = TRANSFER(key_arr, entry%key)

    IF(src==mepos) THEN
       datasize = 1
       IF(ASSOCIATED(entry%value_i4)) THEN
          datatype = 1
       ELSE IF(ASSOCIATED(entry%value_i8)) THEN
          datatype = 2
       ELSE IF(ASSOCIATED(entry%value_r4)) THEN
          datatype = 3
       ELSE IF(ASSOCIATED(entry%value_r8)) THEN
          datatype = 4
       ELSE IF(ASSOCIATED(entry%value_i4_1d)) THEN
          datatype = 5
          datasize = SIZE(entry%value_i4_1d)
       ELSE IF(ASSOCIATED(entry%value_i8_1d)) THEN
          datatype = 6
          datasize = SIZE(entry%value_i8_1d)
       ELSE IF(ASSOCIATED(entry%value_r4_1d)) THEN
          datatype = 7
          datasize = SIZE(entry%value_r4_1d)
       ELSE IF(ASSOCIATED(entry%value_r8_1d)) THEN
          datatype = 8
          datasize = SIZE(entry%value_r8_1d)
       ELSE IF(ASSOCIATED(entry%value_str)) THEN
          datatype = 9
       ELSE
          STOP "swarm_message_entry_mpi_bcast: no value associated"
       END IF
    END IF
    CALL mp_bcast(datatype, src, group)
    CALL mp_bcast(datasize, src, group)

    SELECT CASE(datatype)
    CASE(1)
       IF(src/=mepos) ALLOCATE(entry%value_i4)
       CALL mp_bcast(entry%value_i4, src, group)
    CASE(2)
       IF(src/=mepos) ALLOCATE(entry%value_i8)
       CALL mp_bcast(entry%value_i8, src, group)
    CASE(3)
       IF(src/=mepos) ALLOCATE(entry%value_r4)
       CALL mp_bcast(entry%value_r4, src, group)
    CASE(4)
       IF(src/=mepos) ALLOCATE(entry%value_r8)
       CALL mp_bcast(entry%value_r8, src, group)
    CASE(5)
       IF(src/=mepos) ALLOCATE(entry%value_i4_1d(datasize))
       CALL mp_bcast(entry%value_i4_1d, src, group)
    CASE(6)
       IF(src/=mepos) ALLOCATE(entry%value_i8_1d(datasize))
       CALL mp_bcast(entry%value_i8_1d, src, group)
    CASE(7)
       IF(src/=mepos) ALLOCATE(entry%value_r4_1d(datasize))
       CALL mp_bcast(entry%value_r4_1d, src, group)
    CASE(8)
       IF(src/=mepos) ALLOCATE(entry%value_r8_1d(datasize))
       CALL mp_bcast(entry%value_r8_1d, src, group)
    CASE(9)
       ALLOCATE(value_str_arr(default_string_length))
       IF(src==mepos) value_str_arr=TRANSFER(entry%value_str, value_str_arr)
       CALL mp_bcast(value_str_arr, src, group)
       IF(src/=mepos) THEN
          ALLOCATE(entry%value_str)
          entry%value_str = TRANSFER(value_str_arr, entry%value_str)
       END IF
    CASE DEFAULT
       STOP "swarm_message_entry_mpi_bcast: unkown datatype"
    END SELECT

  END SUBROUTINE swarm_message_entry_mpi_bcast

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_send(entry, group, dest, tag)
    TYPE(message_entry_type),INTENT(IN) :: entry
    INTEGER, INTENT(IN)                 :: group, dest, tag

    INTEGER, DIMENSION(key_length) :: key_arr
    INTEGER, DIMENSION(:), ALLOCATABLE :: value_str_arr

    key_arr = TRANSFER(entry%key, key_arr)
    CALL mp_send(key_arr, dest, tag, group)

    IF(ASSOCIATED(entry%value_i4)) THEN
       CALL mp_send(1, dest, tag, group)
       CALL mp_send(entry%value_i4, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_i8)) THEN
       CALL mp_send(2, dest, tag, group)
       CALL mp_send(entry%value_i8, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_r4)) THEN
       CALL mp_send(3, dest, tag, group)
       CALL mp_send(entry%value_r4, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_r8)) THEN
       CALL mp_send(4, dest, tag, group)
       CALL mp_send(entry%value_r8, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_i4_1d)) THEN
       CALL mp_send(5, dest, tag, group)
       CALL mp_send(SIZE(entry%value_i4_1d), dest, tag, group)
       CALL mp_send(entry%value_i4_1d, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_i8_1d)) THEN
       CALL mp_send(6, dest, tag, group)
       CALL mp_send(SIZE(entry%value_i8_1d), dest, tag, group)
       CALL mp_send(entry%value_i8_1d, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_r4_1d)) THEN
       CALL mp_send(7, dest, tag, group)
       CALL mp_send(SIZE(entry%value_r4_1d), dest, tag, group)
       CALL mp_send(entry%value_r4_1d, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_r8_1d)) THEN
       CALL mp_send(8, dest, tag, group)
       CALL mp_send(SIZE(entry%value_r8_1d), dest, tag, group)
       CALL mp_send(entry%value_r8_1d, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_str)) THEN
       CALL mp_send(9, dest, tag, group)
       ALLOCATE(value_str_arr(default_string_length))
       value_str_arr = TRANSFER(entry%value_str, value_str_arr)
       CALL mp_send(value_str_arr, dest, tag, group)
    ELSE
       CALL abort()
       STOP "swarm_message_entry_mpi_send: no value associated"
    END IF
  END SUBROUTINE swarm_message_entry_mpi_send


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_recv(entry, group, src, tag)
    TYPE(message_entry_type),INTENT(INOUT) :: entry
    INTEGER, INTENT(IN)                    :: group
    INTEGER, INTENT(INOUT)                 :: src, tag

    INTEGER, DIMENSION(key_length) :: key_arr
    INTEGER, DIMENSION(:), ALLOCATABLE :: value_str_arr
    INTEGER :: datatype, s

    CALL mp_recv(key_arr, src, tag, group)
    entry%key = TRANSFER(key_arr, entry%key)

    CALL mp_recv(datatype, src, tag, group)

    SELECT CASE(datatype)
    CASE(1)
       ALLOCATE(entry%value_i4)
       CALL mp_recv(entry%value_i4, src, tag, group)
    CASE(2)
       ALLOCATE(entry%value_i8)
       CALL mp_recv(entry%value_i8, src, tag, group)
    CASE(3)
       ALLOCATE(entry%value_r4)
       CALL mp_recv(entry%value_r4, src, tag, group)
    CASE(4)
       ALLOCATE(entry%value_r8)
       CALL mp_recv(entry%value_r8, src, tag, group)

    CASE(5)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_i4_1d(s))
       CALL mp_recv(entry%value_i4_1d, src, tag, group)
    CASE(6)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_i8_1d(s))
       CALL mp_recv(entry%value_i8_1d, src, tag, group)
    CASE(7)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_r4_1d(s))
       CALL mp_recv(entry%value_r4_1d, src, tag, group)
    CASE(8)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_r8_1d(s))
       CALL mp_recv(entry%value_r8_1d, src, tag, group)
    CASE(9)
       ALLOCATE(value_str_arr(default_string_length))
       ALLOCATE(entry%value_str)
       CALL mp_recv(value_str_arr, src, tag, group)
       entry%value_str = TRANSFER(value_str_arr, entry%value_str)
    CASE DEFAULT
       STOP "swarm_message_entry_mpi_recv: unkown datatype"
    END SELECT
  END SUBROUTINE swarm_message_entry_mpi_recv


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_mpi_send(msg, group, dest, tag)
    TYPE(swarm_message_type),INTENT(IN) :: msg
    INTEGER, INTENT(IN)                 :: group, dest, tag

    TYPE(message_entry_type), POINTER :: curr_entry
    CALL mp_send(swarm_message_length(msg), dest, tag, group)
    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      CALL swarm_message_entry_mpi_send(curr_entry, group, dest, tag)
      curr_entry => curr_entry%next
    END DO
  END SUBROUTINE swarm_message_mpi_send

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_mpi_recv(msg, group, src, tag)
    TYPE(swarm_message_type),INTENT(INOUT) :: msg
    INTEGER, INTENT(IN)                    :: group
    INTEGER, INTENT(INOUT)                 :: src, tag

    TYPE(message_entry_type), POINTER      :: new_entry
    INTEGER :: length, i

    IF(ASSOCIATED(msg%root)) STOP "swarm_message_mpi_recv: message not empty"
    CALL mp_recv(length, src, tag, group)
    DO i=1, length
       ALLOCATE(new_entry)
       CALL swarm_message_entry_mpi_recv(new_entry, group, src, tag)
       new_entry%next => msg%root
       msg%root => new_entry
    END DO

  END SUBROUTINE swarm_message_mpi_recv

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_mpi_bcast(msg, src, group)
    TYPE(swarm_message_type),INTENT(INOUT) :: msg
    INTEGER, INTENT(IN)                    :: src, group

    INTEGER                                :: num_pe, mepos
    TYPE(message_entry_type), POINTER      :: curr_entry
    INTEGER :: length, i

    CALL mp_environ(num_pe, mepos, group)

    IF(mepos/=src .AND. ASSOCIATED(msg%root)) STOP "swarm_message_mpi_bcast: message not empty"
    length = swarm_message_length(msg)
    CALL mp_bcast(length, src, group)

    IF(mepos==src) curr_entry => msg%root

    DO i=1, length
       IF(mepos/=src) ALLOCATE(curr_entry)

       CALL swarm_message_entry_mpi_bcast(curr_entry, src, group, mepos)

       IF(mepos==src) THEN
          curr_entry=>curr_entry%next
       ELSE
          curr_entry%next => msg%root
          msg%root => curr_entry
       END IF
    END DO

  END SUBROUTINE swarm_message_mpi_bcast


! *****************************************************************************
! *****************************************************************************
  FUNCTION swarm_message_length(msg) RESULT(l)
    INTEGER                              :: l
    TYPE(swarm_message_type), INTENT(IN) :: msg
    TYPE(message_entry_type), POINTER    :: curr_entry
    l = 0
    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      l = l + 1
      curr_entry => curr_entry%next
    END DO
  END FUNCTION swarm_message_length


! *****************************************************************************
! *****************************************************************************
  FUNCTION swarm_message_haskey(msg, key) RESULT(res)
    TYPE(swarm_message_type), INTENT(IN) :: msg
    CHARACTER(LEN=*), INTENT(IN)         :: key
    LOGICAL                              :: res

    TYPE(message_entry_type), POINTER    :: curr_entry
    res = .FALSE.
    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         res = .TRUE.
         EXIT
      END IF
      curr_entry => curr_entry%next
    END DO
  END FUNCTION swarm_message_haskey
! *****************************************************************************
! *****************************************************************************

!  END SUBROUTINE swarm_mpi_recv_message_entry(group, dest, tag, entry)
  
! *****************************************************************************
! *****************************************************************************

#include "swarm_message_str.F"
#include "swarm_message_i4.F"
#include "swarm_message_i8.F"
#include "swarm_message_r4.F"
#include "swarm_message_r8.F"
#include "swarm_message_i4_1d.F"
#include "swarm_message_i8_1d.F"
#include "swarm_message_r4_1d.F"
#include "swarm_message_r8_1d.F"


END MODULE swarm_message


