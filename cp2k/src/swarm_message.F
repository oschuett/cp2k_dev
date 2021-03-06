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

  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object
  USE cp_parser_types,                 ONLY: cp_parser_type
  USE kinds,                           ONLY: default_string_length

#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

 TYPE swarm_message_type
   PRIVATE
   TYPE(message_entry_type), POINTER :: root => Null()
 END TYPE swarm_message_type

 INTEGER, PARAMETER  :: key_length      = 20

 TYPE message_entry_type
   CHARACTER(LEN=key_length)                        :: key
   TYPE(message_entry_type),            POINTER     :: next        => Null()
   CHARACTER(LEN=default_string_length),POINTER     :: value_str   => Null()
   INTEGER(KIND=int_4),                 POINTER     :: value_i4    => Null()
   INTEGER(KIND=int_8),                 POINTER     :: value_i8    => Null()
   REAL(KIND=real_4),                   POINTER     :: value_r4    => Null()
   REAL(KIND=real_8),                   POINTER     :: value_r8    => Null()
   INTEGER(KIND=int_4), DIMENSION(:),   POINTER     :: value_i4_1d => Null()
   INTEGER(KIND=int_8), DIMENSION(:),   POINTER     :: value_i8_1d => Null()
   REAL(KIND=real_4),   DIMENSION(:),   POINTER     :: value_r4_1d => Null()
   REAL(KIND=real_8),   DIMENSION(:),   POINTER     :: value_r8_1d => Null()
 END TYPE message_entry_type

 INTERFACE swarm_message_add
    MODULE PROCEDURE swarm_message_add_str
    MODULE PROCEDURE swarm_message_add_i4,    swarm_message_add_i8
    MODULE PROCEDURE swarm_message_add_r4,    swarm_message_add_r8
    MODULE PROCEDURE swarm_message_add_i4_1d, swarm_message_add_i8_1d
    MODULE PROCEDURE swarm_message_add_r4_1d, swarm_message_add_r8_1d
 END INTERFACE swarm_message_add

 INTERFACE swarm_message_get
    MODULE PROCEDURE swarm_message_get_str_subrt
    MODULE PROCEDURE swarm_message_get_i4_subrt,    swarm_message_get_i8_subrt
    MODULE PROCEDURE swarm_message_get_r4_subrt,    swarm_message_get_r8_subrt
    MODULE PROCEDURE swarm_message_get_i4_1d_subrt, swarm_message_get_i8_1d_subrt
    MODULE PROCEDURE swarm_message_get_r4_1d_subrt, swarm_message_get_r8_1d_subrt
 END INTERFACE swarm_message_get

 PUBLIC :: swarm_message_get_str
 PUBLIC :: swarm_message_get_i4,    swarm_message_get_i8
 PUBLIC :: swarm_message_get_r4,    swarm_message_get_r8
 PUBLIC :: swarm_message_get_i4_1d, swarm_message_get_i8_1d
 PUBLIC :: swarm_message_get_r4_1d, swarm_message_get_r8_1d

 PUBLIC :: swarm_message_type, swarm_message_add, swarm_message_get
 PUBLIC :: swarm_message_mpi_send, swarm_message_mpi_recv, swarm_message_mpi_bcast
 PUBLIC :: swarm_message_file_write, swarm_message_file_read
 PUBLIC :: swarm_message_length, swarm_message_haskey, swarm_message_equal
 PUBLIC :: swarm_message_free

 CONTAINS


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_free(msg)
     TYPE(swarm_message_type),INTENT(INOUT) :: msg
     TYPE(message_entry_type), POINTER :: entry, old_entry

     entry => msg%root
     DO WHILE(ASSOCIATED(entry))
       IF(ASSOCIATED(entry%value_str))      DEALLOCATE(entry%value_str)
       IF(ASSOCIATED(entry%value_i4))       DEALLOCATE(entry%value_i4)
       IF(ASSOCIATED(entry%value_i8))       DEALLOCATE(entry%value_i8)
       IF(ASSOCIATED(entry%value_r4))       DEALLOCATE(entry%value_r4)
       IF(ASSOCIATED(entry%value_r8))       DEALLOCATE(entry%value_r8)
       IF(ASSOCIATED(entry%value_i4_1d))    DEALLOCATE(entry%value_i4_1d)
       IF(ASSOCIATED(entry%value_i8_1d))    DEALLOCATE(entry%value_i8_1d)
       IF(ASSOCIATED(entry%value_r4_1d))    DEALLOCATE(entry%value_r4_1d)
       IF(ASSOCIATED(entry%value_r8_1d))    DEALLOCATE(entry%value_r8_1d)
       old_entry => entry
       entry => entry%next
       DEALLOCATE(old_entry)
     END DO

     NULLIFY(msg%root)

     IF(swarm_message_length(msg) > 0) STOP "swarm_message_free failed"
  END SUBROUTINE swarm_message_free

! *****************************************************************************
! *****************************************************************************
  FUNCTION swarm_message_equal(msg1, msg2) result(res)
     TYPE(swarm_message_type),INTENT(IN) :: msg1, msg2
     LOGICAL                             :: res

     res =  swarm_message_equal_oneway(msg1, msg2) .AND. &
            swarm_message_equal_oneway(msg2, msg1)

  END FUNCTION swarm_message_equal

! *****************************************************************************
! *****************************************************************************
  FUNCTION swarm_message_equal_oneway(msg1, msg2) result(res)
     TYPE(swarm_message_type),INTENT(IN) :: msg1, msg2
     TYPE(message_entry_type), POINTER   :: entry1, entry2
     LOGICAL                             :: res
     LOGICAL                             :: found

     res = .FALSE.

     !loop over entries of msg1
     entry1 => msg1%root
     DO WHILE(ASSOCIATED(entry1))

        ! finding matching entry in msg2
        entry2 => msg2%root
        found = .FALSE.
        DO WHILE(ASSOCIATED(entry2))
           IF(TRIM(entry2%key) == TRIM(entry1%key)) THEN
              found = .TRUE.
              EXIT
           END IF
           entry2 => entry2%next
        END DO
        IF(.NOT. found) RETURN

        !compare the two entries
        IF(ASSOCIATED(entry1%value_str)) THEN
           IF(.NOT.ASSOCIATED(entry2%value_str)) RETURN
           If(TRIM(entry1%value_str) /= TRIM(entry2%value_str)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_i4)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_i4)) RETURN
           IF(entry1%value_i4 /= entry2%value_i4) RETURN

        ELSE IF(ASSOCIATED(entry1%value_i8)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_i8)) RETURN
           IF(entry1%value_i8 /= entry2%value_i8) RETURN

        ELSE IF(ASSOCIATED(entry1%value_r4)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_r4)) RETURN
           IF(ABS(entry1%value_r4-entry2%value_r4)>1e-5) RETURN

        ELSE IF(ASSOCIATED(entry1%value_r8)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_r8)) RETURN
           IF(ABS(entry1%value_r8-entry2%value_r8)>1e-10) RETURN

        ELSE IF(ASSOCIATED(entry1%value_i4_1d)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_i4_1d)) RETURN
           IF(ANY(entry1%value_i4_1d /= entry2%value_i4_1d)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_i8_1d)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_i8_1d)) RETURN
           IF(ANY(entry1%value_i8_1d /= entry2%value_i8_1d)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_r4_1d)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_r4_1d)) RETURN
           IF(ANY(ABS(entry1%value_r4_1d-entry2%value_r4_1d)>1e-5)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_r8_1d)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_r8_1d)) RETURN
           IF(ANY(ABS(entry1%value_r8_1d-entry2%value_r8_1d)>1e-10)) RETURN
        ELSE
           STOP "swarm_message_equal_oneway: no value ASSOCIATED"
        END IF

        entry1 => entry1%next
     END DO

     ! if we reach this point no differences were found
     res = .TRUE.
  END FUNCTION swarm_message_equal_oneway


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_bcast(entry, src, group, mepos)
    TYPE(message_entry_type),INTENT(INOUT) :: entry
    INTEGER, INTENT(IN)                    :: src, group, mepos

    INTEGER, DIMENSION(key_length) :: key_arr
    INTEGER                        :: datatype, datasize
    INTEGER, DIMENSION(default_string_length) :: value_str_arr

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
          STOP "swarm_message_entry_mpi_bcast: no value ASSOCIATED"
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
    INTEGER, DIMENSION(default_string_length):: value_str_arr

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
       value_str_arr = TRANSFER(entry%value_str, value_str_arr)
       CALL mp_send(value_str_arr, dest, tag, group)
    ELSE
       CALL abort()
       STOP "swarm_message_entry_mpi_send: no value ASSOCIATED"
    END IF
  END SUBROUTINE swarm_message_entry_mpi_send


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_recv(entry, group, src, tag)
    TYPE(message_entry_type),INTENT(INOUT) :: entry
    INTEGER, INTENT(IN)                    :: group
    INTEGER, INTENT(INOUT)                 :: src, tag

    INTEGER, DIMENSION(key_length) :: key_arr
    INTEGER, DIMENSION(default_string_length) :: value_str_arr
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
          curr_entry => curr_entry%next
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
  SUBROUTINE swarm_message_file_write(msg, unit)
    TYPE(swarm_message_type),INTENT(IN) :: msg
    INTEGER, INTENT(IN)                 :: unit

    TYPE(message_entry_type), POINTER :: curr_entry

    WRITE(unit,"(A)") "BEGIN SWARM_MESSAGE"
    WRITE(unit,"(AI10)") "msg_length: ", swarm_message_length(msg)

    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      CALL swarm_message_entry_file_write(curr_entry, unit)
      curr_entry => curr_entry%next
    END DO

    WRITE(unit,"(A)") "END SWARM_MESSAGE"
    WRITE(unit,"()")
  END SUBROUTINE swarm_message_file_write


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_entry_file_write(entry, unit)
    TYPE(message_entry_type),INTENT(IN) :: entry
    INTEGER, INTENT(IN)                 :: unit

    INTEGER                             :: i

    WRITE(unit,"(AA)") "key: ",  entry%key
    IF(ASSOCIATED(entry%value_i4)) THEN
       WRITE(unit,"(A)") "datatype: i4"
       WRITE(unit,"(AI10)") "value: ", entry%value_i4

    ELSE IF(ASSOCIATED(entry%value_i8)) THEN
       WRITE(unit,"(A)") "datatype: i8"
       WRITE(unit,"(AI20)") "value: ", entry%value_i8

    ELSE IF(ASSOCIATED(entry%value_r4)) THEN
       WRITE(unit,"(A)") "datatype: r4"
       WRITE(unit,"(AE30.20)") "value: ", entry%value_r4

    ELSE IF(ASSOCIATED(entry%value_r8)) THEN
       WRITE(unit,"(A)") "datatype: r8"
       WRITE(unit,"(AE30.20)") "value: ", entry%value_r8

    ELSE IF(ASSOCIATED(entry%value_str)) THEN
       WRITE(unit,"(A)") "datatype: str"
       WRITE(unit,"(AA)") "value: ", entry%value_str

    ELSE IF(ASSOCIATED(entry%value_i4_1d)) THEN
       WRITE(unit,"(A)") "datatype: i4_1d"
       WRITE(unit,"(AI10)") "size: ", SIZE(entry%value_i4_1d)
       WRITE(unit,*) entry%value_i4_1d

    ELSE IF(ASSOCIATED(entry%value_i8_1d)) THEN
       WRITE(unit,"(A)") "datatype: i8_1d"
       WRITE(unit,"(AI20)") "size: ", SIZE(entry%value_i8_1d)
       WRITE(unit,*) entry%value_i8_1d

    ELSE IF(ASSOCIATED(entry%value_r4_1d)) THEN
       WRITE(unit,"(A)") "datatype: r4_1d"
       WRITE(unit,"(AI8)") "size: ", SIZE(entry%value_r4_1d)
       DO i=1, SIZE(entry%value_r4_1d)
         WRITE(unit,"(XE30.20)", advance="no") entry%value_r4_1d(i)
       END DO
       WRITE(unit,"(X)")

    ELSE IF(ASSOCIATED(entry%value_r8_1d)) THEN
       WRITE(unit,"(A)") "datatype: r8_1d"
       WRITE(unit,"(AI8)") "size: ", SIZE(entry%value_r8_1d)
       DO i=1, SIZE(entry%value_r8_1d)
         WRITE(unit,"(XE30.20)", advance="no") entry%value_r8_1d(i)
       END DO
       WRITE(unit,"(X)")

    ELSE
       CALL abort()
       STOP "swarm_message_entry_file_write: no value ASSOCIATED"
    END IF
  END SUBROUTINE swarm_message_entry_file_write


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_file_read(msg, parser, at_end, error)
    TYPE(swarm_message_type),INTENT(OUT)       :: msg
    TYPE(cp_parser_type), POINTER              :: parser
    LOGICAL, INTENT(INOUT)                     :: at_end
    TYPE(cp_error_type)                        :: error

    TYPE(message_entry_type), POINTER          :: new_entry
    INTEGER                                    :: length, i
    CHARACTER(LEN=default_string_length)       :: key

    !IF(ASSOCIATED(msg%root)) STOP "swarm_message_file_read: message not empty"

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN
    IF(TRIM(parser%input_line)/="BEGIN SWARM_MESSAGE") STOP "swarm_message_file_read failed (1)"

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN
    READ (parser%input_line,*) key, length
    IF(TRIM(key)/="msg_length:") STOP "swarm_message_file_read failed (2)"

    DO i=1, length
       ALLOCATE(new_entry)
       CALL swarm_message_entry_file_read(new_entry, parser, at_end, error)
       new_entry%next => msg%root
       msg%root => new_entry
    END DO

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN
    IF(TRIM(parser%input_line)/="END SWARM_MESSAGE") STOP "swarm_message_file_read failed (3)"

  END SUBROUTINE swarm_message_file_read

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE swarm_message_entry_file_read(entry, parser, at_end, error)
    TYPE(message_entry_type),INTENT(INOUT)     :: entry
    TYPE(cp_parser_type), POINTER              :: parser
    LOGICAL, INTENT(INOUT)                     :: at_end
    TYPE(cp_error_type)                        :: error
    CHARACTER(LEN=default_string_length)       :: label, datatype

    INTEGER                                    :: i, arr_size
    LOGICAL                                    :: is_scalar

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN
    READ (parser%input_line,*) label, entry%key
    IF(TRIM(label)/="key:") STOP "swarm_message_entry_file_read failed (1)"

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN
    READ (parser%input_line,*) label, datatype
    IF(TRIM(label)/="datatype:") STOP "swarm_message_entry_file_read failed (2)"

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN

    is_scalar=.TRUE.
    IF(TRIM(datatype)=="i4") THEN
       ALLOCATE(entry%value_i4)
       READ (parser%input_line,*) label, entry%value_i4
    ELSE IF(TRIM(datatype)=="i8") THEN
       ALLOCATE(entry%value_i8)
       READ (parser%input_line,*) label, entry%value_i8
    ELSE IF(TRIM(datatype)=="r4") THEN
       ALLOCATE(entry%value_r4)
       READ (parser%input_line,*) label, entry%value_r4
    ELSE IF(TRIM(datatype)=="r8") THEN
       ALLOCATE(entry%value_r8)
       READ (parser%input_line,*) label, entry%value_r8
    ELSE IF(TRIM(datatype)=="str") THEN
       ALLOCATE(entry%value_str)
       READ (parser%input_line,*) label, entry%value_str
    ELSE
       is_scalar = .FALSE.
    END IF

    IF(is_scalar) THEN
       IF(TRIM(label)/="value:") STOP "swarm_message_entry_file_read failed (3)"
       RETURN
    END IF

    ! musst be an array-datatype
    READ (parser%input_line,*) label, arr_size
    IF(TRIM(label)/="size:") STOP "swarm_message_entry_file_read failed (4)"

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN

    IF(TRIM(datatype)=="i4_1d") THEN
       ALLOCATE(entry%value_i4_1d(arr_size))
       READ (parser%input_line,*) entry%value_i4_1d

    ELSE IF(TRIM(datatype)=="i8_1d") THEN
       ALLOCATE(entry%value_i8_1d(arr_size))
       READ (parser%input_line,*) entry%value_i8_1d

    ELSE IF(TRIM(datatype)=="r4_1d") THEN
       ALLOCATE(entry%value_r4_1d(arr_size))
       READ (parser%input_line,*) entry%value_r4_1d

    ELSE IF(TRIM(datatype)=="r8_1d") THEN
       ALLOCATE(entry%value_r8_1d(arr_size))
       READ (parser%input_line,*) entry%value_r8_1d
    ELSE
       STOP "swarm_message_entry_file_read: unkown datatype"
    END IF
  END SUBROUTINE swarm_message_entry_file_read

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


