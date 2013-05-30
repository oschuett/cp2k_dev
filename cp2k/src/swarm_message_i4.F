 SUBROUTINE swarm_message_add_i4(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
    INTEGER(KIND=int_4), INTENT(IN), TARGET                :: value

   TYPE(message_entry_type), POINTER :: new_entry

   ALLOCATE(new_entry)
   new_entry%key = key
   new_entry%value_i4 => value

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_i4

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE swarm_message_get_i4(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
    INTEGER(KIND=int_4), INTENT(OUT)                   :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         value = curr_entry%value_i4
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   STOP "swarm_message_get: key not found"
 END SUBROUTINE swarm_message_get_i4


