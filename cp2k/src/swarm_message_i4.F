 SUBROUTINE swarm_message_add_i4(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
    INTEGER(KIND=int_4), INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_i4)
   new_entry%value_i4 = value

   WRITE (*,*) "swarm_message_add_i4: key=",key, " value=",new_entry%value_i4

   !TODO: check if name already exists
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
   WRITE (*,*) "swarm_message_get_i4: key=",key

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_i4)) &
            STOP "swarm_message_get_i4: value not associated"
            !TODO: check other datatypes
         value = curr_entry%value_i4
         WRITE (*,*) "swarm_message_get_i4: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   WRITE (*,*) "swarm_message_get: key not found: ",key
   !CALL ABORT()
   STOP 1
   !TODO: call mp_abort
 END SUBROUTINE swarm_message_get_i4


