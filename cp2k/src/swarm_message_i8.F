 SUBROUTINE swarm_message_add_i8(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
    INTEGER(KIND=int_8)                                    :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) THEN
      WRITE (*,*) "swarm_message_add_i8: key already exists: ",key
      CALL ABORT() !TODO: call mp_abort
   END IF

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_i8)
   new_entry%value_i8 = value

   !WRITE (*,*) "swarm_message_add_i8: key=",key, " value=",new_entry%value_i8

   !TODO: check if name already exists
   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_i8

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE swarm_message_get_i8(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
    INTEGER(KIND=int_8)                                :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_i8: key=",key

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_i8)) &
            STOP "swarm_message_get_i8: value not associated"
            !TODO: check other datatypes, call mp_abort
         value = curr_entry%value_i8
         !WRITE (*,*) "swarm_message_get_i8: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL mp_abort("swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_i8


