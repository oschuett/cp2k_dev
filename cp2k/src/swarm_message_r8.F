 SUBROUTINE swarm_message_add_r8(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
    REAL(KIND=real_8), INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) THEN
      WRITE (*,*) "swarm_message_add_r8: key already exists: ",key
      CALL ABORT() !TODO: call mp_abort
   END IF

   ALLOCATE(new_entry)
   new_entry%key = key

   new_entry%value_r8 = value

   !WRITE (*,*) "swarm_message_add_r8: key=",key, " value=",new_entry%value_r8

   !TODO: check if name already exists
   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_r8

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE swarm_message_get_r8(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
    REAL(KIND=real_8) , INTENT(OUT)         :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_r8: key=",key

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ALLOCATED(curr_entry%value_r8)) &
            STOP "swarm_message_get_r8: value not associated"
            !TODO: check other datatypes, call mp_abort
         value = curr_entry%value_r8
         !WRITE (*,*) "swarm_message_get_r8: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL mp_abort("swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_r8


