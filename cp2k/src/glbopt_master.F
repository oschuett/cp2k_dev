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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE glbopt_types,                    ONLY: GLBOPT_CMD_SHUTDOWN,&
                                             GLBOPT_CMD_MDGOPT,&
                                             GLBOPT_CMD_SETIFRAME,&
                                             glbopt_command_type,&
                                             glbopt_report_type
  USE glbopt_history,                  ONLY: glbopt_history_type,&
                                             glbopt_history_init,&
                                             glbopt_history_finalize,&
                                             glbopt_history_add,&
                                             glbopt_history_lookup
  USE hash_functions,                  ONLY: b3hs_hash_key_jenkins
  USE kinds,                           ONLY: dp,&
                                             int_4,default_path_length,default_string_length
  USE message_passing,                 ONLY: mp_any_source,mp_any_tag,&
                                             mp_recv,&
                                             mp_send
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_type
  USE physcon,                         ONLY: kelvin, femtoseconds

  USE glbopt_minhop,                   ONLY: glbopt_minhop_type,&
                                             glbopt_minhop_init,&
                                             glbopt_minhop_finalize,&
                                             glbopt_minhop_steer
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
 USE timings_types,                   ONLY: timer_env_type
 USE timings,                         ONLY: add_timer_env, rm_timer_env, timer_env_create
 USE input_constants,                 ONLY: dump_atomic,&
                                             dump_dcd,&
                                             dump_pdb,&
                                             dump_xmol
 USE cp_parser_types,                  ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
 USE cp_parser_methods,                ONLY: parser_get_next_line, parser_get_object

 USE glbopt_comlog,                    ONLY: glbopt_comlog_write, glbopt_comlog_read
 USE cp_parser_types,                  ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
 USE cp_parser_methods,                ONLY: parser_get_next_line, parser_get_object
 USE global_types,                     ONLY: global_environment_type
 USE cp_external_control,              ONLY: external_control


#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_master'


 PUBLIC :: glbopt_master_type
 PUBLIC :: glbopt_master_init, glbopt_master_finalize
 PUBLIC :: glbopt_master_steer


 TYPE glbopt_master_type
   PRIVATE
   TYPE(glbopt_minhop_type), ALLOCATABLE :: minhop
   REAL(KIND=dp)                                       :: Elowest = Huge(1.0_dp)
   REAL(KIND=dp)                                       :: Emin = Tiny(1.0_dp)
   INTEGER                                             :: iw = 0
   INTEGER                                             :: i_iteration = 0
   INTEGER                                             :: max_iter = 0
   LOGICAL                                             :: should_stop = .FALSE.
   INTEGER                                             :: n_walkers = -1
   INTEGER                                             :: comlog_unit
   TYPE(cp_error_type)                                 :: error
   TYPE(section_vals_type), POINTER                    :: glbopt_section => Null()
   TYPE(cp_para_env_type), POINTER                     :: para_env => Null()
   TYPE(cmd_linked_list_type), POINTER                 :: command_queue => Null()
   TYPE(global_environment_type), POINTER              :: globenv => Null()
 END TYPE glbopt_master_type

 TYPE cmd_linked_list_type
   TYPE(glbopt_command_type)           :: cmd
   TYPE(cmd_linked_list_type), POINTER :: next
 END TYPE cmd_linked_list_type

 CONTAINS



! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_init(master, para_env, globenv, root_section, n_walkers, error)
    TYPE(glbopt_master_type)                 :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: n_walkers
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: output_unit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(timer_env_type), POINTER            :: timer_env


    NULLIFY(timer_env)
    CALL timer_env_create(timer_env)
    CALL add_timer_env(timer_env)
    !CALL add_mp_perf_env(f_env%mp_perf_env) TODO: this as well?

    master%glbopt_section => section_vals_get_subs_vals(root_section, "GLOBAL_OPT", error=error)

    logger => cp_error_get_logger(error)
    master%n_walkers = n_walkers
    master%error = error
    master%para_env => para_env
    master%globenv => globenv
    master%iw = cp_print_key_unit_nr(logger, master%glbopt_section,&
          "PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    master%comlog_unit = cp_print_key_unit_nr(logger,master%glbopt_section,"PRINT%COMMUNICATION_LOG",&
            !extension=my_ext, file_position=my_pos, file_action=my_act,&
            !file_form=my_form, middle_name=TRIM(my_middle), is_new_file=new_file,&
            extension=".xyz", error=error)

    CALL section_vals_val_get(master%glbopt_section,"E_MIN", r_val=master%Emin, error=error)
    CALL section_vals_val_get(master%glbopt_section,"MAX_ITER", i_val=master%max_iter, error=error)

    ALLOCATE(master%minhop)
    CALL glbopt_minhop_init(master%minhop, master%glbopt_section, n_walkers, master%iw, error)

    CALL replay_comlog(master)
 END SUBROUTINE glbopt_master_init




! *****************************************************************************
! *****************************************************************************
 SUBROUTINE replay_comlog(master)
    TYPE(glbopt_master_type) :: master
    CHARACTER(LEN=default_path_length) :: filename
    TYPE(cp_parser_type), POINTER :: parser
    LOGICAL :: explicit
    TYPE(glbopt_report_type) :: report_log
    TYPE(glbopt_command_type) :: cmd_log, cmd_now
    TYPE(cmd_linked_list_type), POINTER :: cmd_ll
    LOGICAL :: at_end
    TYPE(glbopt_command_type), DIMENSION(master%n_walkers) :: last_command
    INTEGER, DIMENSION(master%n_walkers) :: last_iframe
    INTEGER                             :: i

    last_iframe(:) = 0
    last_command(:)%cmd_id = -1

    NULLIFY(parser)

    ! Initialize parser for trajectory
    CALL section_vals_val_get(master%glbopt_section,"REPLAY_COMMUNICATION_LOG",&
       c_val=filename, explicit=explicit, error=master%error)

    IF(.NOT. explicit) RETURN
    IF(master%iw>0) WRITE(master%iw,'(A,A)') &
       " GLOBAL_OPT| Starting replay of communication-log: ", TRIM(filename)

    CALL parser_create(parser, filename, para_env=master%para_env, error=master%error)

    at_end = .FALSE.
    DO
      CALL glbopt_comlog_read(parser, report_log, at_end, master%error)
      IF(at_end) EXIT
      CALL glbopt_comlog_read(parser, cmd_log, at_end, master%error)
      IF(at_end) EXIT
      CALL glbopt_master_steer(master, report_log, cmd_now)
      !TODO: maybe we should just exit the loop instead of stopping?
      IF(.NOT. commands_equal(cmd_now, cmd_log)) STOP "wrong behaviour"

      last_command(cmd_log%walker_id) = cmd_log
      last_iframe(report_log%walker_id) = report_log%iframe
    END DO

    IF(master%iw>0) WRITE(master%iw,'(A,A)') &
       " GLOBAL_OPT| Reached end of communication log. Queueing last commands."

    DO i=1, master%n_walkers
      ALLOCATE(cmd_ll)
      cmd_ll%cmd%walker_id = i
      cmd_ll%cmd%cmd_id = GLBOPT_CMD_SETIFRAME
      cmd_ll%cmd%iframe = last_iframe(i)
      cmd_ll%next => master%command_queue
      master%command_queue => cmd_ll

      IF(last_command(i)%cmd_id > 0) THEN
         ALLOCATE(cmd_ll)
         cmd_ll%cmd = last_command(i)
         cmd_ll%next => master%command_queue
         master%command_queue => cmd_ll
      END IF
    END DO

    CALL parser_release(parser, error=master%error)
 END SUBROUTINE replay_comlog


! *****************************************************************************
! *****************************************************************************
 FUNCTION commands_equal(cmd1, cmd2) RESULT(res)
    TYPE(glbopt_command_type)                :: cmd1, cmd2
    LOGICAL :: res
    res =           cmd1%walker_id == cmd2%walker_id
    res = res .AND. cmd1%cmd_id == cmd2%cmd_id
    res = res .AND. ABS(cmd1%temperature - cmd2%temperature) < 1e-10
    res = res .AND. ALL(ABS(cmd1%positions -cmd2%positions) < 1e-10)
    !WRITE (*,*) cmd1%walker_id, cmd2%walker_id
    !WRITE (*,*) cmd1%cmd_id, cmd2%cmd_id
    !WRITE (*,*) ABS(cmd1%temperature - cmd2%temperature)
 END FUNCTION commands_equal


! *****************************************************************************
! *****************************************************************************
  FUNCTION command_queue_pop(master, walker_id) RESULT(cmd)
    TYPE(glbopt_master_type)   :: master
    INTEGER                    :: walker_id
    TYPE(glbopt_command_type)  :: cmd
    TYPE(cmd_linked_list_type), POINTER :: cmd_ll, cmd_ll_prev

    cmd%cmd_id = -1 ! default result = nothing found

    NULLIFY(cmd_ll_prev)
    cmd_ll => master%command_queue
    DO WHILE(ASSOCIATED(cmd_ll))
       IF(cmd_ll%cmd%walker_id /= walker_id) THEN
          cmd_ll_prev => cmd_ll
          cmd_ll => cmd_ll%next
          CYCLE
       END IF

       ! we found a queued command for this walker
       cmd = cmd_ll%cmd
       IF(ASSOCIATED(cmd_ll_prev)) THEN
          cmd_ll_prev%next => cmd_ll%next
       ELSE
           master%command_queue => cmd_ll%next
       END IF
       DEALLOCATE(cmd_ll)
       EXIT
    END DO

  END FUNCTION command_queue_pop


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_steer(master, report, cmd)
    TYPE(glbopt_master_type)   :: master
    TYPE(glbopt_report_type)   :: report
    TYPE(glbopt_command_type)  :: cmd

    ! First check if there are queued commands for this walker
    cmd = command_queue_pop(master, report%walker_id)
    IF(cmd%cmd_id > 0) THEN
       IF(master%iw>0) WRITE (master%iw,'(A,A,A,I9,X,A)') ' GLOBAL_OPT| ',&
         REPEAT("*", 9), " Sending out queued command to walker: ",&
         cmd%walker_id,REPEAT("*", 9)
       RETURN
    END IF

    ! There are no queued commands. Do the normal processing.
    cmd%walker_id = report%walker_id
    master%Elowest = MIN(master%Elowest, report%Epot)
    master%i_iteration = master%i_iteration + 1
    IF(master%iw>0) THEN
       WRITE (master%iw,'(A,A,X,I8,A,A)') ' GLOBAL_OPT| ',REPEAT("*", 13),&
          master%i_iteration,'. Global Optimisation Iteration ',REPEAT("*", 13)
       WRITE (master%iw,'(A,42X,I8)')&
              " GLOBAL_OPT| Reporting walker ", report%walker_id
       WRITE (master%iw,'(A,16X,E15.8)')&
              " GLOBAL_OPT| Reported potential Energy [Hartree] ", report%Epot
       WRITE (master%iw,'(A,9X,E15.8)')&
        " GLOBAL_OPT| Lowest reported potential Energy [Hartree] ",master%Elowest
    END IF

    IF(master%Elowest<master%Emin .AND. .NOT. master%should_stop) THEN
       IF(master%iw>0)  WRITE (master%iw, *)  " GLOBAL_OPT| Reached Emin > Epot. Quitting."
       master%should_stop = .TRUE.
    ENDIF

    IF(master%i_iteration>=master%max_iter .AND. .NOT. master%should_stop) THEN
       IF(master%iw>0)  WRITE (master%iw, *)  " GLOBAL_OPT| Reached MAX_ITER. Quitting."
       master%should_stop = .TRUE.
    ENDIF

    IF(.NOT. master%should_stop) &
      CALL external_control(master%should_stop,"GLBOPT",master%globenv,error=master%error)

    CALL glbopt_comlog_write(master%comlog_unit, master%i_iteration, report)
    CALL glbopt_minhop_steer(master%minhop, report, cmd)
    CALL glbopt_comlog_write(master%comlog_unit, master%i_iteration, cmd)

    IF(master%should_stop) THEN
       cmd%cmd_id = GLBOPT_CMD_SHUTDOWN  !overwrite command
       IF(master%iw>0) WRITE(master%iw,'(A,25X,I8)') &
         " GLOBAL_OPT| Sending shutdown command to walker",cmd%walker_id
    ENDIF

 END SUBROUTINE


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_finalize(master)
    TYPE(glbopt_master_type)                 :: master
    TYPE(cp_logger_type), POINTER            :: logger

    IF(master%iw > 0) WRITE(master%iw,"(A)") " GLOBAL_OPT| Shutting down the master."

    CALL glbopt_minhop_finalize(master%minhop)
    DEALLOCATE(master%minhop)

    logger => cp_error_get_logger(master%error)
    CALL cp_print_key_finished_output(master%iw, logger,&
            master%glbopt_section, "PRINT%MASTER_RUN_INFO", error=master%error)
    CALL cp_print_key_finished_output(master%comlog_unit, logger,&
            master%glbopt_section, "PRINT%COMMUNICATION_LOG", error=master%error)

    CALL rm_timer_env() !pops the top-most timer
 END SUBROUTINE glbopt_master_finalize

END MODULE glbopt_master

