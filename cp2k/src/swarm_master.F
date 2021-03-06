!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE swarm_master
  USE cp_external_control,             ONLY: external_control
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE glbopt_master,                  ONLY:  glbopt_master_finalize,&
                                             glbopt_master_init,&
                                             glbopt_master_steer,&
                                             glbopt_master_type
  USE paramopt_master,                ONLY:  paramopt_master_finalize,&
                                             paramopt_master_init,&
                                             paramopt_master_steer,&
                                             paramopt_master_type
  USE swarm_message,                   ONLY: swarm_message_type,&
                                             swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_length,&
                                             swarm_message_file_write,&
                                             swarm_message_file_read,&
                                             swarm_message_equal,&
                                             swarm_message_free,&
                                             swarm_message_get_str
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE timings,                         ONLY: add_timer_env,&
                                             rm_timer_env,&
                                             timer_env_create
  USE timings_types,                   ONLY: timer_env_type
  USE input_constants,                 ONLY: swarm_do_glbopt,&
                                             swarm_do_paramopt
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_master'


 PUBLIC :: swarm_master_type
 PUBLIC :: swarm_master_init, swarm_master_finalize
 PUBLIC :: swarm_master_steer

 TYPE swarm_message_p_type
   TYPE(swarm_message_type), ALLOCATABLE               :: p
 END TYPE swarm_message_p_type


 TYPE swarm_master_type
   PRIVATE
   INTEGER                                             :: behavior = -1
   TYPE(glbopt_master_type), ALLOCATABLE               :: glbopt
   TYPE(paramopt_master_type), ALLOCATABLE             :: paramopt
   INTEGER                                             :: iw = 0
   INTEGER                                             :: i_iteration = 0
   INTEGER                                             :: max_iter = 0
   LOGICAL                                             :: should_stop = .FALSE.
   INTEGER                                             :: n_workers = -1
   INTEGER                                             :: comlog_unit
   TYPE(cp_error_type)                                 :: error
   TYPE(section_vals_type), POINTER                    :: swarm_section => Null()
   TYPE(cp_para_env_type), POINTER                     :: para_env => Null()
   TYPE(swarm_message_p_type), DIMENSION(:), ALLOCATABLE :: queued_commands
   TYPE(global_environment_type), POINTER              :: globenv => Null()
 END TYPE swarm_master_type


 CONTAINS



! *****************************************************************************
! *****************************************************************************
! TODO: input_path and root_section not really need
 SUBROUTINE swarm_master_init(master, para_env, globenv, root_section, input_path, n_workers, error)
    TYPE(swarm_master_type)                 :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: n_workers
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(timer_env_type), POINTER            :: timer_env

    NULLIFY(timer_env)
    CALL timer_env_create(timer_env)
    CALL add_timer_env(timer_env)
    !CALL add_mp_perf_env(f_env%mp_perf_env) TODO: this as well?

    master%swarm_section => section_vals_get_subs_vals(root_section, "SWARM", error=error)

    logger => cp_error_get_logger(error)
    master%n_workers = n_workers
    master%error = error
    master%para_env => para_env
    master%globenv => globenv
    ALLOCATE(master%queued_commands(master%n_workers))
    master%iw = cp_print_key_unit_nr(logger, master%swarm_section,&
          "PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    CALL section_vals_val_get(master%swarm_section,"BEHAVIOR", i_val=master%behavior, error=error)

    ! uses logger%iter_info%project_name to construct filename
    master%comlog_unit = cp_print_key_unit_nr(logger,master%swarm_section,"PRINT%COMMUNICATION_LOG",&
            !middle_name="comlog", extension=".xyz", &
            extension=".comlog", &
            file_action="WRITE", file_position="REWIND", error=error)

    CALL section_vals_val_get(master%swarm_section,"MAX_ITER", i_val=master%max_iter, error=error)

    SELECT CASE(master%behavior)
      CASE(swarm_do_glbopt)
         ALLOCATE(master%glbopt)
         CALL glbopt_master_init(master%glbopt, para_env, master%swarm_section, input_path, n_workers, master%iw, error)
      CASE(swarm_do_paramopt)
         ALLOCATE(master%paramopt)
         CALL paramopt_master_init(master%paramopt, para_env, master%swarm_section, input_path, n_workers, master%iw, error)
      CASE DEFAULT
         STOP "swarm_master_init: got unkown behavior"
    END SELECT

    CALL replay_comlog(master)
 END SUBROUTINE swarm_master_init




! *****************************************************************************
! *****************************************************************************
 SUBROUTINE replay_comlog(master)
    TYPE(swarm_master_type)                 :: master

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: i
    TYPE(swarm_message_type), &
      DIMENSION(master%n_workers)            :: last_commands
    TYPE(swarm_message_type)                 :: cmd_log
    TYPE(swarm_message_type), POINTER        :: cmd_now
    TYPE(cp_parser_type), POINTER            :: parser
    !TYPE(cmd_linked_list_type), POINTER      :: cmd_ll
    LOGICAL                                  :: at_end, explicit
    TYPE(swarm_message_type)                 :: report_log

    INTEGER                                  :: worker_id
    NULLIFY(parser)

    ! Initialize parser for trajectory
    CALL section_vals_val_get(master%swarm_section,"REPLAY_COMMUNICATION_LOG",&
       c_val=filename, explicit=explicit, error=master%error)

    IF(.NOT. explicit) RETURN
    IF(master%iw>0) WRITE(master%iw,'(A,A)') &
       " SWARM| Starting replay of communication-log: ", TRIM(filename)

    CALL parser_create(parser, filename, para_env=master%para_env, error=master%error)

    at_end = .FALSE.
    DO
      CALL swarm_message_file_read(report_log, parser, at_end, master%error)
      IF(at_end) EXIT

      CALL swarm_message_file_read(cmd_log, parser, at_end, master%error)
      IF(at_end) EXIT

      ALLOCATE(cmd_now)
      CALL swarm_master_steer(master, report_log, cmd_now)

      !TODO: maybe we should just exit the loop instead of stopping?
      IF(TRIM(swarm_message_get_str(cmd_log, "command")) /= "shutdown") THEN
        IF(.NOT. commands_equal(cmd_now, cmd_log, master%iw)) STOP "wrong behaviour"
      END IF

      CALL swarm_message_free(cmd_log)
      CALL swarm_message_free(report_log)
      CALL swarm_message_get(cmd_now, "worker_id", worker_id)
      CALL swarm_message_free(last_commands(worker_id))
      last_commands(worker_id) = cmd_now
      DEALLOCATE(cmd_now)
    END DO

    CALL swarm_message_free(report_log) !don't worry about double-frees
    CALL swarm_message_free(cmd_log)
    !TODO: free queued commands

    IF(master%iw>0) WRITE(master%iw,'(A,A)') &
       " SWARM| Reached end of communication log. Queueing last commands."

    DO i=1, master%n_workers
       master%queued_commands(i)%p = last_commands(i)
    END DO

    CALL parser_release(parser, error=master%error)
 END SUBROUTINE replay_comlog


! *****************************************************************************
! *****************************************************************************
 FUNCTION commands_equal(cmd1, cmd2, iw) RESULT(res)
    TYPE(swarm_message_type)                 :: cmd1, cmd2
    LOGICAL                                  :: res
    INTEGER                                  :: iw

    res = swarm_message_equal(cmd1, cmd2)
    IF(.NOT. res .AND. iw>0) THEN
       WRITE(iw, *) "Command 1:"
       CALL swarm_message_file_write(cmd1, iw)
       WRITE(iw, *) "Command 2:"
       CALL swarm_message_file_write(cmd2, iw)
    END IF
 END FUNCTION commands_equal



! *****************************************************************************
! *****************************************************************************
 SUBROUTINE swarm_master_steer(master, report, cmd)
    TYPE(swarm_master_type), INTENT(INOUT)        :: master
    TYPE(swarm_message_type), INTENT(IN)          :: report
    TYPE(swarm_message_type), INTENT(OUT)         :: cmd

    INTEGER                                       :: worker_id
    LOGICAL                                       :: should_stop
    should_stop = .FALSE.

    ! First check if there are queued commands for this worker
    CALL swarm_message_get(report, "worker_id", worker_id)

    IF(ALLOCATED(master%queued_commands(worker_id)%p)) THEN
       cmd = master%queued_commands(worker_id)%p
       DEALLOCATE(master%queued_commands(worker_id)%p)
       IF(master%iw>0) WRITE (master%iw,'(A,A,A,I9,X,A)') ' SWARM| ',&
          REPEAT("*", 9), " Sending out queued command to worker: ",&
          worker_id,REPEAT("*", 9)
       RETURN
    END IF

    ! There are no queued commands. Do the normal processing.
    master%i_iteration = master%i_iteration + 1

    IF(master%iw>0) WRITE (master%iw,'(A,A,X,I8,A,A)') ' SWARM| ',REPEAT("*", 15),&
         master%i_iteration,'. Master / Worker Communication  ',REPEAT("*", 15)

    IF(master%i_iteration>=master%max_iter .AND. .NOT. master%should_stop) THEN
       IF(master%iw>0)  WRITE (master%iw, *)  "SWARM| Reached MAX_ITER. Quitting."
       master%should_stop = .TRUE.
    ENDIF

    IF(.NOT. master%should_stop) THEN
      CALL external_control(master%should_stop,"SWARM",master%globenv,error=master%error)
      IF(master%should_stop .AND. master%iw>0) &
         WRITE (master%iw, *)  "SWARM| Received stop from external_control. Quitting."
    END IF

    !IF(unit > 0) &

    CALL swarm_message_file_write(report, master%comlog_unit) ! master%i_iteration


    CALL swarm_message_add(cmd, "worker_id", worker_id)

    IF(master%should_stop) THEN
       CALL swarm_message_add(cmd, "command", "shutdown")
       IF(master%iw>0) WRITE(master%iw,'(A,25X,I8)') &
         " SWARM| Sending shutdown command to worker",worker_id
    ELSE
       SELECT CASE(master%behavior)
         CASE(swarm_do_glbopt)
            CALL glbopt_master_steer(master%glbopt, report, cmd, should_stop)
         CASE(swarm_do_paramopt)
            CALL paramopt_master_steer(master%paramopt, report, cmd, should_stop)
         CASE DEFAULT
            STOP "swarm_master_steer: got unkown behavior"
       END SELECT

       IF(should_stop) THEN
          !TODO: does overwritting entries in swarm_message actually work?
          CALL swarm_message_add(cmd, "command", "shutdown") !overwrite command
          IF(master%iw>0) WRITE(master%iw,'(A,25X,I8)') &
              " SWARM| Sending shutdown command to worker",worker_id
           master%should_stop = .TRUE.
       ENDIF
    ENDIF

    CALL swarm_message_file_write(cmd, master%comlog_unit) ! master%i_iteration

 END SUBROUTINE swarm_master_steer


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE swarm_master_finalize(master)
    TYPE(swarm_master_type)                 :: master

    TYPE(cp_logger_type), POINTER            :: logger

    IF(master%iw > 0) WRITE(master%iw,"(A)") " SWARM| Shutting down the master."

    SELECT CASE(master%behavior)
      CASE(swarm_do_glbopt)
         CALL glbopt_master_finalize(master%glbopt)
         DEALLOCATE(master%glbopt)
      CASE(swarm_do_paramopt)
         CALL paramopt_master_finalize(master%paramopt)
         DEALLOCATE(master%paramopt)
      CASE DEFAULT
         STOP "swarm_master_finalize: got unkown behavior"
    END SELECT

    logger => cp_error_get_logger(master%error)
    CALL cp_print_key_finished_output(master%iw, logger,&
            master%swarm_section, "PRINT%MASTER_RUN_INFO", error=master%error)
    CALL cp_print_key_finished_output(master%comlog_unit, logger,&
            master%swarm_section, "PRINT%COMMUNICATION_LOG", error=master%error)

    CALL rm_timer_env() !pops the top-most timer
 END SUBROUTINE swarm_master_finalize

END MODULE swarm_master

