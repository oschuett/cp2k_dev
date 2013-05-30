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
  USE swarm_comlog,                   ONLY: swarm_comlog_read,&
                                             swarm_comlog_write
  USE glbopt_master,                  ONLY:  glbopt_master_finalize,&
                                             glbopt_master_init,&
                                             glbopt_master_steer,&
                                             glbopt_master_type
  USE swarm_types,                    ONLY:  SWARM_CMD_SHUTDOWN,&
                                             swarm_command_type,&
                                             swarm_report_type
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
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_master'


 PUBLIC :: swarm_master_type
 PUBLIC :: swarm_master_init, swarm_master_finalize
 PUBLIC :: swarm_master_steer


 TYPE swarm_master_type
   PRIVATE
   TYPE(glbopt_master_type), ALLOCATABLE               :: glbopt
   REAL(KIND=dp)                                       :: Elowest = HUGE(1.0_dp)
   REAL(KIND=dp)                                       :: Emin = TINY(1.0_dp)
   INTEGER                                             :: iw = 0
   INTEGER                                             :: i_iteration = 0
   INTEGER                                             :: max_iter = 0
   LOGICAL                                             :: should_stop = .FALSE.
   INTEGER                                             :: n_workers = -1
   INTEGER                                             :: comlog_unit
   TYPE(cp_error_type)                                 :: error
   TYPE(section_vals_type), POINTER                    :: swarm_section => Null()
   TYPE(cp_para_env_type), POINTER                     :: para_env => Null()
   TYPE(cmd_linked_list_type), POINTER                 :: command_queue => Null()
   TYPE(global_environment_type), POINTER              :: globenv => Null()
 END TYPE swarm_master_type

 TYPE cmd_linked_list_type
   TYPE(swarm_command_type)           :: cmd
   TYPE(cmd_linked_list_type), POINTER :: next
 END TYPE cmd_linked_list_type

 CONTAINS



! *****************************************************************************
! *****************************************************************************
 SUBROUTINE swarm_master_init(master, para_env, globenv, root_section, n_workers, error)
    TYPE(swarm_master_type)                 :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: root_section
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
    master%iw = cp_print_key_unit_nr(logger, master%swarm_section,&
          "PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    ! uses logger%iter_info%project_name to construct filename
    master%comlog_unit = cp_print_key_unit_nr(logger,master%swarm_section,"PRINT%COMMUNICATION_LOG",&
            middle_name="comlog", extension=".xyz", &
            file_action="WRITE", file_position="REWIND", error=error)

    CALL section_vals_val_get(master%swarm_section,"E_MIN", r_val=master%Emin, error=error)
    CALL section_vals_val_get(master%swarm_section,"MAX_ITER", i_val=master%max_iter, error=error)

    ALLOCATE(master%glbopt)
    CALL glbopt_master_init(master%glbopt, master%swarm_section, n_workers, master%iw, error)

    CALL replay_comlog(master)
 END SUBROUTINE swarm_master_init




! *****************************************************************************
! *****************************************************************************
 SUBROUTINE replay_comlog(master)
    TYPE(swarm_master_type)                 :: master

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: i
    TYPE(swarm_command_type), &
      DIMENSION(master%n_workers)            :: last_command
    TYPE(swarm_command_type)                :: cmd_log, cmd_now
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cmd_linked_list_type), POINTER      :: cmd_ll
    LOGICAL                                  :: at_end, explicit
    INTEGER, DIMENSION(master%n_workers)     :: last_iframe
    TYPE(swarm_report_type)                 :: report_log

    last_iframe(:) = 0
    last_command(:)%cmd_id = -1

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
      CALL swarm_comlog_read(parser, report_log, at_end, master%error)
      IF(at_end) EXIT
      CALL swarm_comlog_read(parser, cmd_log, at_end, master%error)
      IF(at_end) EXIT
      CALL swarm_master_steer(master, report_log, cmd_now)
      !TODO: maybe we should just exit the loop instead of stopping?
      IF(.NOT. commands_equal(cmd_now, cmd_log, master%iw)) STOP "wrong behaviour"

      last_command(cmd_log%worker_id) = cmd_log
      last_iframe(report_log%worker_id) = report_log%iframe
    END DO

    IF(master%iw>0) WRITE(master%iw,'(A,A)') &
       " SWARM| Reached end of communication log. Queueing last commands."

    DO i=1, master%n_workers
      !TODO send iframe along with every command
      !ALLOCATE(cmd_ll)
      !cmd_ll%cmd%worker_id = i
      !cmd_ll%cmd%cmd_id = SWARM_CMD_SETIFRAME
      !cmd_ll%cmd%iframe = last_iframe(i)
      !cmd_ll%next => master%command_queue
      !master%command_queue => cmd_ll

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
 FUNCTION commands_equal(cmd1, cmd2, iw) RESULT(res)
    TYPE(swarm_command_type)                :: cmd1, cmd2
    LOGICAL                                  :: res
    INTEGER                                  :: iw

    res =           cmd1%worker_id == cmd2%worker_id
    res = res .AND. cmd1%cmd_id == cmd2%cmd_id
    res = res .AND. ABS(cmd1%temperature - cmd2%temperature) < 1e-10
    res = res .AND. ALL(ABS(cmd1%positions -cmd2%positions) < 1e-10)

    IF(.NOT.res .AND. iw>0)THEN
      WRITE (iw,*) "commands_equal failed"
      WRITE (iw,*) "worker_id",cmd1%worker_id, cmd2%worker_id
      WRITE (iw,*) "cmd_id", cmd1%cmd_id, cmd2%cmd_id
      WRITE (iw,*) "temperature", cmd1%temperature, cmd2%temperature
      WRITE (iw,*) "cmd1%positions: ", cmd1%positions
      WRITE (iw,*) "cmd2%positions: ", cmd2%positions
    END IF
 END FUNCTION commands_equal


! *****************************************************************************
! *****************************************************************************
  FUNCTION command_queue_pop(master, worker_id) RESULT(cmd)
    TYPE(swarm_master_type)                 :: master
    INTEGER                                  :: worker_id
    TYPE(swarm_command_type)                :: cmd

    TYPE(cmd_linked_list_type), POINTER      :: cmd_ll, cmd_ll_prev

    cmd%cmd_id = -1 ! default result = nothing found

    NULLIFY(cmd_ll_prev)
    cmd_ll => master%command_queue
    DO WHILE(ASSOCIATED(cmd_ll))
       IF(cmd_ll%cmd%worker_id /= worker_id) THEN
          cmd_ll_prev => cmd_ll
          cmd_ll => cmd_ll%next
          CYCLE
       END IF

       ! we found a queued command for this worker
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
 SUBROUTINE swarm_master_steer(master, report, cmd)
    TYPE(swarm_master_type)                 :: master
    TYPE(swarm_report_type)                 :: report
    TYPE(swarm_command_type)                :: cmd

! First check if there are queued commands for this worker

    cmd = command_queue_pop(master, report%worker_id)
    IF(cmd%cmd_id > 0) THEN
       IF(master%iw>0) WRITE (master%iw,'(A,A,A,I9,X,A)') ' SWARM| ',&
         REPEAT("*", 9), " Sending out queued command to worker: ",&
         cmd%worker_id,REPEAT("*", 9)
       RETURN
    END IF

    ! There are no queued commands. Do the normal processing.
    cmd%worker_id = report%worker_id
    master%Elowest = MIN(master%Elowest, report%Epot)
    master%i_iteration = master%i_iteration + 1
    IF(master%iw>0) THEN
       WRITE (master%iw,'(A,A,X,I8,A,A)') ' SWARM| ',REPEAT("*", 13),&
          master%i_iteration,'. Global Optimisation Iteration ',REPEAT("*", 13)
       WRITE (master%iw,'(A,42X,I8)')&
              " SWARM| Reporting worker ", report%worker_id
       WRITE (master%iw,'(A,16X,E15.8)')&
              " SWARM| Reported potential Energy [Hartree] ", report%Epot
       WRITE (master%iw,'(A,9X,E15.8)')&
        " SWARM| Lowest reported potential Energy [Hartree] ",master%Elowest
    END IF

    IF(master%Elowest<master%Emin .AND. .NOT. master%should_stop) THEN
       IF(master%iw>0)  WRITE (master%iw, *)  "SWARM| Reached Emin > Epot. Quitting."
       master%should_stop = .TRUE.
    ENDIF

    IF(master%i_iteration>=master%max_iter .AND. .NOT. master%should_stop) THEN
       IF(master%iw>0)  WRITE (master%iw, *)  "SWARM| Reached MAX_ITER. Quitting."
       master%should_stop = .TRUE.
    ENDIF

    IF(.NOT. master%should_stop) THEN
      CALL external_control(master%should_stop,"SWARM",master%globenv,error=master%error)
      IF(master%should_stop .AND. master%iw>0) &
         WRITE (master%iw, *)  "SWARM| Received stop from external_control. Quitting."
    END IF

    CALL swarm_comlog_write(master%comlog_unit, master%i_iteration, report)
    CALL glbopt_master_steer(master%glbopt, report, cmd)
    CALL swarm_comlog_write(master%comlog_unit, master%i_iteration, cmd)

    IF(master%should_stop) THEN
       cmd%cmd_id = SWARM_CMD_SHUTDOWN  !overwrite command
       IF(master%iw>0) WRITE(master%iw,'(A,25X,I8)') &
         " SWARM| Sending shutdown command to worker",cmd%worker_id
    ENDIF

 END SUBROUTINE


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE swarm_master_finalize(master)
    TYPE(swarm_master_type)                 :: master

    TYPE(cp_logger_type), POINTER            :: logger

    IF(master%iw > 0) WRITE(master%iw,"(A)") " SWARM| Shutting down the master."

    CALL glbopt_master_finalize(master%glbopt)
    DEALLOCATE(master%glbopt)

    logger => cp_error_get_logger(master%error)
    CALL cp_print_key_finished_output(master%iw, logger,&
            master%swarm_section, "PRINT%MASTER_RUN_INFO", error=master%error)
    CALL cp_print_key_finished_output(master%comlog_unit, logger,&
            master%swarm_section, "PRINT%COMMUNICATION_LOG", error=master%error)

    CALL rm_timer_env() !pops the top-most timer
 END SUBROUTINE swarm_master_finalize

END MODULE swarm_master

