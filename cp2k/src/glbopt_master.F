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
   INTEGER                                             :: translog_unit
   TYPE(cp_error_type)                                 :: error
   TYPE(section_vals_type), POINTER                    :: glbopt_section => Null()
   TYPE(cp_para_env_type), POINTER                     :: para_env => Null()
 END TYPE glbopt_master_type

 CONTAINS



! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_init(master, para_env, root_section, n_walkers, error)
    TYPE(glbopt_master_type)                 :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
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
    master%iw = cp_print_key_unit_nr(logger, master%glbopt_section,&
          "PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    master%translog_unit = cp_print_key_unit_nr(logger,master%glbopt_section,"PRINT%TRANSACTION_LOG",&
            !extension=my_ext, file_position=my_pos, file_action=my_act,&
            !file_form=my_form, middle_name=TRIM(my_middle), is_new_file=new_file,&
            extension=".xyz", error=error)

    CALL section_vals_val_get(master%glbopt_section,"E_MIN", r_val=master%Emin, error=error)
    CALL section_vals_val_get(master%glbopt_section,"MAX_ITER", i_val=master%max_iter, error=error)

    ALLOCATE(master%minhop)
    CALL glbopt_minhop_init(master%minhop, master%glbopt_section, n_walkers, master%iw, error)

    CALL do_replay(master)

 END SUBROUTINE glbopt_master_init







! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_steer(master, report, cmd)
    TYPE(glbopt_master_type)   :: master
    TYPE(glbopt_report_type)   :: report
    TYPE(glbopt_command_type)  :: cmd

    master%Elowest = MIN(master%Elowest, report%Epot)
    master%i_iteration = master%i_iteration + 1
    IF(master%iw>0) WRITE(master%iw,'(A,I3,A,E15.8,A,E15.8)') &
       " GLOBAL_OPT| Walker ", report%walker_id, " reported minima with Epot=",&
       report%Epot," Elowest=",master%Elowest

    IF(master%Elowest<master%Emin .AND. .NOT. master%should_stop) THEN
       CALL print(master, "Reached Emin > Epot. Quitting.")
       master%should_stop = .TRUE.
    ENDIF

    IF(master%i_iteration>=master%max_iter .AND. .NOT. master%should_stop) THEN
       CALL print(master, "Reached MAX_ITER. Quitting.")
       master%should_stop = .TRUE.
    ENDIF

    CALL translog_write_report(master, report)
    CALL glbopt_minhop_steer(master%minhop, report, cmd)
    CALL translog_write_command(master, cmd)

    IF(master%should_stop) THEN
       !CALL print(master, "Sending shutdown command")
       cmd%cmd_id = GLBOPT_CMD_SHUTDOWN  !overwrite command
    ENDIF

 END SUBROUTINE


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE glbopt_master_finalize(master)
    TYPE(glbopt_master_type)                 :: master
    TYPE(cp_logger_type), POINTER            :: logger

    CALL glbopt_minhop_finalize(master%minhop)
    DEALLOCATE(master%minhop)

    logger => cp_error_get_logger(master%error)
    CALL cp_print_key_finished_output(master%iw, logger,&
            master%glbopt_section, "PRINT%MASTER_RUN_INFO", error=master%error)
    CALL cp_print_key_finished_output(master%translog_unit, logger,&
            master%glbopt_section, "PRINT%TRANSACTION_LOG", error=master%error)

    CALL rm_timer_env() !pops the top-most timer
 END SUBROUTINE glbopt_master_finalize

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE do_replay(master)
    TYPE(glbopt_master_type)                 :: master
    CHARACTER(LEN=default_path_length)       :: filename
    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL                                  :: explicit
    TYPE(glbopt_report_type)                 :: report_log
    TYPE(glbopt_command_type)                :: cmd_log, cmd_now
    LOGICAL                                  :: at_end
    NULLIFY(parser)

    ! Initialize parser for trajectory
    CALL section_vals_val_get(master%glbopt_section,"REPLAY_TRANSACTION_LOG",&
       c_val=filename, explicit=explicit, error=master%error)

    IF(.NOT. explicit) RETURN
    IF(master%iw>0) WRITE(master%iw,'(A,A)') &
       " GLOBAL_OPT| Starting replay of transaction-log: ", TRIM(filename)

    CALL parser_create(parser, filename, para_env=master%para_env, error=master%error)

    at_end = .FALSE.
    DO
       CALL translog_read_report(master, parser, report_log, at_end)
       IF(at_end) EXIT
       CALL translog_read_command(master, parser, cmd_log, at_end)
       IF(at_end) EXIT
       CALL glbopt_master_steer(master, report_log, cmd_now)
       IF(.NOT. commands_equal(cmd_now, cmd_log)) STOP "wrong behaviour"
    END DO

    IF(master%iw>0) WRITE(master%iw,'(A,A)') &
       " GLOBAL_OPT| Replay completed. Starting normal operations."

    !  store last command for each walker
    !  send out commands for each walker
    !  evtl wail for 5 sec command (come_back_later)
    CALL parser_release(parser, error=master%error)
 END SUBROUTINE do_replay

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
 SUBROUTINE translog_write_report(master, report)
    TYPE(glbopt_master_type)   :: master
    TYPE(glbopt_report_type)   :: report
    INTEGER :: i

    IF(master%translog_unit <= 0) RETURN

    WRITE(master%translog_unit,"(I8)") SIZE(report%positions)/3 ! number of atoms
    WRITE(master%translog_unit,"(A,I8,A,A,I8,A,F20.10)")&
        'i = ',master%i_iteration,&
        ', type = "glbopt_report"',&
        ", walker_id = ",report%walker_id,&
        ", Epot = ",report%Epot! info line, free text

    DO i=1, SIZE(report%positions), 3
      WRITE(master%translog_unit,"(3F20.10)") report%positions(i:i+2)
    END DO


 END SUBROUTINE

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE translog_read_report(master, parser, report, at_end)
    TYPE(glbopt_master_type)                 :: master
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(glbopt_report_type), INTENT(OUT)    :: report
    LOGICAL                                  :: at_end

    INTEGER                                  :: natoms, i, val1
    CHARACTER(LEN=default_string_length)     :: key1, key2, key3, key4, val2, eq

    CALL parser_get_next_line(parser, 1, at_end, master%error)
    IF(at_end .OR. LEN_TRIM(parser%input_line)==0) RETURN
    CALL parser_get_object(parser, natoms, error=master%error)

    CALL parser_get_next_line(parser, 1, at_end, master%error)
    IF(at_end .OR. LEN_TRIM(parser%input_line)==0) RETURN

    READ (parser%input_line,*) key1, eq, val1, key2, eq, val2, key3, eq, report%walker_id, key4, eq, report%Epot
    IF(TRIM(key1)/="i" .OR. TRIM(key2)/="type"&
     .OR. TRIM(val2)/="glbopt_report" .OR. TRIM(key3)/="walker_id"&
     .OR. TRIM(key4)/="Epot") STOP "translog_read_report failed"

    ALLOCATE(report%positions(3*natoms))
    report%positions = 0
    DO i=0, natoms-1
        CALL parser_get_next_line(parser, 1, at_end, master%error)
        IF(at_end .OR. LEN_TRIM(parser%input_line)==0) RETURN
        READ (parser%input_line,*) report%positions(3*i + 1),&
                                   report%positions(3*i + 2),&
                                   report%positions(3*i + 3)
     END DO
 END SUBROUTINE translog_read_report

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE translog_write_command(master, cmd)
    TYPE(glbopt_master_type)                 :: master
    TYPE(glbopt_command_type)                :: cmd
    INTEGER :: i

    IF(master%translog_unit <= 0) RETURN

    WRITE(master%translog_unit,"(I8)") SIZE(cmd%positions)/3 ! number of atoms
    WRITE(master%translog_unit,"(A,I8,A,A,I8,A,I3,A,F20.10)")&
        'i = ',master%i_iteration,&
        ', type = "glbopt_command"',&
        ", walker_id = ",cmd%walker_id,&
        ", cmd_id = ",cmd%cmd_id,&
        ", temperature = ",cmd%temperature! info line, free text

    DO i=1, SIZE(cmd%positions), 3
      WRITE(master%translog_unit,"(3F20.10)") cmd%positions(i:i+2)
    END DO

 END SUBROUTINE

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE translog_read_command(master, parser, cmd, at_end)
    TYPE(glbopt_master_type)                 :: master
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(glbopt_command_type), INTENT(OUT)   :: cmd
    LOGICAL                                  :: at_end

    INTEGER                                  :: natoms, i, val1
    CHARACTER(LEN=default_string_length)     :: key1, key2, key3, key4, key5, val2, eq

    CALL parser_get_next_line(parser, 1, at_end, master%error)
    IF(at_end .OR. LEN_TRIM(parser%input_line)==0) RETURN
    CALL parser_get_object(parser,natoms, error=master%error)

    CALL parser_get_next_line(parser, 1, at_end, master%error)
    IF(at_end .OR. LEN_TRIM(parser%input_line)==0) RETURN
    READ (parser%input_line,*) key1, eq, val1, key2, eq, val2,&
      key3, eq, cmd%walker_id, key4, eq, cmd%cmd_id, key5, eq, cmd%temperature

    IF(TRIM(key1)/="i" .OR. TRIM(key2)/="type" .OR. TRIM(val2)/="glbopt_command" &
      .OR. TRIM(key3)/="walker_id" .OR. TRIM(key4)/="cmd_id" &
      .OR. TRIM(key5)/="temperature") STOP "translog_read_command failed"

    ALLOCATE(cmd%positions(3*natoms))
    cmd%positions = 0
    DO i=0, natoms-1
        CALL parser_get_next_line(parser, 1, at_end, master%error)
        IF(at_end .OR. LEN_TRIM(parser%input_line)==0) RETURN
        READ (parser%input_line,*) cmd%positions(3*i + 1),&
                                   cmd%positions(3*i + 2),&
                                   cmd%positions(3*i + 3)
     END DO
 END SUBROUTINE translog_read_command

! *****************************************************************************
! *****************************************************************************
   SUBROUTINE print(master, message)
      TYPE(glbopt_master_type), INTENT(INOUT)         :: master
      CHARACTER(LEN=*)                         :: message
      IF (master%iw>0) &
        WRITE (master%iw,*) "GLOBAL_OPT| "//TRIM(message)
   END SUBROUTINE print

END MODULE glbopt_master

