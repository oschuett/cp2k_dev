!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_comlog
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

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_comlog'


 PUBLIC :: glbopt_comlog_write,  glbopt_comlog_read

  INTERFACE glbopt_comlog_write
     MODULE PROCEDURE write_report, write_command
  END INTERFACE

  INTERFACE glbopt_comlog_read
     MODULE PROCEDURE read_report, read_command
  END INTERFACE
 CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE write_report(unit, i_iteration, report)
    INTEGER                    :: unit
    INTEGER                    :: i_iteration
    TYPE(glbopt_report_type)   :: report
    INTEGER :: i

    IF(unit <= 0) RETURN

    WRITE(unit,"(I8)") SIZE(report%positions)/3 ! number of atoms
    WRITE(unit,"(A,I8,A,A,I8,A,I8A,F20.10)")&
        'i = ',i_iteration,&
        ', type = "glbopt_report"',&
        ", walker_id = ",report%walker_id,&
        ", iframe = ",report%iframe,&
        ", Epot = ",report%Epot! info line, free text

    DO i=1, SIZE(report%positions), 3
      WRITE(unit,"(3F20.10)") report%positions(i:i+2)
    END DO


 END SUBROUTINE write_report

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE read_report(parser, report, at_end, error)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(glbopt_report_type), INTENT(OUT)    :: report
    LOGICAL, INTENT(INOUT)                   :: at_end
    TYPE(cp_error_type)                                 :: error

    INTEGER                                  :: natoms, i, val1
    CHARACTER(LEN=default_string_length)     :: key1, key2, key3, key4, key5, val2, eq

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN
    CALL parser_get_object(parser, natoms, error=error)

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN

    READ (parser%input_line,*) key1, eq, val1, key2, eq, val2,&
      key3, eq, report%walker_id, key4, eq, report%iframe, key5, eq, report%Epot
    IF(TRIM(key1)/="i" .OR. TRIM(key2)/="type"&
     .OR. TRIM(val2)/="glbopt_report" .OR. TRIM(key3)/="walker_id"&
     .OR. TRIM(key4)/="iframe" .OR. TRIM(key5)/="Epot")&
        STOP "comlog_read_report failed"

    ALLOCATE(report%positions(3*natoms))
    report%positions = 0
    DO i=0, natoms-1
        CALL parser_get_next_line(parser, 1, at_end, error)
        at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
        IF(LEN_TRIM(parser%input_line) < 60) THEN
           WRITE(*,*) "LEN_TRIM(parser%input_line)=",LEN_TRIM(parser%input_line)
           WRITE(*,*) parser%input_line
        END IF
        IF(at_end) RETURN
        READ (parser%input_line,*) report%positions(3*i + 1),&
                                   report%positions(3*i + 2),&
                                   report%positions(3*i + 3)
     END DO
 END SUBROUTINE read_report

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE write_command(unit, i_iteration, cmd)
    INTEGER                                  :: unit
    INTEGER                                  :: i_iteration
    TYPE(glbopt_command_type)                :: cmd
    INTEGER :: i

    IF(unit <= 0) RETURN

    WRITE(unit,"(I8)") SIZE(cmd%positions)/3 ! number of atoms
    WRITE(unit,"(A,I8,A,A,I8,A,I3,A,I8,A,F20.10)")&
        'i = ',i_iteration,&
        ', type = "glbopt_command"',&
        ", walker_id = ",cmd%walker_id,&
        ", cmd_id = ",cmd%cmd_id,&
        ", iframe = ",cmd%iframe,&
        ", temperature = ",cmd%temperature! info line, free text

    DO i=1, SIZE(cmd%positions), 3
      WRITE(unit,"(3F20.10)") cmd%positions(i:i+2)
    END DO

 END SUBROUTINE write_command

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE read_command(parser, cmd, at_end, error)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(glbopt_command_type), INTENT(OUT)   :: cmd
    LOGICAL, INTENT(INOUT)                   :: at_end
    TYPE(cp_error_type)                      :: error

    INTEGER                                  :: natoms, i, val1
    CHARACTER(LEN=default_string_length)     :: key1, key2, key3, key4, key5, key6, val2, eq

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN
    CALL parser_get_object(parser,natoms, error=error)

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
    IF(at_end) RETURN
    READ (parser%input_line,*) key1, eq, val1, key2, eq, val2,&
      key3, eq, cmd%walker_id, key4, eq, cmd%cmd_id,&
      key5, eq, cmd%iframe, key6, eq, cmd%temperature

    IF(TRIM(key1)/="i" .OR. TRIM(key2)/="type" .OR. TRIM(val2)/="glbopt_command" &
      .OR. TRIM(key3)/="walker_id" .OR. TRIM(key4)/="cmd_id" &
      .OR. TRIM(key5)/="iframe".OR. TRIM(key6)/="temperature") &
         STOP "comlog_read_command failed"

    ALLOCATE(cmd%positions(3*natoms))
    cmd%positions = 0
    DO i=0, natoms-1
        CALL parser_get_next_line(parser, 1, at_end, error)
        at_end = at_end .OR. LEN_TRIM(parser%input_line)==0
        IF(at_end) RETURN
        READ (parser%input_line,*) cmd%positions(3*i + 1),&
                                   cmd%positions(3*i + 2),&
                                   cmd%positions(3*i + 3)
     END DO
 END SUBROUTINE read_command

END MODULE glbopt_comlog

