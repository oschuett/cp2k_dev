!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE hash_functions
  USE kinds,                           ONLY: int_4

 IMPLICIT NONE
 PRIVATE
 
 PUBLIC :: b3hs_hash_key_jenkins
 CONTAINS

! *****************************************************************************
!> \brief Take from: http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2005-10/msg00164.html
!>        Converted to take a character-array rathern than an assumed-length string as key.        
!> \par History
!> \author Ole
! *****************************************************************************
FUNCTION b3hs_hash_key_jenkins (key, range) RESULT (code)
  CHARACTER, DIMENSION(:), INTENT(in) :: key
  INTEGER, INTENT(in)      :: range
  INTEGER                  :: code


  INTEGER                  :: len_key
  INTEGER(int_4)           :: a
  INTEGER(int_4)           :: b
  INTEGER(int_4)           :: c
  INTEGER                  :: k


! Hash the key into a code, using the algorithm
! described by Bob Jenkins at:
!  http://burtleburtle.net/bob/hash/doobs.html
!
! Note that range should be a power of 2, and
! that the 32-bit algorithm is used


  !len_key = LEN_TRIM(key)
  len_key = SIZE(key)


  a = -1640531527_int_4 ! 0x9E3779B9
  b = a
  c = 305419896_int_4   ! 0x12345678


  k = 1


  char_loop : DO


     IF(len_key < 12) EXIT char_loop


! Pack the key into 32 bits


     a = a + ICHAR(key(k+0))  + ISHFT(ICHAR(key(k+1)), 8) + &
     &       ISHFT(ICHAR(key(k+2)), 16) + ISHFT(ICHAR(key(k+3)), 24)
     b = b + ICHAR(key(k+4))  + ISHFT(ICHAR(key(k+5)), 8) + &
     &       ISHFT(ICHAR(key(k+6)), 16) + ISHFT(ICHAR(key(k+7)), 24)
     c = c + ICHAR(key(k+8))  + ISHFT(ICHAR(key(k+9)), 8) + &
     &       ISHFT(ICHAR(key(k+10)), 16) + ISHFT(ICHAR(key(k+11)), 24)


! Mix it up


     CALL b3hs_hash_key_jenkins_mix_()


     k = k + 12


     len_key = len_key - 12


  END DO char_loop


  c = c + len_key


! Process remaining bits


  SELECT CASE(len_key)
  CASE(11)
     c = c + ISHFT(ICHAR(key(k+10)), 24) + ISHFT(ICHAR(key(k+9)), 16) + &
     &       ISHFT(ICHAR(key(k+8)), 8)
     b = b + ISHFT(ICHAR(key(k+7)), 24) + ISHFT(ICHAR(key(k+6)), 16) + &
     &       ISHFT(ICHAR(key(k+5)), 8) + ICHAR(key(k+4))
     a = a + ISHFT(ICHAR(key(k+3)), 24) + ISHFT(ICHAR(key(k+2)), 16) + &
     &       ISHFT(ICHAR(key(k+1)), 8) + ICHAR(key(k))
  CASE(10)
     c = c + ISHFT(ICHAR(key(k+9)), 16) + ISHFT(ICHAR(key(k+8)), 8)
     b = b + ISHFT(ICHAR(key(k+7)), 24) + ISHFT(ICHAR(key(k+6)), 16) + &
     &       ISHFT(ICHAR(key(k+5)), 8) + ICHAR(key(k+4))
     a = a + ISHFT(ICHAR(key(k+3)), 24) + ISHFT(ICHAR(key(k+2)), 16) + &
     &       ISHFT(ICHAR(key(k+1)), 8) + ICHAR(key(k))
  CASE(9)
     c = c + ISHFT(ICHAR(key(k+8)), 8)
     b = b + ISHFT(ICHAR(key(k+7)), 24) + ISHFT(ICHAR(key(k+6)), 16) + &
     &       ISHFT(ICHAR(key(k+5)), 8) + ICHAR(key(k+4))
     a = a + ISHFT(ICHAR(key(k+3)), 24) + ISHFT(ICHAR(key(k+2)), 16) + &
     &       ISHFT(ICHAR(key(k+1)), 8) + ICHAR(key(k))
  CASE(8)
     b = b + ISHFT(ICHAR(key(k+7)), 24) + ISHFT(ICHAR(key(k+6)), 16) + &
     &       ISHFT(ICHAR(key(k+5)), 8) + ICHAR(key(k+4))
     a = a + ISHFT(ICHAR(key(k+3)), 24) + ISHFT(ICHAR(key(k+2)), 16) + &
     &       ISHFT(ICHAR(key(k+1)), 8) + ICHAR(key(k))
  CASE(7)
     b = b + ISHFT(ICHAR(key(k+6)), 16) + ISHFT(ICHAR(key(k+5)), 8) + &
     &       ICHAR(key(k+4))
     a = a + ISHFT(ICHAR(key(k+3)), 24) + ISHFT(ICHAR(key(k+2)), 16) + &
     &       ISHFT(ICHAR(key(k+1)), 8) + ICHAR(key(k))
  CASE(6)
     b = b + ISHFT(ICHAR(key(k+5)), 8) + ICHAR(key(k+4))
     a = a + ISHFT(ICHAR(key(k+3)), 24) + ISHFT(ICHAR(key(k+2)), 16) + &
     &       ISHFT(ICHAR(key(k+1)), 8) + ICHAR(key(k))
  CASE(5)
     b = b + ICHAR(key(k+4))
     a = a + ISHFT(ICHAR(key(k+3)), 24) + ISHFT(ICHAR(key(k+2)), 16) + &
     &       ISHFT(ICHAR(key(k+1)), 8) + ICHAR(key(k))
  CASE(4)
     a = a + ISHFT(ICHAR(key(k+3)), 24) + ISHFT(ICHAR(key(k+2)), 16) + &
     &       ISHFT(ICHAR(key(k+1)), 8) + ICHAR(key(k))
  CASE(3)
     a = a + ISHFT(ICHAR(key(k+2)), 16) + ISHFT(ICHAR(key(k+1)), 8) + &
     &       ICHAR(key(k))
  CASE(2)
     a = a + ISHFT(ICHAR(key(k+1)), 8) + ICHAR(key(k))
  CASE(1)
     a = a + ICHAR(key(k))
  END SELECT


  CALL b3hs_hash_key_jenkins_mix_()


  code = IAND(c, range - 1) + 1


! Finish


  RETURN


CONTAINS


  SUBROUTINE b3hs_hash_key_jenkins_mix_


! Mix a, b and c


    a = IEOR(a - b - c, ISHFT(c, -13))
    b = IEOR(b - c - a, ISHFT(a, 8))
    c = IEOR(c - a - b, ISHFT(b, -13))


    a = IEOR(a - b - c, ISHFT(c, -12))
    b = IEOR(b - c - a, ISHFT(a, 16))
    c = IEOR(c - a - b, ISHFT(b, -5))


    a = IEOR(a - b - c, ISHFT(c, -3))
    b = IEOR(b - c - a, ISHFT(a, 10))
    c = IEOR(c - a - b, ISHFT(b, -15))


! Finish


    RETURN


  END SUBROUTINE b3hs_hash_key_jenkins_mix_


END FUNCTION b3hs_hash_key_jenkins

END MODULE hash_functions
