CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: MAKE_UPPER_CASE (LINE_STRING,MAXCOL)                   C
C  Purpose: change lowercase characters to uppercase in input line     C
C                                                                      C
C  Author: P.Nicoletti                                Date: 26-NOV-91  C
C  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: A_UP, A_LO, Z_LO, A_DIFF, INT_C, L                 C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE MAKE_UPPER_CASE(LINE_STRING,MAXCOL)
C
      IMPLICIT NONE
C
C passed arguments:
C
C                   input line to change to uppercase
      CHARACTER*(*) LINE_STRING
C
C                   number of characters to look at in LINE_STRING
      INTEGER       MAXCOL
C
C local variables:
C
C                   ICHAR value for UPPERCASE A
      INTEGER       A_UP
C
C                   ICHAR value for lowercase a
      INTEGER       A_LO
C
C                   ICHAR value for lowercase z
      INTEGER       Z_LO
C
C                   ICHAR differnce between lower and uppercase letters
      INTEGER       A_DIFF
C
C                   holds ICHAR value of current character
      INTEGER       INT_C
C
C                   loop index
      INTEGER       L
C
      A_UP = ICHAR('A')
      A_LO = ICHAR('a')
      Z_LO = ICHAR('z')
      A_DIFF = A_LO - A_UP
C
      DO 100 L = 1,MAXCOL
         INT_C = ICHAR(LINE_STRING(L:L))
         IF (A_LO.LE.INT_C .AND. INT_C.LE.Z_LO) THEN
            INT_C = INT_C - A_DIFF
            LINE_STRING(L:L) = CHAR(INT_C)
         END IF
100   CONTINUE
C
      RETURN
      END
