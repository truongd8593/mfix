CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: SET_DOLLAR (LINE,LINE_LEN)                             C
C  Purpose: Append a $ to the end of a string (for NCAR graphics)      C
C                                                                      C
C  Author: P. Nicoletti                               Date: 12-APR-92  C
C  Reviewer:                                                           C
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
C  Local variables: L, LAST_CHAR                                       C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_DOLLAR(LINE,LINE_LEN)
C
      IMPLICIT NONE
C
      CHARACTER*(*) LINE
      INTEGER       LINE_LEN, L, LAST_CHAR
C
      LAST_CHAR = LINE_LEN
      DO L = 1,LINE_LEN-1
         IF (LINE(L:L).NE.' ') LAST_CHAR = L
      END DO
C
      IF (LAST_CHAR.LT.LINE_LEN) THEN
         LINE(LAST_CHAR+1:LAST_CHAR+1) = '$'
      ELSE
         LINE(LAST_CHAR:LAST_CHAR) = '$'
      END IF
C
      RETURN
      END
