CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: STREQS(STRING1, STRING2)                               C
C  Purpose: Transfer string2 to string1                                C
C                                                                      C
C  Author: M. Syamlal                                 Date: 03-NOV-93  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE STREQS(STRING1, STRING2)
C
      IMPLICIT NONE
C
      CHARACTER*(*) STRING1, STRING2
      INTEGER LEN1, L
C
      LEN1 = MIN(LEN(STRING1), LEN(STRING2))
      DO 10 L = 1, LEN1
        STRING1(L:L) = STRING2(L:L)
10    CONTINUE
      RETURN
      END
