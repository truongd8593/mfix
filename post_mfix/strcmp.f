CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: STRCMP                                                 C
C  Purpose: Compare two strings                                        C
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
      LOGICAL FUNCTION STRCMP(STRING1, STRING2)
C
      IMPLICIT NONE
C
      CHARACTER*(*) STRING1, STRING2
      INTEGER LEN1, LEN2, L
C
      STRCMP = .FALSE.
      LEN1 = LEN(STRING1)
      LEN2 = LEN(STRING2)
      IF(LEN1 .NE. LEN2) RETURN
      DO 10 L = 1, LEN1
        IF(STRING1(L:L) .NE. STRING2(L:L))RETURN
10    CONTINUE
      STRCMP = .TRUE.
      RETURN
      END
