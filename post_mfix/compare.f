CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: COMPARE(V1, V2)                                        C
C  Purpose: Return .TRUE. if values V1 and V2 are nearly equal         C
C                                                                      C
C  Author: M. Syamlal                                 Date: 29-JUL-92  C
C  Reviewer: W. Rogers                                Date: 11-DEC-92  C
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
      LOGICAL FUNCTION COMPARE(V1, V2)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'toleranc.inc'
C
C  Local variables
C
C                      Values to be compared
      DOUBLE PRECISION V1, V2
C
      IF(ABS(V1) .LE. SMALL_NUMBER)THEN
        IF(ABS(V2) .LE. SMALL_NUMBER)THEN
          COMPARE = .TRUE.
        ELSE
          COMPARE = .FALSE.
        ENDIF
      ELSE
        IF(ABS(V2/V1 - ONE) .LE. TOL_COM)THEN
          COMPARE = .TRUE.
        ELSE
          COMPARE = .FALSE.
        ENDIF
      ENDIF
      RETURN
      END
