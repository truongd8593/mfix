CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: GET_INDEX(A, ARRAY, LMAX, EXT, L1,AC)                  C
C  Purpose: Get the index for A from ARRAY which is dimensioned LMAX   C
C           An  index value of 0 indicates error                       C
C                                                                      C
C  Author: M. Syamlal                                 Date: 06-DEC-93  C
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
      INTEGER FUNCTION GET_INDEX(A, ARRAY, LMAX, EXT, L1, AC)
C
      IMPLICIT NONE
C
      INTEGER L, LMAX, L1
      REAL A, ARRAY(LMAX), DIFF_LAST
      LOGICAL EXT
      CHARACTER ANS
      CHARACTER*3 AC
C
      GET_INDEX = 0
      L1        = 0
      IF( A .LT. ARRAY(1) .OR. A .GT. ARRAY(LMAX)) THEN
        IF(.NOT. EXT)THEN
          WRITE(*, '(A, A3, A, G12.5)')' Could not find ', AC, ' = ', A
          WRITE(*, '(A,$)')' Enter Y to extrapolate > '
          READ(*,'(A)')ANS
          IF(ANS .EQ. 'Y' .OR. ANS .EQ. 'y')THEN
            EXT = .TRUE.
          ELSE
            STOP
          ENDIF
        ENDIF
        IF(A .LT. ARRAY(1)) THEN
          IF(LMAX .GE. 2) THEN
            GET_INDEX = 2
            L1        = 2
          ELSE
            GET_INDEX = 1
            L1        = 1
          ENDIF
        ELSEIF( A .GT. ARRAY(LMAX))THEN
          IF(LMAX .GE. 2) THEN
            GET_INDEX = LMAX - 1
            L1        = LMAX - 1
          ELSE
            GET_INDEX = 1
            L1        = 1
          ENDIF
        ENDIF
      ENDIF
      DIFF_LAST = 1E32
      DO 10 L = 1, LMAX
        IF(ARRAY(L) .GE. A) THEN
          IF( ABS(ARRAY(L) - A) .LE. DIFF_LAST) THEN
            GET_INDEX = L
            L1        = MAX((L - 1), 1)
          ELSE
            GET_INDEX = MAX((L - 1), 1)
            L1        = L
          ENDIF
          RETURN
        ELSE
          DIFF_LAST = ABS(ARRAY(L) - A)
        ENDIF
10    CONTINUE
C
      RETURN
      END
