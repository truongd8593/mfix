CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,          C
C                              TIME_NOW, TIME_REAL, NSTEP_1)           C
C                                                                      C
C  Purpose: Synchronize the files enabled by READ_SPX                  C
C                                                                      C
C  Author: M. Syamlal                                 Date: 27-OCT-93  C
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
C  Local variables: TDIFF                                              C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_SAME_TIME (READ_SPX, REC_POINTER,
     &                          AT_EOF, TIME_NOW, TIME_REAL, NSTEP_1)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
C
      LOGICAL COMPARE
C
      INTEGER   REC_POINTER(*) , NSTEP_1
      INTEGER   L
      LOGICAL   READ_SPX(*) , READ_SPX_STORE(N_SPX), AT_EOF(*)
      REAL      TDIFF , TIME_REAL(*), TIME_NOW
C
      TIME_NOW = -ONE
C
C  Store READ_SPX array
C
      DO 10 L = 1, N_SPX
        READ_SPX_STORE(L) = READ_SPX(L)
10    CONTINUE
C
100   CONTINUE
      CALL READ_SPX1(READ_SPX,REC_POINTER,AT_EOF, TIME_REAL, NSTEP_1)
C
C  Compare times
C
      DO 120 L = 1, N_SPX
        IF(READ_SPX(L)) THEN
          IF(AT_EOF(L)) THEN
            TIME_NOW = -ONE
            RETURN
          ENDIF
          TIME_NOW = MAX(TIME_NOW, TIME_REAL(L))
        ENDIF
120   CONTINUE
C
      DO 140 L = 1, N_SPX
        IF(READ_SPX(L)) THEN
          IF( COMPARE( DBLE(TIME_REAL(L)), DBLE(TIME_NOW))) THEN
            READ_SPX(L) = .FALSE.
          ENDIF
        ENDIF
140   CONTINUE
C
      DO 160 L = 1, N_SPX
        IF(READ_SPX(L)) GOTO 100
160   CONTINUE
C
C  Restore READ_SPX array
C
      DO 500 L = 1, N_SPX
        READ_SPX(L) = READ_SPX_STORE(L)
500   CONTINUE
C
      RETURN
C
      END
