CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: GET_MU_s                                               C
C  Purpose: calculate MU_s from the RES file                           C
C                                                                      C
C  Author: P. Nicoletti                               Date: 13-MAR-92  C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: MMAX, IJKMAX2, EP_s, RUN_NAME, TIME           C
C  Variables modified: M, IJK                                          C
C                                                                      C
C  Local variables: ARRAY                                              C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_MU_s
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'visc_s.inc'
      INCLUDE 'run.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'xforms.inc'
C
C  Functions
C
      DOUBLE PRECISION G_0
C
      DOUBLE PRECISION    ARRAY(DIMENSION_3,DIMENSION_M), K_1m
      LOGICAL             INTER
      INTEGER             IER, M, IJK
C
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
C
      INTER = .FALSE.
C
      IF (DO_XFORMS) THEN
         C_E = E_PASS
      ELSE
         IF (C_E.EQ.UNDEFINED) THEN
            WRITE (*,*) ' Enter Coefficient of restitution (e) value'
            READ  (*,*) C_e
         END IF
      END IF
C
      CLOSE(UNIT_OUT)
      IF (.NOT.DO_XFORMS) CALL GET_FILE_NAME(TEMP_FILE)
      OPEN (UNIT=UNIT_OUT,FILE=TEMP_FILE,STATUS='UNKNOWN')
      CALL READ_RES1
      WRITE (UNIT_OUT,*) RUN_NAME
      WRITE (UNIT_OUT,*) ' '
      WRITE (UNIT_OUT,*) ' data from restart file'
      WRITE (UNIT_OUT,*) ' TIME = ' , TIME
C
      DO M = 1, MMAX
        CALL CALC_MU_s(M, IER)
      ENDDO
C
      DO 200 M = 1,MMAX
         IF (DO_XFORMS) THEN
            CALL CHECK_INTER(INTER)
            IF (INTER) RETURN
         END IF
         DO 100 IJK = 1,IJKMAX2
            IF (EP_s(IJK,M) .NE. 0.0) THEN
               K_1m = 2.D0 * (ONE + C_e) * RO_s(M) * G_0(IJK, M, M)
               ARRAY(IJK,M) = P_s(IJK,M) / K_1m / EP_s(IJK,M)**2
            ELSE
               ARRAY(IJK,M) = 0.0
            END IF
100      CONTINUE
200   CONTINUE
C
      DO 300 M = 1,MMAX
         IF (DO_XFORMS) THEN
            CALL CHECK_INTER(INTER)
            IF (INTER) RETURN
         END IF
         WRITE (UNIT_OUT,*) ' M    = ' , M
         WRITE (UNIT_OUT,*) ' '
         CALL OUT_ARRAY(ARRAY(1,M),'THETA')
         CALL OUT_ARRAY(MU_s(1,M),'MU_s')
300   CONTINUE
C
      CLOSE (UNIT=UNIT_OUT)
      RETURN
      END
