CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: GET_GRANULAR_QTY                                       C
C  Purpose: calculate MU_s/EP_s from the RES file                      C
C                                                                      C
C  Author: M. Syamlal                                 Date: 24-AUG-94  C
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
C  Local variables: FILE_NAME, ARRAY                                   C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_GRANULAR_QTY
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
C  Local variables
C
      INTEGER      MM, M, IER, IJK
C
C                file version ID
      CHARACTER  VERSION*512
      INTEGER     REC_POINTER(N_SPX), L, NEXT_REC, NUM_REC, NSTEP_1
      LOGICAL     READ_SPX(N_SPX), AT_EOF(N_SPX)
      CHARACTER   ANSWER
      REAL        TIME_REAL(N_SPX), TIME_NOW
      REAL        TIME_IN_RES
      DOUBLE PRECISION K_1m, Tmp(DIMENSION_3)
C
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
C
      CALL READ_RES1
      TIME_IN_RES = TIME
C
      IF (DO_XFORMS) THEN
         C_E = E_PASS
      ELSE
         IF (C_E .EQ. UNDEFINED) THEN
            WRITE (*,'(A,$)')
     &        ' Enter Coefficient of restitution value (e): '
            READ  (*,*) C_e
         END IF
      END IF
C
      IF(.NOT.DO_XFORMS)THEN
         CALC_GT = .FALSE.
         CALC_GV = .FALSE.
         WRITE (*,'(A,$)')
     &     ' Do you need granular temperature ? (Y/N) '
         READ (*,'(A)') ANSWER
         IF(ANSWER .EQ. 'Y' .OR. ANSWER .EQ. 'y') CALC_GT = .TRUE.
      ENDIF
      IF (CALC_GT) THEN
        IF (.NOT.DO_XFORMS) CALL GET_FILE_NAME(TEMP_FILE2)
        OPEN (UNIT=70,FILE=TEMP_FILE2,STATUS='NEW',RECL=128,
     &      ACCESS='DIRECT',FORM='UNFORMATTED',ERR=1000)
C
        VERSION = 'SP1 = 01.00'
        WRITE (70,REC=1) VERSION
        WRITE (70,REC=2)RUN_NAME,ID_MONTH,ID_DAY,ID_YEAR,
     &          ID_HOUR,ID_MINUTE,ID_SECOND
        WRITE (70,REC=3) 4, -1
      ENDIF
C
      IF(.NOT.DO_XFORMS)THEN
         CALC_GV = .FALSE.
         WRITE (*,'(A,$)')
     &     ' Do you need granular viscosity ? (Y/N) '
         READ (*,'(A)') ANSWER
         IF(ANSWER .EQ. 'Y' .OR. ANSWER .EQ. 'y') CALC_GV = .TRUE.
      ENDIF
      IF (CALC_GV) THEN
        IF (.NOT.DO_XFORMS) CALL GET_FILE_NAME(TEMP_FILE)
        OPEN (UNIT=71,FILE=TEMP_FILE,STATUS='NEW',RECL=128,
     &      ACCESS='DIRECT',FORM='UNFORMATTED',ERR=1000)
C
        VERSION = 'SP1 = 01.00'
        WRITE (71,REC=1) VERSION
        WRITE (71,REC=2)RUN_NAME,ID_MONTH,ID_DAY,ID_YEAR,
     &          ID_HOUR,ID_MINUTE,ID_SECOND
        WRITE (71,REC=3) 4, -1
      ENDIF
C
      IF(CALC_GT .OR. CALC_GV) THEN
        IF(MMAX .GT. 1) THEN
          IF (.NOT.DO_XFORMS) THEN
             WRITE(*, '(A,$)')' Enter solids phase number: '
             READ(*, *)M_USE
          END IF
        ELSE
          M_USE = 1
        ENDIF
      ELSE
        RETURN
      ENDIF
      MM = M_USE
C
10    IF (.NOT.DO_XFORMS) THEN
         WRITE (*,'(A,$)') ' Enter TIME_START and TIME_END: '
         READ (*, *)TIME_START, TIME_END
      END IF
C
      DO L = 1, N_SPX
         REC_POINTER(L) = 4
      END DO
C
      DO L = 1, N_SPX
         READ_SPX(L)    = .FALSE.
         AT_EOF(L)      = .FALSE.
      END DO
      READ_SPX(1)       = .TRUE.
      READ_SPX(4)       = .TRUE.
      IF(MMAX .GT. 1)READ_SPX(5) = .TRUE.
C
      IF(TIME_START .LT. TIME_IN_RES) THEN
        CALL SEEK_TIME(READ_SPX, TIME_START, REC_POINTER, TIME_NOW)
        IF(TIME_NOW .LT. ZERO) THEN
          WRITE(*,*)' Could not find record for TIME_START'
          GOTO 10
        ENDIF
      ENDIF
C
100   CONTINUE
      IF(TIME_START .LT. TIME_IN_RES) THEN
        CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,
     &                    TIME_NOW, TIME_REAL, NSTEP_1)
      ELSE
        CALL READ_RES1
        TIME_NOW = TIME_IN_RES
      ENDIF
C
      IF (TIME_NOW .LT. ZERO) GOTO 1000
      IF (TIME_NOW .LT. TIME_START) GOTO 100
      IF(MMAX .EQ. 1) THEN
        DO IJK = 1, IJKMAX2
          ROP_s(IJK, 1) = RO_s(1) * (1. - EP_g(IJK))
        END DO
      ENDIF
C
      DO M = 1, MMAX
        CALL CALC_MU_s(M, IER)
      ENDDO
c
      TIME = TIME_NOW
      NSTEP = NSTEP_1
      IF(CALC_GT) THEN
         DO IJK = 1, IJKMAX2
           IF( FLUID_AT(IJK) ) THEN
             K_1m = 2.D0 * (ONE + C_e) * RO_s(MM) * G_0(IJK, MM, MM)
             IF(K_1m .NE. 0.0 .AND. EP_s(IJK,MM) .NE. 0.0) THEN
               Tmp(IJK) = P_s(IJK,MM)/K_1m/EP_s(IJK, MM)**2
             ELSE
               Tmp(IJK) = 0.0
             ENDIF
           ELSE
             Tmp(IJK) = 0.0
           ENDIF
         END DO
         READ  (70,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (70,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL OUT_BIN_R(70,Tmp,IJKMAX2,NEXT_REC)
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(70,REC=3) NEXT_REC, NUM_REC
      ENDIF
C
      IF(CALC_GV) THEN
         DO IJK = 1, IJKMAX2
           IF( FLUID_AT(IJK) ) THEN
             Tmp(IJK) = MU_s(IJK, MM)
           ELSE
             Tmp(IJK) = 0.0
           ENDIF
         END DO
         READ  (71,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (71,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL OUT_BIN_R(71,MU_s(1,MM),IJKMAX2,NEXT_REC)
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(71,REC=3) NEXT_REC, NUM_REC
      ENDIF
C
      IF(TIME_START .LT. TIME_IN_RES) GOTO 100
C
1000  CLOSE (UNIT=70)
      CLOSE (UNIT=71)
      RETURN
      END
