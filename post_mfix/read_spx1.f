CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: READ_SPX1(READ_SPX1,REC_POINTER,AT_EOF,                C
C                         TIME_REAL, NSTEP_1)                          c
C  Purpose: read in the time-dependent restart records (REAL)          C
C                                                                      C
C  Author: P. Nicoletti                               Date: 13-DEC-91  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: TIME, NSTEP, EP_g, RO_g, P_g, P_star, U_g       C
C                        V_g, W_g, U_s, V_s, W_s, ROP_s, T_g, T_s1     C
C                        T_s2, IJKMAX2, MMAX                           C
C                                                                      C
C  Local variables:  LC, NEXT_REC, TIME_REAL                           C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE READ_SPX1(READ_SPX,REC_POINTER,AT_EOF,
     &                     TIME_REAL,NSTEP_1)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'run.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'post3d.inc'
C
C passed arguments
C
C             flag whether to read a particular SPx file this time step
      LOGICAL READ_SPX(*)
C
C             pointers to next record to read in each file
      INTEGER REC_POINTER(*)
C
      LOGICAL AT_EOF(*)
      INTEGER NSTEP_1
      REAL    TIME_REAL(*)

C
C                      Dummy variable for reading T_s2
      DOUBLE PRECISION Tmp(DIMENSION_3)
C
C local variables
C
C             loop counters
      INTEGER LC,N
C
C             Pointer to the next record
      INTEGER NEXT_REC
C
C
C ".SP1" FILE         EP_g    [ ROP_g , RO_g must be calculated ...
C                                        not written out ]
C
      IF (READ_SPX(1).AND..NOT.AT_EOF(1)) THEN
         IF(.NOT.SPX_OPEN(1)) THEN
           WRITE(*,*)' SP1 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(1)
         IF (NEXT_REC.EQ.LAST_REC(1)) THEN
            AT_EOF(1) = .TRUE.
            RETURN
         END IF
         AT_EOF(1) = .FALSE.
         READ (UNIT_SPX+1,REC=NEXT_REC) TIME_REAL(1) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL IN_BIN_R(UNIT_SPX+1,EP_g,IJKMAX2,NEXT_REC)
         REC_POINTER(1) = NEXT_REC
         NSTEP_1 = NSTEP
      END IF
C
C ".SP2" FILE         P_g , P_star
C
      IF (READ_SPX(2).AND..NOT.AT_EOF(2)) THEN
         IF(.NOT.SPX_OPEN(2)) THEN
           WRITE(*,*)' SP2 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(2)
         IF (NEXT_REC.EQ.LAST_REC(2)) THEN
            AT_EOF(2) = .TRUE.
            RETURN
         END IF
         AT_EOF(2) = .FALSE.
         READ (UNIT_SPX+2,REC=NEXT_REC) TIME_REAL(2) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL IN_BIN_R(UNIT_SPX+2,P_g,IJKMAX2,NEXT_REC)
         CALL IN_BIN_R(UNIT_SPX+2,P_star,IJKMAX2,NEXT_REC)
         REC_POINTER(2) = NEXT_REC
      END IF
C
C ".SP3" FILE         U_g , V_g , W_g
C
      IF (READ_SPX(3).AND..NOT.AT_EOF(3)) THEN
         IF(.NOT.SPX_OPEN(3)) THEN
           WRITE(*,*)' SP3 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(3)
         IF (NEXT_REC.EQ.LAST_REC(3)) THEN
            AT_EOF(3) = .TRUE.
            RETURN
         END IF
         AT_EOF(3) = .FALSE.
         READ (UNIT_SPX+3,REC=NEXT_REC) TIME_REAL(3) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL IN_BIN_R(UNIT_SPX+3,U_g,IJKMAX2,NEXT_REC)
         CALL IN_BIN_R(UNIT_SPX+3,V_g,IJKMAX2,NEXT_REC)
         CALL IN_BIN_R(UNIT_SPX+3,W_g,IJKMAX2,NEXT_REC)
         REC_POINTER(3) = NEXT_REC
      END IF
C
C ".SP4" FILE         U_s , V_s , W_s
C
      IF (READ_SPX(4).AND..NOT.AT_EOF(4)) THEN
         IF(.NOT.SPX_OPEN(4)) THEN
           WRITE(*,*)' SP4 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(4)
         IF (NEXT_REC.EQ.LAST_REC(4)) THEN
            AT_EOF(4) = .TRUE.
            RETURN
         END IF
         AT_EOF(4) = .FALSE.
         READ (UNIT_SPX+4,REC=NEXT_REC) TIME_REAL(4) , NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 100 LC = 1,MMAX
            CALL IN_BIN_R(UNIT_SPX+4,U_s(1,LC),IJKMAX2,NEXT_REC)
            CALL IN_BIN_R(UNIT_SPX+4,V_s(1,LC),IJKMAX2,NEXT_REC)
            CALL IN_BIN_R(UNIT_SPX+4,W_s(1,LC),IJKMAX2,NEXT_REC)
100      CONTINUE
         REC_POINTER(4) = NEXT_REC
      END IF
C
C ".SP5" FILE         ROP_s
C
      IF (READ_SPX(5).AND..NOT.AT_EOF(5)) THEN
         IF(.NOT.SPX_OPEN(5)) THEN
           WRITE(*,*)' SP5 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(5)
         IF (NEXT_REC.EQ.LAST_REC(5)) THEN
            AT_EOF(5) = .TRUE.
            RETURN
         END IF
         AT_EOF(5) = .FALSE.
         READ (UNIT_SPX+5,REC=NEXT_REC) TIME_REAL(5) , NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 200 LC = 1,MMAX
            CALL IN_BIN_R(UNIT_SPX+5,ROP_s(1,LC),IJKMAX2,NEXT_REC)
200      CONTINUE
         REC_POINTER(5) = NEXT_REC
      END IF
C
C ".SP6" FILE         T_g  , T_s1 , T_s2
C
      IF (READ_SPX(6).AND..NOT.AT_EOF(6)) THEN
         IF(.NOT.SPX_OPEN(6)) THEN
           WRITE(*,*)' SP6 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(6)
         IF (NEXT_REC.EQ.LAST_REC(6)) THEN
            AT_EOF(6) = .TRUE.
            RETURN
         END IF
         AT_EOF(6) = .FALSE.
         READ (UNIT_SPX+6,REC=NEXT_REC) TIME_REAL(6) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL IN_BIN_R(UNIT_SPX+6,T_g,IJKMAX2,NEXT_REC)
         IF(VERSION_NUMBER .LE. 1.)THEN
           CALL IN_BIN_R (UNIT_SPX+6,T_s(1,1),IJKMAX2,NEXT_REC)
           IF(MMAX .GE. 2) THEN
             CALL IN_BIN_R (UNIT_SPX+6,T_s(1,2),IJKMAX2,NEXT_REC)
           ELSE
             CALL IN_BIN_R (UNIT_SPX+6,Tmp,IJKMAX2,NEXT_REC)
           ENDIF
         ELSE
           DO 220 LC = 1,MMAX
             CALL IN_BIN_R(UNIT_SPX+6,T_s(1,LC),IJKMAX2,NEXT_REC)
220        CONTINUE
         ENDIF
         REC_POINTER(6) = NEXT_REC
      END IF
C
C
C ".SP7" FILE         X_g, X_s
C
      IF (READ_SPX(7).AND..NOT.AT_EOF(7)) THEN
         IF(.NOT.SPX_OPEN(7)) THEN
           WRITE(*,*)' SP7 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(7)
         IF (NEXT_REC.EQ.LAST_REC(7)) THEN
            AT_EOF(7) = .TRUE.
            RETURN
         END IF
         AT_EOF(7) = .FALSE.
         READ (UNIT_SPX+7,REC=NEXT_REC) TIME_REAL(7) , NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 250 N = 1, NMAX(0)
           CALL IN_BIN_R(UNIT_SPX+7,X_g(1,N),IJKMAX2,NEXT_REC)
250      CONTINUE
         DO 300 LC = 1,MMAX
           DO 270 N = 1, NMAX(LC)
             CALL IN_BIN_R(UNIT_SPX+7,X_s(1,LC, N),IJKMAX2,NEXT_REC)
270        CONTINUE
300      CONTINUE
         REC_POINTER(7) = NEXT_REC
      END IF
C
C ".SP8" FILE         THETA_m
C
      IF (READ_SPX(8).AND..NOT.AT_EOF(8)) THEN
         IF(.NOT.SPX_OPEN(8)) THEN
           WRITE(*,*)' SP8 file is not open'
           STOP
         ENDIF
         NEXT_REC = REC_POINTER(8)
         IF (NEXT_REC.EQ.LAST_REC(8)) THEN
            AT_EOF(8) = .TRUE.
            RETURN
         END IF
         AT_EOF(8) = .FALSE.
         READ (UNIT_SPX+8,REC=NEXT_REC) TIME_REAL(8),NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 400 LC = 1,MMAX
            CALL IN_BIN_R(UNIT_SPX+8,THETA_m(1,LC),IJKMAX2,NEXT_REC)
400      CONTINUE
         REC_POINTER(8) = NEXT_REC
      END IF

      RETURN
      END
