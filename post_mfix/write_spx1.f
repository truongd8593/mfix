CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: WRITE_SPX1                                             C
C  Purpose: write out the time-dependent restart records (REAL)        C
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
C  Variables referenced: TIME, NSTEP, EP_g, P_g, P_star, U_g           C
C                        V_g, W_g, U_s, V_s, W_s, ROP_s, T_g, T_s     C
C                        IJKMAX2, MMAX                           C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables:  LC, N, NEXT_REC, NUM_REC                          C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE WRITE_SPX1(L)
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
      INCLUDE 'output.inc'
C
C passed arguments
C
C             flag whether to write a particular SPx file 
      INTEGER L
C
C local variables
C
C             loop counters
      INTEGER LC, N
C
C             Pointer to the next record
      INTEGER NEXT_REC
C
C              Number of records written each time step
      INTEGER  NUM_REC
C
C ".SP1" FILE         EP_g    [ ROP_g, RO_g  must be calculated ...
C                                        not written out ]
C
      IF (L .EQ. 1) THEN
         READ  (UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (UNIT_SPX+L,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL OUT_BIN_R(UNIT_SPX+L,EP_g,IJKMAX2,NEXT_REC)
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         CALL FLUSH(UNIT_SPX+L)
C
C ".SP2" FILE         P_g , P_star
C
      ELSEIF (L .EQ. 2) THEN
         READ  (UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (UNIT_SPX+L,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL OUT_BIN_R(UNIT_SPX+L,P_g,IJKMAX2,NEXT_REC)
         CALL OUT_BIN_R(UNIT_SPX+L,P_star,IJKMAX2,NEXT_REC)
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         CALL FLUSH(UNIT_SPX+L)
C
C ".SP3" FILE         U_g , V_g , W_g
C
      ELSEIF (L .EQ. 3) THEN
         READ  (UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (UNIT_SPX+L,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL OUT_BIN_R(UNIT_SPX+L,U_g,IJKMAX2,NEXT_REC)
         CALL OUT_BIN_R(UNIT_SPX+L,V_g,IJKMAX2,NEXT_REC)
         CALL OUT_BIN_R(UNIT_SPX+L,W_g,IJKMAX2,NEXT_REC)
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         CALL FLUSH(UNIT_SPX+L)
C
C ".SP4" FILE         U_s , V_s , W_s
C
      ELSEIF (L .EQ. 4) THEN
         READ  (UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (UNIT_SPX+L,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 100 LC = 1,MMAX
            CALL OUT_BIN_R(UNIT_SPX+L,U_s(1,LC),IJKMAX2,NEXT_REC)
            CALL OUT_BIN_R(UNIT_SPX+L,V_s(1,LC),IJKMAX2,NEXT_REC)
            CALL OUT_BIN_R(UNIT_SPX+L,W_s(1,LC),IJKMAX2,NEXT_REC)
100      CONTINUE
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         CALL FLUSH(UNIT_SPX+L)
C
C ".SP5" FILE         ROP_s
C
      ELSEIF (L .EQ. 5) THEN
         READ  (UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (UNIT_SPX+L,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 200 LC = 1,MMAX
            CALL OUT_BIN_R(UNIT_SPX+L,ROP_s(1,LC),IJKMAX2,NEXT_REC)
200      CONTINUE
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(UNIT_SPX+L, REC=3) NEXT_REC, NUM_REC
         CALL FLUSH(UNIT_SPX+L)
C
C ".SP6" FILE         T_g  , T_s
C
      ELSEIF (L .EQ. 6) THEN
         READ  (UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (UNIT_SPX+L,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL OUT_BIN_R(UNIT_SPX+L,T_g,IJKMAX2,NEXT_REC)
         DO 220 LC = 1,MMAX
           CALL OUT_BIN_R(UNIT_SPX+L,T_s(1,LC),IJKMAX2,NEXT_REC)
220      CONTINUE
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         CALL FLUSH(UNIT_SPX+L)
C
C ".SP7" FILE         X_g, X_s
C
      ELSEIF (L .EQ. 7) THEN
         READ  (UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (UNIT_SPX+L,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 250 N = 1, NMAX(0)
           CALL OUT_BIN_R(UNIT_SPX+L,X_g(1,N),IJKMAX2,NEXT_REC)
250      CONTINUE
         DO 300 LC = 1,MMAX
           DO 270 N = 1, NMAX(LC)
             CALL OUT_BIN_R(UNIT_SPX+L,X_s(1,LC, N),IJKMAX2,NEXT_REC)
270        CONTINUE
300      CONTINUE
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(UNIT_SPX+L, REC=3) NEXT_REC, NUM_REC
         CALL FLUSH(UNIT_SPX+L)
C
C ".SP8" FILE         THETA_m
C
      ELSEIF (L .EQ. 8) THEN
         READ  (UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (UNIT_SPX+L,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         DO 320 LC = 1,MMAX
           CALL OUT_BIN_R(UNIT_SPX+L,THETA_m(1,LC),IJKMAX2,NEXT_REC)
320      CONTINUE
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(UNIT_SPX+L,REC=3) NEXT_REC, NUM_REC
         CALL FLUSH(UNIT_SPX+L)
C

      END IF
      RETURN
      END
