CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: WRITE_RES1                                             C
C  Purpose: write out the time-dependent restart records               C
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
C  Variables referenced: TIME, NSTEP, EP_g, P_g, P_star, RO_g, ROP_g   C
C                        T_g, T_s, U_g, V_g, W_g, ROP_s, U_s    C
C                        V_s, W_s, IJKMAX2                             C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: LC, N, NEXT_REC                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE WRITE_RES1
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
C             loop counter
      INTEGER LC, N
C
C             pointer to first time-dependent record in restart file
      INTEGER NEXT_REC
C
      READ (UNIT_RES, REC=3)NEXT_REC
      WRITE (UNIT_RES,REC=NEXT_REC) TIME, DT, NSTEP
      NEXT_REC = NEXT_REC + 1
C
      CALL OUT_BIN_512 (UNIT_RES,EP_g,IJKMAX2,NEXT_REC)
      CALL OUT_BIN_512 (UNIT_RES,P_g,IJKMAX2,NEXT_REC)
      CALL OUT_BIN_512 (UNIT_RES,P_star,IJKMAX2,NEXT_REC)
      CALL OUT_BIN_512 (UNIT_RES,RO_g,IJKMAX2,NEXT_REC)
      CALL OUT_BIN_512 (UNIT_RES,ROP_g,IJKMAX2,NEXT_REC)
      CALL OUT_BIN_512 (UNIT_RES,T_g,IJKMAX2,NEXT_REC)
      DO 50 N = 1, NMAX(0)
        CALL OUT_BIN_512 (UNIT_RES,X_g(1,N),IJKMAX2,NEXT_REC)
50    CONTINUE
      CALL OUT_BIN_512 (UNIT_RES,U_g,IJKMAX2,NEXT_REC)
      CALL OUT_BIN_512 (UNIT_RES,V_g,IJKMAX2,NEXT_REC)
      CALL OUT_BIN_512 (UNIT_RES,W_g,IJKMAX2,NEXT_REC)
C
      DO 100 LC = 1,MMAX
         CALL OUT_BIN_512 (UNIT_RES,ROP_s(1,LC),IJKMAX2,NEXT_REC)
         CALL OUT_BIN_512 (UNIT_RES,T_s(1,LC),IJKMAX2,NEXT_REC)
         CALL OUT_BIN_512 (UNIT_RES,U_s(1,LC),IJKMAX2,NEXT_REC)
         CALL OUT_BIN_512 (UNIT_RES,V_s(1,LC),IJKMAX2,NEXT_REC)
         CALL OUT_BIN_512 (UNIT_RES,W_s(1,LC),IJKMAX2,NEXT_REC)
         CALL OUT_BIN_512 (UNIT_RES,THETA_m(1,LC),IJKMAX2,NEXT_REC)
         DO 70 N = 1, NMAX(LC)
           CALL OUT_BIN_512 (UNIT_RES,X_s(1,LC,N),IJKMAX2,NEXT_REC)
70       CONTINUE
100   CONTINUE
C
      CALL FLUSH(UNIT_RES)
C
      RETURN
      END
