CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: READ_RES1                                              C
C  Purpose: read in the time-dependent restart records                 C
C                                                                      C
C  Author: P. Nicoletti                               Date: 03-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: IJKMAX2, MMAX, DT                             C
C  Variables modified: TIME, NSTEP, EP_g, P_g, P_star, RO_g            C
C                      ROP_g, T_g, T_s,  U_g, V_g, W_g, ROP_s    C
C                      U_s, V_s, W_s                                   C
C                                                                      C
C  Local variables: TIME_READ, LC, NEXT_REC                            C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE READ_RES1
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
C
C             loop counter
      INTEGER LC 
C 
C                      Local species index
      INTEGER          N
C
C             pointer to the next record
      INTEGER NEXT_REC
C
C                file version id
      CHARACTER  VERSION*512
C
C                version number
      REAL       VERSION_NUMBER
C
C                      Dummy array
      DOUBLE PRECISION Tmp(DIMENSION_3)
C
      READ (UNIT_RES,REC=1) VERSION
      READ(VERSION(6:512),*)VERSION_NUMBER
C
      READ (UNIT_RES, REC=3) NEXT_REC
      IF(VERSION_NUMBER .GE. 1.12)THEN
        READ (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP
      ELSE
        READ (UNIT_RES, REC=NEXT_REC) TIME, NSTEP
      ENDIF

      NEXT_REC = NEXT_REC + 1
C
      CALL IN_BIN_512 (UNIT_RES,EP_g,IJKMAX2,NEXT_REC)
      CALL IN_BIN_512 (UNIT_RES,P_g,IJKMAX2,NEXT_REC)
      CALL IN_BIN_512 (UNIT_RES,P_star,IJKMAX2,NEXT_REC)
      CALL IN_BIN_512 (UNIT_RES,RO_g,IJKMAX2,NEXT_REC)
      CALL IN_BIN_512 (UNIT_RES,ROP_g,IJKMAX2,NEXT_REC)
      CALL IN_BIN_512 (UNIT_RES,T_g,IJKMAX2,NEXT_REC)
      IF(VERSION_NUMBER .LT. 1.15)THEN
        CALL IN_BIN_512 (UNIT_RES,T_s(1,1),IJKMAX2,NEXT_REC)
        IF(MMAX .GE. 2) THEN
          CALL IN_BIN_512 (UNIT_RES,T_s(1,2),IJKMAX2,NEXT_REC)
        ELSE
          CALL IN_BIN_512 (UNIT_RES,Tmp,IJKMAX2,NEXT_REC)
        ENDIF
      ENDIF
      IF(VERSION_NUMBER .GE. 1.05)THEN
        DO 50 N = 1, NMAX(0)
          CALL IN_BIN_512 (UNIT_RES,X_g(1,N),IJKMAX2,NEXT_REC)
50      CONTINUE
      ENDIF
      CALL IN_BIN_512 (UNIT_RES,U_g,IJKMAX2,NEXT_REC)
      CALL IN_BIN_512 (UNIT_RES,V_g,IJKMAX2,NEXT_REC)
      CALL IN_BIN_512 (UNIT_RES,W_g,IJKMAX2,NEXT_REC)
C
      DO 100 LC = 1,MMAX
         CALL IN_BIN_512 (UNIT_RES,ROP_s(1,LC),IJKMAX2,NEXT_REC)
         IF(VERSION_NUMBER .GE. 1.15)THEN
           CALL IN_BIN_512 (UNIT_RES,T_s(1,LC),IJKMAX2,NEXT_REC)
         ENDIF
         CALL IN_BIN_512 (UNIT_RES,U_s(1,LC),IJKMAX2,NEXT_REC)
         CALL IN_BIN_512 (UNIT_RES,V_s(1,LC),IJKMAX2,NEXT_REC)
         CALL IN_BIN_512 (UNIT_RES,W_s(1,LC),IJKMAX2,NEXT_REC)
         IF(VERSION_NUMBER .GE. 1.2)THEN
           CALL IN_BIN_512 (UNIT_RES,THETA_m(1,LC),IJKMAX2,NEXT_REC)
         ENDIF
         IF(VERSION_NUMBER .GE. 1.05)THEN
           DO 70 N = 1, NMAX(LC)
             CALL IN_BIN_512 (UNIT_RES,X_s(1,LC,N),IJKMAX2,NEXT_REC)
70         CONTINUE
         ENDIF
100   CONTINUE
C
      RETURN
      END
