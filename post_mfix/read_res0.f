CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: READ_RES0                                              C
C  Purpose: read the initial restart records (namelist data)           C
C                                                                      C
C  Author: P. Nicoletti                               Date: 13-DEC-91  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables modified : RUN_NAME  ,  ID_MONTH  ,  ID_DAY , ID_YEAR     C
C                        ID_HOUR, ID_MINUTE, ID_SECOND, IMAX, JMAX     C
C                        KMAX, IMAX1, JMAX1, KMAX1, IMAX2, JMAX2,KMAX2 C
C                        IJMAX2, IJKMAX2, MMAX, DT, XLENGTH, YLENGTH   C
C                        ZLENGTH, DX, DY, DZ, RUN_NAME, DESCRIPTION    C
C                        UNITS, RUN_TYPE, CORDINATES, D_p, RO_s,       C
C                        EP_star, MU_g, MW_AVG, IC_X_w, IC_X_e, IC_Y_s C
C                        IC_Y_n, IC_Z_b, IC_Z_t, IC_I_w, IC_I_e        C
C                        IC_J_s, IC_J_n, IC_K_b, IC_K_t, IC_EP_g       C
C                        IC_P_g, IC_T_g, IC_T_s,  IC_U_g      C
C                        IC_V_g, IC_W_g, IC_ROP_s, IC_U_s, IC_V_s      C
C                        IC_W_s, BC_X_w, BC_X_e, BC_Y_s, BC_Y_n        C
C                        BC_Z_b, BC_Z_t, BC_I_w, BC_I_e, BC_J_s        C
C                        BC_K_b, BC_K_t, BC_EP_g, BC_P_g, BC_T_g       C
C                        BC_T_s,  BC_U_g, BC_V_g, BC_W_g      C
C                        BC_RO_g, BC_ROP_g, BC_VOLFLOW_g,BC_MASSFLOW_g C
C                        BC_ROP_s, BC_U_s, BC_V_s, BC_VOLFLOW_s        C
C                        BC_MASSFLOW_s, BC_TYPE, FLAG                  C
C  Variables referenced: None                                          C
C                                                                      C
C  Local variables: LC, L, NEXT_RECA, VERSION                          C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE READ_RES0
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'run.inc'
      INCLUDE 'ic.inc'
      INCLUDE 'bc.inc'
      INCLUDE 'is.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'output.inc'
      INCLUDE 'scales.inc'
      INCLUDE 'ur_facs.inc'
      INCLUDE 'toleranc.inc'
      INCLUDE 'leqsol.inc'
C
C                loop counters
      INTEGER    LC, L , N, M
C
C                Pointer to the next record
      INTEGER    NEXT_RECA
C
C                file version id
      CHARACTER  VERSION*512
C
C                version number
      REAL       VERSION_NUMBER
C
C                      Temporary arrays
      DOUBLE PRECISION IC_Tmp(DIMENSION_IC), BC_Tmp(DIMENSION_BC)
C
      INTEGER    DIM_IC , DIM_BC , DIM_C , DIM_IS
C
C  1) Check to ensure that this subroutine was updated.  
C  2) Initialize missing constants from earlier versions.
C  3) Add new read statements at the end of the file. 
C
      READ (UNIT_RES,REC=1) VERSION
      READ(VERSION(6:512),*)VERSION_NUMBER
      IF(VERSION_NUMBER .GT. 1.2) THEN
        WRITE(*,*)' Update Subroutine read_res0'
        CALL SLUMBER
        STOP
      ENDIF
C
C  Initialize required constants missing from earlier versions
      P_ref   = ZERO
      P_scale = ONE
      DIM_IC = 5
      DIM_BC = 5
      DIM_C  = 5
      DIM_IS = 5
      c_e = 1.0
      c_f = 0.0
      phi = 0.0
      phi_w = 0.0
C
C
      READ (UNIT_RES,REC=2) RUN_NAME,ID_MONTH,ID_DAY,ID_YEAR,ID_HOUR,
     &                 ID_MINUTE,ID_SECOND
      READ (UNIT_RES,REC=3) NEXT_RECA
      IF (VERSION.EQ.'RES = 01.00') THEN
         READ (UNIT_RES,REC=4) IMIN1, JMIN1, KMIN1,
     &                 IMAX, JMAX, KMAX, IMAX1 , JMAX1 , KMAX1 ,
     &                 IMAX2 , JMAX2 , KMAX2 , IJMAX2 , IJKMAX2 ,
     &                 MMAX , DT , XLENGTH ,
     &                 YLENGTH , ZLENGTH
      ELSE IF (VERSION.EQ.'RES = 01.01' .OR. 
     &                        VERSION.EQ.'RES = 01.02') THEN
         READ (UNIT_RES,REC=4) IMIN1, JMIN1, KMIN1,
     &                 IMAX, JMAX, KMAX, IMAX1 , JMAX1 , KMAX1 ,
     &                 IMAX2 , JMAX2 , KMAX2 , IJMAX2 , IJKMAX2 ,
     &                 MMAX , DIM_IC, DIM_BC , DT , XLENGTH ,
     &                 YLENGTH , ZLENGTH
      ELSE IF (VERSION.EQ.'RES = 01.03') THEN
         READ (UNIT_RES,REC=4) IMIN1, JMIN1, KMIN1,
     &                 IMAX, JMAX, KMAX, IMAX1 , JMAX1 , KMAX1 ,
     &                 IMAX2 , JMAX2 , KMAX2 , IJMAX2 , IJKMAX2 ,
     &                 MMAX , DIM_IC, DIM_BC , DT , XMIN , XLENGTH ,
     &                 YLENGTH , ZLENGTH
      ELSE IF (VERSION.EQ.'RES = 01.04') THEN
         READ (UNIT_RES,REC=4) IMIN1, JMIN1, KMIN1,
     &                 IMAX, JMAX, KMAX, IMAX1 , JMAX1 , KMAX1 ,
     &                 IMAX2 , JMAX2 , KMAX2 , IJMAX2 , IJKMAX2 ,
     &                 MMAX , DIM_IC, DIM_BC , DIM_C ,
     &                 DT , XMIN , XLENGTH , YLENGTH , ZLENGTH
      ELSE IF (VERSION.EQ.'RES = 01.05') THEN
         READ (UNIT_RES,REC=4) IMIN1, JMIN1, KMIN1,
     &                 IMAX, JMAX, KMAX, IMAX1 , JMAX1 , KMAX1 ,
     &                 IMAX2 , JMAX2 , KMAX2 , IJMAX2 , IJKMAX2 ,
     &                 MMAX , DIM_IC, DIM_BC , DIM_C , DIM_IS,
     &                 DT , XMIN , XLENGTH , YLENGTH , ZLENGTH
      ELSE
         READ (UNIT_RES,REC=4) IMIN1, JMIN1, KMIN1,
     &                 IMAX, JMAX, KMAX, IMAX1 , JMAX1 , KMAX1 ,
     &                 IMAX2 , JMAX2 , KMAX2 , IJMAX2 , IJKMAX2 ,
     &                 MMAX , DIM_IC, DIM_BC , DIM_C , DIM_IS,
     &                 DT , XMIN , XLENGTH , YLENGTH , ZLENGTH,
     &                 C_e, C_f, Phi, Phi_w
      END IF

C
C CHECK DIMENSIONS
C
      IF (IMAX2.GT.DIMENSION_I) GOTO 900
      IF (JMAX2.GT.DIMENSION_J) GOTO 900
      IF (KMAX2.GT.DIMENSION_K) GOTO 900
      IF (IJKMAX2.GT.DIMENSION_3) GOTO 900
      IF (MMAX.GT.DIMENSION_M)  GOTO 900
      IF (DIM_IC.GT.DIMENSION_IC) GOTO 900
      IF (DIM_BC.GT.DIMENSION_BC) GOTO 900
      IF (DIM_C.GT.DIMENSION_C) GOTO 900
      IF (DIM_IS.GT.DIMENSION_IS) GOTO 900
      DO 10 M = 0, MMAX
        NMAX(M) = 1
10    CONTINUE
C
      NEXT_RECA = 5
C
      IF (VERSION_NUMBER .GE. 1.04) THEN
         CALL IN_BIN_512 (UNIT_RES,C,DIM_C,NEXT_RECA)
         next_reca = 1 + next_reca    ! work around for -O3 compiler bug
         next_reca = next_reca - 1
         DO LC  = 1,DIM_C
            READ(UNIT_RES,REC=NEXT_RECA) C_NAME(LC)
            NEXT_RECA = NEXT_RECA + 1
         END DO
        IF (VERSION_NUMBER .LT. 1.12) THEN
          CALL IN_BIN_512I (UNIT_RES,NMAX,MMAX+1,NEXT_RECA)
        ELSE
          READ (UNIT_RES,REC=NEXT_RECA)(NMAX(L),L=0,MMAX)
          NEXT_RECA = NEXT_RECA + 1
        ENDIF
      END IF
C
      IF(NMAX(0) .GT. DIMENSION_N_g) GOTO 900
      DO 20 M = 1, MMAX
        IF(NMAX(M) .GT. DIMENSION_N_s) GOTO 900
20    CONTINUE
C
      CALL IN_BIN_512 (UNIT_RES,DX,IMAX2,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,DY,JMAX2,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,DZ,KMAX2,NEXT_RECA)
C
      READ (UNIT_RES,REC=NEXT_RECA) RUN_NAME , DESCRIPTION , UNITS ,
     &                         RUN_TYPE ,COORDINATES
      NEXT_RECA = NEXT_RECA + 1
      IF (VERSION .EQ. 'RES = 01.00' .OR. 
     &    VERSION .EQ. 'RES = 01.01') THEN
         READ (UNIT_RES,REC=NEXT_RECA) (D_p(L),L=1,MMAX), 
     &                 (RO_s(L),L=1,MMAX),
     &                 EP_star , MU_g0 , MW_AVG
      ELSE IF (VERSION.EQ.'RES = 01.02') THEN
         READ (UNIT_RES,REC=NEXT_RECA) (D_p(L),L=1,MMAX), 
     &                 (RO_s(L),L=1,MMAX),
     &                 EP_star , RO_g0 , MU_g0 , MW_AVG
      ELSE IF (VERSION.EQ.'RES = 01.03') THEN
         READ (UNIT_RES,REC=NEXT_RECA) (D_p(L),L=1,MMAX), 
     &                 (RO_s(L),L=1,MMAX),
     &                 EP_star , RO_g0 , MU_g0 , MW_AVG
      ELSEIF (VERSION_NUMBER .GE. 1.04) THEN
         READ (UNIT_RES,REC=NEXT_RECA) (D_p(L),L=1,MMAX), 
     &                 (RO_s(L),L=1,MMAX),
     &                 EP_star , RO_g0 , MU_g0 , MW_AVG
      END IF
      NEXT_RECA = NEXT_RECA + 1
C
      IF (VERSION_NUMBER .GE. 1.04) THEN
         CALL IN_BIN_512(UNIT_RES,MW_g,NMAX(0),NEXT_RECA)
         DO LC = 1,MMAX
            READ(UNIT_RES,REC=NEXT_RECA)(MW_s(LC,N),N=1,NMAX(LC))
            NEXT_RECA = NEXT_RECA + 1
         END DO
      END IF
C
      CALL IN_BIN_512 (UNIT_RES,IC_X_w,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,IC_X_e,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,IC_Y_s,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,IC_Y_n,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,IC_Z_b,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,IC_Z_t,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,IC_I_w,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,IC_I_e,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,IC_J_s,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,IC_J_n,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,IC_K_b,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,IC_K_t,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,IC_EP_g,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,IC_P_g,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,IC_T_g,DIM_IC,NEXT_RECA)
      IF(VERSION_NUMBER .LT. 1.15) THEN
        CALL IN_BIN_512 (UNIT_RES,IC_T_s(1, 1),DIM_IC,NEXT_RECA)
        IF(MMAX .GE. 2)THEN
          CALL IN_BIN_512 (UNIT_RES,IC_T_s(1,2),DIM_IC,NEXT_RECA)
        ELSE
          CALL IN_BIN_512 (UNIT_RES,IC_Tmp,DIM_IC,NEXT_RECA)
        ENDIF
      ENDIF
C
      IF (VERSION_NUMBER .GE. 1.04) THEN
         DO N = 1,NMAX(0)
            CALL IN_BIN_512(UNIT_RES,IC_X_g(1,N),DIM_IC,NEXT_RECA)
         END DO
      END IF
C
      CALL IN_BIN_512 (UNIT_RES,IC_U_g,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,IC_V_g,DIM_IC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,IC_W_g,DIM_IC,NEXT_RECA)
      DO 100 LC = 1,MMAX
         CALL IN_BIN_512 (UNIT_RES,IC_ROP_s(1,LC),DIM_IC,
     &                     NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES,IC_U_s(1,LC),DIM_IC,
     &                     NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES,IC_V_s(1,LC),DIM_IC,
     &                     NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES,IC_W_s(1,LC),DIM_IC,
     &                     NEXT_RECA)
         IF(VERSION_NUMBER .GE. 1.15) THEN
           CALL IN_BIN_512 (UNIT_RES,IC_T_s(1,LC),DIM_IC,
     &                     NEXT_RECA)
         ENDIF
C
      IF (VERSION_NUMBER .GE. 1.04) THEN
         DO N = 1,NMAX(LC)
            CALL IN_BIN_512 (UNIT_RES,IC_X_s(1,LC,N),DIM_IC,
     &                     NEXT_RECA)
         END DO
      END IF
C
100   CONTINUE
C
      CALL IN_BIN_512 (UNIT_RES,BC_X_w,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_X_e,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_Y_s,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_Y_n,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_Z_b,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_Z_t,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,BC_I_w,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,BC_I_e,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,BC_J_s,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,BC_J_n,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,BC_K_b,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512I(UNIT_RES,BC_K_t,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_EP_g,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_P_g,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_T_g,DIM_BC,NEXT_RECA)
      IF(VERSION_NUMBER .LT. 1.15) THEN
        CALL IN_BIN_512 (UNIT_RES,BC_T_s(1,1),DIM_BC,NEXT_RECA)
        IF(MMAX .GE. 2)THEN
          CALL IN_BIN_512 (UNIT_RES,BC_T_s(1,2),DIM_BC,NEXT_RECA)
        ELSE
          CALL IN_BIN_512 (UNIT_RES,BC_Tmp,DIM_BC,NEXT_RECA)
        ENDIF
      ENDIF
C
      IF (VERSION_NUMBER .GE. 1.04) THEN
         DO N = 1,NMAX(0)
            CALL IN_BIN_512(UNIT_RES,BC_X_g(1,N),DIM_BC,NEXT_RECA)
         END DO
      END IF
C
      CALL IN_BIN_512 (UNIT_RES,BC_U_g,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_V_g,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_W_g,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_RO_g,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_ROP_g,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_VOLFLOW_g,DIM_BC,NEXT_RECA)
      CALL IN_BIN_512 (UNIT_RES,BC_MASSFLOW_g,DIM_BC,NEXT_RECA)
      DO 200 LC = 1,MMAX
         CALL IN_BIN_512 (UNIT_RES,BC_ROP_s(1,LC),DIM_BC,
     &                     NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES,BC_U_s(1,LC),DIM_BC,
     &                     NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES,BC_V_s(1,LC),DIM_BC,
     &                     NEXT_RECA)
C
C Note : previous versions did not write out BC_W_s
C
        IF (VERSION_NUMBER .GE. 1.04) THEN
            CALL IN_BIN_512 (UNIT_RES,BC_W_s(1,LC),DIM_BC,
     &                     NEXT_RECA)

            IF(VERSION_NUMBER .GE. 1.15) THEN
              CALL IN_BIN_512 (UNIT_RES,BC_T_s(1,LC),DIM_BC,
     &                     NEXT_RECA)
            ENDIF

            DO N = 1,NMAX(LC)
              CALL IN_BIN_512 (UNIT_RES,BC_X_s(1,LC,N),DIM_BC,
     &                     NEXT_RECA)
            END DO
        END IF
C
         CALL IN_BIN_512 (UNIT_RES,BC_VOLFLOW_s(1,LC),DIM_BC,
     &                     NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES,BC_MASSFLOW_s(1,LC),DIM_BC,
     &                     NEXT_RECA)
200   CONTINUE
C
      IF (VERSION.EQ.'RES = 01.00') THEN
         L = 10
      ELSE
         L = DIM_BC
      END IF
      DO 300 LC = 1,L
         READ (UNIT_RES,REC=NEXT_RECA) BC_TYPE(LC)
         NEXT_RECA = NEXT_RECA + 1
300   CONTINUE
C
      CALL IN_BIN_512I(UNIT_RES,FLAG,IJKMAX2,NEXT_RECA)
C
      IF (VERSION_NUMBER .GE. 1.04) THEN
         CALL IN_BIN_512(UNIT_RES,IS_X_w,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512(UNIT_RES,IS_X_e,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512(UNIT_RES,IS_Y_s,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512(UNIT_RES,IS_Y_n,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512(UNIT_RES,IS_Z_b,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512(UNIT_RES,IS_Z_t,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512I(UNIT_RES,IS_I_w,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512I(UNIT_RES,IS_I_e,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512I(UNIT_RES,IS_J_s,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512I(UNIT_RES,IS_J_n,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512I(UNIT_RES,IS_K_b,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512I(UNIT_RES,IS_K_t,DIM_IS,NEXT_RECA)
         CALL IN_BIN_512(UNIT_RES,IS_PC(1,1),DIM_IS,NEXT_RECA)
         CALL IN_BIN_512(UNIT_RES,IS_PC(1,2),DIM_IS,NEXT_RECA)
         IF(VERSION_NUMBER .GE. 1.07)THEN
           DO 340 LC = 1, MMAX
             CALL IN_BIN_512 (UNIT_RES,IS_VEL_s(1,LC),DIM_IS,
     &                     NEXT_RECA)
340        CONTINUE
         ENDIF
         DO LC = 1,DIM_IS
            READ(UNIT_RES,REC=NEXT_RECA) IS_TYPE(LC)
            NEXT_RECA = NEXT_RECA + 1
         END DO
      END IF
C
C  Additions from new versions of .RES file
C
      IF(VERSION_NUMBER .GE. 1.08)THEN
        READ (UNIT_RES,REC=NEXT_RECA)CYCLIC_X, CYCLIC_Y, CYCLIC_Z,
     &    CYCLIC_X_PD, CYCLIC_Y_PD, CYCLIC_Z_PD, DELP_X, DELP_Y, DELP_Z,
     &    U_g0, U_s0, V_g0, V_s0, W_g0, W_s0
        NEXT_RECA = NEXT_RECA + 1
      ENDIF
C
      IF (VERSION_NUMBER .GE. 1.09) THEN
        READ (UNIT_RES,REC=NEXT_RECA) TIME, TSTOP, ENERGY_EQ,
     &  RES_DT, OUT_DT, NLOG, L_scale0, NO_I, NO_J, NO_K, CALL_USR
        NEXT_RECA = NEXT_RECA + 1
        DO 360 LC = 1,N_SPX
          READ (UNIT_RES,REC=NEXT_RECA) SPX_DT(LC)
          NEXT_RECA = NEXT_RECA + 1
360     CONTINUE
        DO 361 LC = 0, MMAX
          READ (UNIT_RES,REC=NEXT_RECA) SPECIES_EQ(LC)
          NEXT_RECA = NEXT_RECA + 1
361     CONTINUE
        CALL IN_BIN_512 (UNIT_RES,USR_DT,DIMENSION_USR,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,USR_X_w,DIMENSION_USR,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,USR_X_e,DIMENSION_USR,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,USR_Y_s,DIMENSION_USR,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,USR_Y_n,DIMENSION_USR,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,USR_Z_b,DIMENSION_USR,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,USR_Z_t,DIMENSION_USR,NEXT_RECA)
        DO 365 LC = 1, DIMENSION_USR
          READ (UNIT_RES,REC=NEXT_RECA)
     &     USR_FORMAT(LC), USR_EXT(LC), USR_TYPE(LC), USR_VAR(LC)
          NEXT_RECA = NEXT_RECA + 1
365     CONTINUE
        CALL IN_BIN_512 (UNIT_RES,IC_P_star,DIM_IC,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,IC_L_scale,DIM_IC,NEXT_RECA)
        DO 366 LC = 1, DIM_IC
          READ (UNIT_RES,REC=NEXT_RECA)IC_TYPE (LC)
          NEXT_RECA = NEXT_RECA + 1
366     CONTINUE
        CALL IN_BIN_512 (UNIT_RES,BC_DT_0,DIM_BC,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,BC_Jet_g0,DIM_BC,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,BC_DT_h,DIM_BC,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,BC_Jet_gh,DIM_BC,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,BC_DT_l,DIM_BC,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,BC_Jet_gl,DIM_BC,NEXT_RECA)
      ENDIF

      IF (VERSION_NUMBER .GE. 1.10) THEN
        READ (UNIT_RES,REC=NEXT_RECA) MU_gmax
        NEXT_RECA = NEXT_RECA + 1
      ENDIF

      IF (VERSION_NUMBER .GE. 1.11) THEN
        READ (UNIT_RES,REC=NEXT_RECA) V_ex, MODEL_B
        NEXT_RECA = NEXT_RECA + 1
      ENDIF

      IF (VERSION_NUMBER .GE. 1.12) THEN
        READ (UNIT_RES,REC=NEXT_RECA) P_ref, P_scale, UR_FAC,
     &                                TOL_RESID, DT_MAX, DT_MIN, DT_FAC,
     &                                CLOSE_PACKED, GRAVITY, MU_s0
        NEXT_RECA = NEXT_RECA + 1
        READ (UNIT_RES,REC=NEXT_RECA)LEQ_IT, LEQ_METHOD
        NEXT_RECA = NEXT_RECA + 1
        CALL IN_BIN_512 (UNIT_RES,BC_hw_g,DIM_BC,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,BC_Uw_g,DIM_BC,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,BC_Vw_g,DIM_BC,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,BC_Ww_g,DIM_BC,NEXT_RECA)
        DO 400 LC = 1,MMAX
          CALL IN_BIN_512 (UNIT_RES,BC_hw_s(1,LC),DIM_BC,
     &                     NEXT_RECA)
          CALL IN_BIN_512 (UNIT_RES,BC_Uw_s(1,LC),DIM_BC,
     &                     NEXT_RECA)
          CALL IN_BIN_512 (UNIT_RES,BC_Vw_s(1,LC),DIM_BC,
     &                     NEXT_RECA)
          CALL IN_BIN_512 (UNIT_RES,BC_Ww_s(1,LC),DIM_BC,
     &                     NEXT_RECA)
400     CONTINUE
      ENDIF

      DO 410 LC = 0, MMAX
        MOMENTUM_X_EQ(LC) = .TRUE.
        MOMENTUM_Y_EQ(LC) = .TRUE.
        MOMENTUM_Z_EQ(LC) = .TRUE.
410   CONTINUE
      TOL_DIVERGE = 1.E+4
      IF (VERSION_NUMBER .GE. 1.13) THEN
        READ (UNIT_RES,REC=NEXT_RECA) MOMENTUM_X_EQ, MOMENTUM_Y_EQ,
     &  MOMENTUM_Z_EQ, TOL_DIVERGE, DISCRETIZE, FULL_LOG
        NEXT_RECA = NEXT_RECA + 1
      ENDIF

      IF (VERSION_NUMBER .GE. 1.14) THEN
        READ (UNIT_RES,REC=NEXT_RECA) DETECT_STALL
        NEXT_RECA = NEXT_RECA + 1
      ENDIF

      IF (VERSION_NUMBER .GE. 1.15) THEN
        READ (UNIT_RES,REC=NEXT_RECA) K_g0, K_s0, C_pg0, C_ps0,
     &               TOL_RESID_T, TOL_RESID_X
        NEXT_RECA = NEXT_RECA + 1
        CALL IN_BIN_512 (UNIT_RES,IC_GAMA_Rg,DIM_IC,NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES,IC_T_Rg,DIM_IC,NEXT_RECA)
        DO 420 LC = 1, MMAX
          CALL IN_BIN_512 (UNIT_RES,IC_GAMA_Rs(1,LC),DIM_IC,
     &                     NEXT_RECA)
          CALL IN_BIN_512 (UNIT_RES,IC_T_Rs(1,LC),DIM_IC,
     &                     NEXT_RECA)
420     CONTINUE
      ENDIF

      IF (VERSION_NUMBER .GE. 1.2) THEN
        READ (UNIT_RES,REC=NEXT_RECA) NORM_g, NORM_s 
        NEXT_RECA = NEXT_RECA + 1
      ENDIF
C
C  Add new read statements above this line.  Remember to update NEXT_RECA.
C  Remember to update the version number check near begining of this subroutine.
C------------------------------------------------------------------------------
C
      READ (UNIT_RES,REC=3) NEXT_RECA
C
C  Since the value of UNDEFINED was changed ...
C
      IF(RO_g0 .GE. 1E30) RO_g0 = UNDEFINED
      IF(MU_g0 .GE. 1E30) MU_g0 = UNDEFINED
      IF(MW_AVG .GE. 1E30) MW_AVG = UNDEFINED
      IF(C_E .GE. 1E30) C_E = UNDEFINED
C
      RETURN
C
C HERE IF DIMENSION ERROR
C
900   CONTINUE
      WRITE (*,*) ' '
      WRITE (*,*) ' **************************************'
      WRITE (*,*) ' From: READ_RES0'
      WRITE (*,*) ' DIMENSION ERROR ---'
      WRITE (*,*) ' '
      WRITE (*,*) ' DIMENSION_I  = ' , DIMENSION_I,
     &            ' IMAX2        = ' , IMAX2
      WRITE (*,*) ' DIMENSION_J  = ' , DIMENSION_J,
     &            ' JMAX2        = ' , JMAX2
      WRITE (*,*) ' DIMENSION_K  = ' , DIMENSION_K,
     &            ' KMAX2        = ' , KMAX2
      WRITE (*,*) ' DIMENSION_3  = ' , DIMENSION_3,
     &            ' IJKMAX2      = ' , IJKMAX2
      WRITE (*,*) ' DIMENSION_M  = ' , DIMENSION_M,
     &            ' MMAX         = ' , MMAX
      WRITE (*,*) ' DIMENSION_IC = ' , DIMENSION_IC,
     &            ' DIM_IC       = ' , DIM_IC
      WRITE (*,*) ' DIMENSION_BC = ' , DIMENSION_BC,
     &            ' DIM_BC       = ' , DIM_BC
      WRITE (*,*) ' DIMENSION_IS = ' , DIMENSION_IS,
     &            ' DIM_IS       = ' , DIM_IS
      WRITE (*,*) ' DIMENSION_C  = ' , DIMENSION_C,
     &            ' DIM_C        = ' , DIM_C
      WRITE (*,*) ' DIMENSION_N_g= ' , DIMENSION_N_g,
     &            ' NMAX(0)      = ' , NMAX(0)
      WRITE (*,*) ' DIMENSION_N_s= ' , DIMENSION_N_s
      DO 50 M = 1, MMAX
        WRITE (*,'(A, I2, A, I4)') ' NMAX(',M,') = ',NMAX(M)
50    CONTINUE
      WRITE (*,*) ' '
C
      STOP
      END
