CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: WRITE_RES0                                             C
C  Purpose: write out the initial restart records (namelist data)      C
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
C  Variables referenced: RUN_NAME  ,  ID_MONTH  ,  ID_DAY , ID_YEAR    C
C                        ID_HOUR, ID_MINUTE, ID_SECOND, IMAX, JMAX     C
C                        KMAX, IMAX1, JMAX1, KMAX1, IMAX2, JMAX2,KMAX2 C
C                        IJMAX2, IJKMAX2, MMAX, DT, XLENGTH, YLENGTH   C
C                        ZLENGTH, DX, DY, DZ, RUN_NAME, DESCRIPTION    C
C                        UNITS, RUN_TYPE, CORDINATES, D_p, RO_s,       C
C                        EP_star, MU_g0, MW_AVG, IC_X_w, IC_X_e, IC_Y_sC
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
C                        BC_MASSFLOW_s, BC_TYPE, FLAG, RO_g0           C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: LC, L, N, NEXT_RECA, VERSION                       C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE WRITE_RES0
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'run.inc'
      INCLUDE 'ic.inc'
      INCLUDE 'is.inc'
      INCLUDE 'bc.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'output.inc'
      INCLUDE 'scales.inc'
      INCLUDE 'ur_facs.inc'
      INCLUDE 'leqsol.inc'
      INCLUDE 'toleranc.inc'
C
C                loop counters
      INTEGER    LC, L, N
C
C                Pointer to the next record
      INTEGER    NEXT_RECA
C
C                file version id
      CHARACTER  VERSION*512
C
      NEXT_RECA = 5
C
C     Add new data entries at the end of the file and identify version no.
C------------------------------------------------------------------------
      VERSION = 'RES = 01.2'
C------------------------------------------------------------------------
C
      WRITE (UNIT_RES,REC=1) VERSION
      WRITE (UNIT_RES,REC=2) RUN_NAME,ID_MONTH,ID_DAY,ID_YEAR,ID_HOUR,
     &                 ID_MINUTE,ID_SECOND
      WRITE (UNIT_RES,REC=3) NEXT_RECA
      WRITE (UNIT_RES,REC=4) IMIN1, JMIN1, KMIN1,
     &                 IMAX, JMAX, KMAX, IMAX1 , JMAX1 , KMAX1 ,
     &                 IMAX2 , JMAX2 , KMAX2 , IJMAX2 , IJKMAX2 ,
     &                 MMAX , DIMENSION_IC , DIMENSION_BC , DIMENSION_C,
     &                 DIMENSION_IS ,
     &                 DT , XMIN, XLENGTH , YLENGTH , ZLENGTH,
     &                 C_e, C_f, Phi, Phi_w
      CALL OUT_BIN_512 (UNIT_RES,C,DIMENSION_C,NEXT_RECA)
      next_reca = 1 + next_reca    ! work around for -O3 compiler bug
      next_reca = next_reca - 1
c      write(*,*)next_reca 
      DO 50 LC = 1, DIMENSION_C
c        write(*,*)next_reca,LC
        WRITE (UNIT_RES,REC=NEXT_RECA)C_NAME(LC)
        NEXT_RECA = NEXT_RECA + 1
50    CONTINUE

      WRITE (UNIT_RES,REC=NEXT_RECA)(NMAX(L),L=0,MMAX)
      NEXT_RECA = NEXT_RECA + 1

c      CALL OUT_BIN_512I (UNIT_RES, NMAX, MMAX+1, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,DX,IMAX2,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,DY,JMAX2,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,DZ,KMAX2,NEXT_RECA)
      WRITE (UNIT_RES,REC=NEXT_RECA) RUN_NAME , DESCRIPTION , UNITS ,
     &                         RUN_TYPE ,COORDINATES
      NEXT_RECA = NEXT_RECA + 1
      WRITE (UNIT_RES,REC=NEXT_RECA) (D_p(L),L=1,MMAX), 
     &                 (RO_s(L),L=1,MMAX),
     &                 EP_star , RO_g0, MU_g0 , MW_AVG
      NEXT_RECA = NEXT_RECA + 1
      CALL OUT_BIN_512 (UNIT_RES,MW_g,NMAX(0),NEXT_RECA)
      DO 60 LC = 1, MMAX
        WRITE (UNIT_RES,REC=NEXT_RECA)(MW_s(LC, N), N=1,NMAX(LC))
        NEXT_RECA = NEXT_RECA + 1
60    CONTINUE
      CALL OUT_BIN_512 (UNIT_RES,IC_X_w,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_X_e,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_Y_s,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_Y_n,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_Z_b,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_Z_t,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IC_I_w,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IC_I_e,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IC_J_s,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IC_J_n,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IC_K_b,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IC_K_t,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_EP_g,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_P_g,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_T_g,DIMENSION_IC,NEXT_RECA)
      DO 70 N = 1, NMAX(0)
        CALL OUT_BIN_512 (UNIT_RES,IC_X_g(1,N),DIMENSION_IC,NEXT_RECA)
70    CONTINUE
      CALL OUT_BIN_512 (UNIT_RES,IC_U_g,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_V_g,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_W_g,DIMENSION_IC,NEXT_RECA)
      DO 100 LC = 1,MMAX
         CALL OUT_BIN_512 (UNIT_RES,IC_ROP_s(1,LC),DIMENSION_IC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,IC_U_s(1,LC),DIMENSION_IC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,IC_V_s(1,LC),DIMENSION_IC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,IC_W_s(1,LC),DIMENSION_IC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,IC_T_s(1,LC),DIMENSION_IC,
     &                     NEXT_RECA)
        DO 80 N = 1, NMAX(LC)
          CALL OUT_BIN_512 (UNIT_RES,IC_X_s(1,LC,N),DIMENSION_IC,
     &                     NEXT_RECA)
80      CONTINUE
100   CONTINUE
C
      CALL OUT_BIN_512 (UNIT_RES,BC_X_w,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_X_e,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_Y_s,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_Y_n,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_Z_b,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_Z_t,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,BC_I_w,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,BC_I_e,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,BC_J_s,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,BC_J_n,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,BC_K_b,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,BC_K_t,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_EP_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_P_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_T_g,DIMENSION_BC,NEXT_RECA)
      DO 170 N = 1, NMAX(0)
        CALL OUT_BIN_512 (UNIT_RES,BC_X_g(1,N),DIMENSION_BC,NEXT_RECA)
170   CONTINUE
      CALL OUT_BIN_512 (UNIT_RES,BC_U_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_V_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_W_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_RO_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_ROP_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_VOLFLOW_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_MASSFLOW_g,DIMENSION_BC,NEXT_RECA)
      DO 200 LC = 1,MMAX
         CALL OUT_BIN_512 (UNIT_RES,BC_ROP_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,BC_U_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,BC_V_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,BC_W_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,BC_T_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
         DO 180 N = 1, NMAX(LC)
           CALL OUT_BIN_512 (UNIT_RES,BC_X_s(1,LC,N),DIMENSION_BC,
     &                     NEXT_RECA)
180      CONTINUE
         CALL OUT_BIN_512 (UNIT_RES,BC_VOLFLOW_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,BC_MASSFLOW_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
200   CONTINUE
C
      DO 300 LC = 1,DIMENSION_BC
         WRITE (UNIT_RES,REC=NEXT_RECA) BC_TYPE(LC)
         NEXT_RECA = NEXT_RECA + 1
300   CONTINUE
C
      CALL OUT_BIN_512I(UNIT_RES,FLAG,IJKMAX2,NEXT_RECA)
C
      CALL OUT_BIN_512 (UNIT_RES,IS_X_w,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IS_X_e,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IS_Y_s,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IS_Y_n,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IS_Z_b,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IS_Z_t,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IS_I_w,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IS_I_e,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IS_J_s,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IS_J_n,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IS_K_b,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512I(UNIT_RES,IS_K_t,DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IS_PC(1,1),DIMENSION_IS,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IS_PC(1,2),DIMENSION_IS,NEXT_RECA)
      DO 340 LC = 1, MMAX
        CALL OUT_BIN_512 (UNIT_RES,IS_VEL_s(1,LC),DIMENSION_IS,
     &                     NEXT_RECA)
340   CONTINUE
C
      DO 350 LC = 1,DIMENSION_IS
         WRITE (UNIT_RES,REC=NEXT_RECA) IS_TYPE(LC)
         NEXT_RECA = NEXT_RECA + 1
350   CONTINUE
C
C  The statements above constitute version 01.07
C  Add new write statements below and identify data file version
C
C     Version 01.08
      WRITE (UNIT_RES,REC=NEXT_RECA)CYCLIC_X, CYCLIC_Y, CYCLIC_Z,
     &  CYCLIC_X_PD, CYCLIC_Y_PD, CYCLIC_Z_PD, DELP_X, DELP_Y, DELP_Z,
     &  U_g0, U_s0, V_g0, V_s0, W_g0, W_s0
      NEXT_RECA = NEXT_RECA + 1
C
C     Version 01.09
      WRITE (UNIT_RES,REC=NEXT_RECA) TIME, TSTOP, ENERGY_EQ,
     &  RES_DT, OUT_DT, NLOG, L_scale0, NO_I, NO_J, NO_K, CALL_USR
      NEXT_RECA = NEXT_RECA + 1
      DO 360 LC = 1,N_SPX
         WRITE (UNIT_RES,REC=NEXT_RECA) SPX_DT(LC)
         NEXT_RECA = NEXT_RECA + 1
360   CONTINUE
      DO 361 LC = 0, MMAX
         WRITE (UNIT_RES,REC=NEXT_RECA) SPECIES_EQ(LC)
         NEXT_RECA = NEXT_RECA + 1
361   CONTINUE
      CALL OUT_BIN_512 (UNIT_RES,USR_DT,DIMENSION_USR,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,USR_X_w,DIMENSION_USR,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,USR_X_e,DIMENSION_USR,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,USR_Y_s,DIMENSION_USR,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,USR_Y_n,DIMENSION_USR,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,USR_Z_b,DIMENSION_USR,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,USR_Z_t,DIMENSION_USR,NEXT_RECA)
      DO 365 LC = 1, DIMENSION_USR
         WRITE (UNIT_RES,REC=NEXT_RECA)
     &     USR_FORMAT(LC), USR_EXT(LC), USR_TYPE(LC), USR_VAR(LC)
         NEXT_RECA = NEXT_RECA + 1
365   CONTINUE
      CALL OUT_BIN_512 (UNIT_RES,IC_P_star,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_L_scale,DIMENSION_IC,NEXT_RECA)
      DO 366 LC = 1, DIMENSION_IC
         WRITE (UNIT_RES,REC=NEXT_RECA)IC_TYPE (LC)
         NEXT_RECA = NEXT_RECA + 1
366   CONTINUE
      CALL OUT_BIN_512 (UNIT_RES,BC_DT_0,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_Jet_g0,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_DT_h,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_Jet_gh,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_DT_l,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_Jet_gl,DIMENSION_BC,NEXT_RECA)
C
C     Version 01.10
      WRITE (UNIT_RES,REC=NEXT_RECA) MU_gmax
      NEXT_RECA = NEXT_RECA + 1
C
C     Version 01.11
      WRITE (UNIT_RES,REC=NEXT_RECA) V_ex, MODEL_B
      NEXT_RECA = NEXT_RECA + 1
C
C     Version 01.12
      WRITE (UNIT_RES,REC=NEXT_RECA) P_ref, P_scale, UR_FAC,
     &                               TOL_RESID, DT_MAX, DT_MIN, DT_FAC,
     &                               CLOSE_PACKED, GRAVITY, MU_s0
      NEXT_RECA = NEXT_RECA + 1
      WRITE (UNIT_RES,REC=NEXT_RECA)LEQ_IT, LEQ_METHOD
      NEXT_RECA = NEXT_RECA + 1
      CALL OUT_BIN_512 (UNIT_RES,BC_hw_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_Uw_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_Vw_g,DIMENSION_BC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,BC_Ww_g,DIMENSION_BC,NEXT_RECA)
      DO 400 LC = 1,MMAX
         CALL OUT_BIN_512 (UNIT_RES,BC_hw_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,BC_Uw_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,BC_Vw_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
         CALL OUT_BIN_512 (UNIT_RES,BC_Ww_s(1,LC),DIMENSION_BC,
     &                     NEXT_RECA)
400   CONTINUE
C
C     Version 01.13
      WRITE (UNIT_RES,REC=NEXT_RECA) MOMENTUM_X_EQ, MOMENTUM_Y_EQ,
     &  MOMENTUM_Z_EQ, TOL_DIVERGE, DISCRETIZE, FULL_LOG
      NEXT_RECA = NEXT_RECA + 1
C
C     Version 01.14
      WRITE (UNIT_RES,REC=NEXT_RECA) DETECT_STALL
      NEXT_RECA = NEXT_RECA + 1
C
C     Version 01.15
      WRITE (UNIT_RES,REC=NEXT_RECA) K_g0, K_s0, C_pg0, C_ps0,
     &               TOL_RESID_T, TOL_RESID_X
      NEXT_RECA = NEXT_RECA + 1
      CALL OUT_BIN_512 (UNIT_RES,IC_GAMA_Rg,DIMENSION_IC,NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES,IC_T_Rg,DIMENSION_IC,NEXT_RECA)
      DO 420 LC = 1, MMAX
        CALL OUT_BIN_512 (UNIT_RES,IC_GAMA_Rs(1,LC),DIMENSION_IC,
     &                     NEXT_RECA)
        CALL OUT_BIN_512 (UNIT_RES,IC_T_Rs(1,LC),DIMENSION_IC,
     &                     NEXT_RECA)
420   CONTINUE
C
C     Version 01.2
      WRITE (UNIT_RES,REC=NEXT_RECA) NORM_g, NORM_s
      NEXT_RECA = NEXT_RECA + 1
C
C  Add new write statements above this line.  Remember to update NEXT_RECA.
C  Remember to change the version number near begining of this subroutine.
C  Also modify READ_RES0.  The routines such as OUT_BIN_512 etc. writes out
C  arrays dimensioned ARRAY(DIM).  So arrays dimensioned ARRAY(DIM1:DIM2)
C  should be passed as ARRAY(DIM1) and array length as DIM2-DIM1+1.
C---------------------------------------------------------------------------
      WRITE (UNIT_RES,REC=3) NEXT_RECA
      CALL FLUSH(UNIT_RES)
      RETURN
      END
