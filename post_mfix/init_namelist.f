CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: INIT_NAMELIST                                          C
C  Purpose: initialize the NAMELIST variables                          C
C                                                                      C
C  Author: P. Nicoletti                               Date: 26-NOV-91  C
C  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 27-JAN-92  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: Initialize Phi and Phi_w                                   C
C  Author: M. Syamlal                                 Date: 11-FEB-93  C
C  Revision Number: 2                                                  C
C  Purpose: Add L_scale0, L_scale                                      C
C  Author: W. Sams                                    Date: 04-MAY-94  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: RUN_NAME, DESCRIPTION, UNITS, RUN_TYPE, TIME    C
C                      TSTOP, DT, RES_DT, SPX_DT, OUT_DT, NLOG         C
C                      COORDINATES, IMAX, DX, XLENGTH, JMAX, DY        C
C                      YLENGTH, KMAX, DZ, ZLENGTH, MMAX, D_p, RO_s     C
C                      L_scale0, L_scale, EP_star, MU_g0, MW_AVG       C
C                      IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t  C
C                      IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t  C
C                      IC_EP_g, IC_P_g, IC_ROP_s, IC_T_g, IC_T_s      C
C                       IC_U_g, IC_U_s, IC_V_g, IC_V_s, IC_W_g C
C                      IC_W_s, BC_X_w, BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b  C
C                      BC_Z_t, BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b  C
C                      BC_K_t, BC_EP_g, BC_P_g, BC_RO_g, BC_ROP_g      C
C                      BC_ROP_s, BC_T_g, BC_T_s,  BC_U_g      C
C                      BC_U_s, BC_V_g,BC_V_s, BC_W_g, BC_W_s, BC_TYPE  C
C                      BC_VOLFLOW_g, BC_VOLFLOW_s, BC_MASSFLOW_g       C
C                      BC_MASSFLOW_s, BC_DT_0, BC_Jet_g0, BC_DT_h      C
C                      BC_Jet_gh, BC_DT_l, BC_Jet_gl, NO_I, NO_J, NO_K C
C                      RO_g0, MU_gmax                                  C
C                                                                      C
C  Local variables:  LC, LCM                                           C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE INIT_NAMELIST
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'run.inc'
      INCLUDE 'output.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'ic.inc'
      INCLUDE 'bc.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'is.inc'
      INCLUDE 'toleranc.inc'
      INCLUDE 'scales.inc'
      INCLUDE 'ur_facs.inc'
      INCLUDE 'leqsol.inc'
      INCLUDE 'residual.inc'
      INCLUDE 'rxns.inc'
C
C                      Coefficient of restitution (old symbol)
      DOUBLE PRECISION e
C
      INCLUDE 'namelist.inc'
C
C              loop counters
      INTEGER  LC , LCM, M, N
C
C INITIALIZE THE RUN CONTROL SECTION
C
      RUN_NAME    = UNDEFINED_C
      DESCRIPTION = UNDEFINED_C
      UNITS       = UNDEFINED_C
      RUN_TYPE    = UNDEFINED_C
      TIME        = UNDEFINED
      TSTOP       = UNDEFINED
      DT          = UNDEFINED
      DT_MAX      = 1.
      DT_MIN      = 1.E-6
      DT_FAC      = 0.9
      DETECT_STALL= .TRUE.
      ENERGY_EQ   = .TRUE.
      GRANULAR_ENERGY   = .FALSE.
      DO 10 M = 0, DIMENSION_M
        MOMENTUM_X_EQ(M) = .TRUE.
        MOMENTUM_Y_EQ(M) = .TRUE.
        MOMENTUM_Z_EQ(M) = .TRUE.
        SPECIES_EQ(M)    = .TRUE.
10    CONTINUE
      CALL_USR    = .FALSE.
      MODEL_B     = .FALSE.
      DISCRETIZE(1)  = 0
      DISCRETIZE(2)  = 0
      DISCRETIZE(3)  = 0
      DISCRETIZE(4)  = 0
      DISCRETIZE(5)  = 0
      DISCRETIZE(6)  = 0
      DISCRETIZE(7)  = 0
      DISCRETIZE(8)  = 0
C
C INITIALIZE THE OUTPUT CONTROL SECTION
C
      NLOG        = 25
      FULL_LOG    = .FALSE.
      RES_DT      = UNDEFINED
      DO 50 LC = 1,N_SPX
         SPX_DT(LC)    = UNDEFINED
50    CONTINUE
      OUT_DT      = UNDEFINED
      DO 100 LC = 1,DIMENSION_USR
         USR_DT(LC)      = UNDEFINED
         USR_TYPE(LC)    = UNDEFINED_C
         USR_VAR(LC)     = UNDEFINED_C
         USR_FORMAT(LC)  = UNDEFINED_C
         USR_EXT(LC)     = UNDEFINED_C
100   CONTINUE

      DO 101 LC = 1, 8
        RESID_STRING(LC) = UNDEFINED_C
101   CONTINUE

C
C INITIALIZE THE GEOMETRY AND DISCRETIZATION SECTION
C
      COORDINATES = UNDEFINED_C
      NO_I        = .FALSE.
      NO_J        = .FALSE.
      NO_K        = .FALSE.
      IMAX        = UNDEFINED_I
      JMAX        = UNDEFINED_I
      KMAX        = UNDEFINED_I
      MMAX        = 1
      XMIN        = ZERO
      XLENGTH     = UNDEFINED
      YLENGTH     = UNDEFINED
      ZLENGTH     = UNDEFINED
      DO 120 LC = 1,DIMENSION_I
         DX(LC) = UNDEFINED
120   CONTINUE
      DO 140 LC = 1,DIMENSION_J
         DY(LC) = UNDEFINED
140   CONTINUE
      DO 160 LC = 1,DIMENSION_K
         DZ(LC) = UNDEFINED
160   CONTINUE
      CYCLIC_X    = .FALSE.
      CYCLIC_Y    = .FALSE.
      CYCLIC_Z    = .FALSE.
      CYCLIC_X_PD = .FALSE.
      CYCLIC_Y_PD = .FALSE.
      CYCLIC_Z_PD = .FALSE.
C
C  Constants
C
      DO 170 LC = 1, DIMENSION_C
        C(LC) = UNDEFINED
        C_NAME(LC) = '....................'
170   CONTINUE
      GRAVITY    = UNDEFINED
      C_e        = UNDEFINED
      C_f        = UNDEFINED
      Phi        = UNDEFINED
      Phi_w      = ZERO
      L_scale0   = ZERO
      MU_gmax    = UNDEFINED
      NORM_g     = UNDEFINED
      NORM_s     = UNDEFINED
      TOL_RESID    = 1.0E-3
      TOL_RESID_T  = 1.0E-4
      TOL_RESID_X  = 1.0E-4
      TOL_DIVERGE= 1.0E+4
      V_ex       = ZERO
      P_ref      = ZERO
      P_scale    = ONE
      MAX_NIT    = 500
      LEQ_IT(1)  = 10
      LEQ_IT(2)  = 10
      LEQ_IT(3)  =  5
      LEQ_IT(4)  =  5
      LEQ_IT(5)  =  5
      LEQ_IT(6)  =  5
      LEQ_IT(7)  =  5
      LEQ_IT(8)  =  5
      LEQ_METHOD(1)  = 2
      LEQ_METHOD(2)  = 1
      LEQ_METHOD(3)  = 1
      LEQ_METHOD(4)  = 1
      LEQ_METHOD(5)  = 1
      LEQ_METHOD(6)  = 2
      LEQ_METHOD(7)  = 2
      LEQ_METHOD(8)  = 2
      UR_FAC(1)  = 0.8           !pressure
      UR_FAC(2)  = 0.5           !rho, ep
      UR_FAC(3)  = 0.5           !U
      UR_FAC(4)  = 0.5           !V
      UR_FAC(5)  = 0.5           !W
      UR_FAC(6)  = 0.8           !T
      UR_FAC(7)  = 0.8           !X
      UR_FAC(8)  = 0.5           !Th
C
C INITIALIZE THE GAS PHASE SECTION
C
      RO_g0   = UNDEFINED
      MU_g0   = UNDEFINED
      K_g0    = UNDEFINED
      DIF_g0    = UNDEFINED
      C_pg0   = UNDEFINED
      MW_AVG = UNDEFINED
      NMAX(0) = 1
      DO 200 N = 1, DIMENSION_N_g
        MW_g(N) = UNDEFINED
200   CONTINUE
C
C INITIALIZE THE SOLID PHASE SECTION
C
      MU_s0   = UNDEFINED
      K_s0    = UNDEFINED
      DIF_s0    = UNDEFINED
      C_ps0   = UNDEFINED
      DO 300 LC   = 1,DIMENSION_M
         D_p(LC)  = UNDEFINED
         RO_s(LC) = UNDEFINED
         NMAX(LC) = 1
         CLOSE_PACKED(LC) = .TRUE.
         DO 250 N = 1, DIMENSION_N_s
           MW_s(LC, N) = UNDEFINED
250      CONTINUE
300   CONTINUE
      EP_star     = UNDEFINED
C
C INITIALIZE THE INITIAL CONDITIONS
C
      DO 400 LC = 1,DIMENSION_IC
         IC_X_w(LC)  = UNDEFINED
         IC_X_e(LC)  = UNDEFINED
         IC_Y_s(LC)  = UNDEFINED
         IC_Y_n(LC)  = UNDEFINED
         IC_Z_b(LC)  = UNDEFINED
         IC_Z_t(LC)  = UNDEFINED
         IC_I_w(LC)  = UNDEFINED_I
         IC_I_e(LC)  = UNDEFINED_I
         IC_J_s(LC)  = UNDEFINED_I
         IC_J_n(LC)  = UNDEFINED_I
         IC_K_b(LC)  = UNDEFINED_I
         IC_K_t(LC)  = UNDEFINED_I
         IC_TYPE(LC) = UNDEFINED_C
         IC_EP_g(LC) = UNDEFINED
         IC_P_g(LC)  = UNDEFINED
         IC_P_star(LC)  = UNDEFINED
         IC_L_scale(LC) = UNDEFINED
         IC_T_g(LC)  = UNDEFINED
         IC_GAMA_Rg(LC)  = ZERO
         IC_T_Rg(LC)  = UNDEFINED
         DO 330 N = 1, DIMENSION_N_g
           IC_X_g(LC, N) = UNDEFINED
330      CONTINUE
         IC_U_g(LC)  = UNDEFINED
         IC_V_g(LC)  = UNDEFINED
         IC_W_g(LC)  = UNDEFINED
         DO 350 LCM = 1,DIMENSION_M
            IC_ROP_s(LC,LCM) = UNDEFINED
            IC_U_s(LC,LCM)   = UNDEFINED
            IC_V_s(LC,LCM)   = UNDEFINED
            IC_W_s(LC,LCM)   = UNDEFINED
            IC_T_s(LC,LCM)   = UNDEFINED
            IC_Theta_m(LC,LCM)   = UNDEFINED
            IC_GAMA_Rs(LC,LCM)  = ZERO
            IC_T_Rs(LC,LCM)  = UNDEFINED
            DO 340 N = 1, DIMENSION_N_s
              IC_X_s(LC, LCM, N) = UNDEFINED
340         CONTINUE
350      CONTINUE
400   CONTINUE
C
C INITIALIZE THE BOUNDARY CONDITIONS
C
      DELP_X      = UNDEFINED
      DELP_Y      = UNDEFINED
      DELP_Z      = UNDEFINED
      U_g0        = UNDEFINED
      V_g0        = UNDEFINED
      W_g0        = UNDEFINED
      DO 410 LCM = 1,DIMENSION_M
        U_s0(LCM) = UNDEFINED
        V_s0(LCM) = UNDEFINED
        W_s0(LCM) = UNDEFINED
410   CONTINUE
      DO 500 LC = 1,DIMENSION_BC
         BC_X_w(LC)   = UNDEFINED
         BC_X_e(LC)   = UNDEFINED
         BC_Y_s(LC)   = UNDEFINED
         BC_Y_n(LC)   = UNDEFINED
         BC_Z_b(LC)   = UNDEFINED
         BC_Z_t(LC)   = UNDEFINED
         BC_I_w(LC)   = UNDEFINED_I
         BC_I_e(LC)   = UNDEFINED_I
         BC_J_s(LC)   = UNDEFINED_I
         BC_J_n(LC)   = UNDEFINED_I
         BC_K_b(LC)   = UNDEFINED_I
         BC_K_t(LC)   = UNDEFINED_I
         BC_EP_g(LC)  = UNDEFINED
         BC_P_g(LC)   = UNDEFINED
         BC_ROP_g(LC) = UNDEFINED
         BC_T_g(LC)   = UNDEFINED
         DO 430 N = 1, DIMENSION_N_g
           BC_X_g(LC, N) = UNDEFINED
           BC_hw_X_g(LC, N) = UNDEFINED
           BC_Xw_g(LC, N) = UNDEFINED
           BC_C_X_g(LC, N) = UNDEFINED
430      CONTINUE
         BC_U_g(LC)   = UNDEFINED
         BC_V_g(LC)   = UNDEFINED
         BC_W_g(LC)   = UNDEFINED
         BC_TYPE(LC)  = UNDEFINED_C
         BC_VOLFLOW_g(LC) = UNDEFINED
         BC_MASSFLOW_g(LC) = UNDEFINED
         BC_DT_0(LC) = UNDEFINED
         BC_DT_h(LC) = UNDEFINED
         BC_DT_l(LC) = UNDEFINED
         BC_Jet_g0(LC) = UNDEFINED
         BC_Jet_gh(LC) = UNDEFINED
         BC_Jet_gl(LC) = UNDEFINED

         BC_hw_g(LC) = UNDEFINED
         BC_Uw_g(LC) = UNDEFINED
         BC_Vw_g(LC) = UNDEFINED
         BC_Ww_g(LC) = UNDEFINED

         BC_hw_T_g(LC) = UNDEFINED
         BC_Tw_g(LC) = UNDEFINED
         BC_C_T_g(LC) = UNDEFINED

Cstart kapil&anuj 01/19/98
         BC_JJ_PS(LC) = UNDEFINED_I
Cend   kapil&anuj 01/19/98

         DO 450 LCM = 1,DIMENSION_M
            BC_ROP_s(LC, LCM) = UNDEFINED
            BC_U_s(LC,LCM)  = UNDEFINED
            BC_V_s(LC,LCM)  = UNDEFINED
            BC_W_s(LC,LCM)  = UNDEFINED
            BC_T_s(LC,LCM)  = UNDEFINED
            BC_VOLFLOW_s(LC,LCM) = UNDEFINED
            BC_MASSFLOW_s(LC,LCM) = UNDEFINED

            BC_hw_s(LC,LCM) = UNDEFINED
            BC_Uw_s(LC,LCM) = UNDEFINED
            BC_Vw_s(LC,LCM) = UNDEFINED
            BC_Ww_s(LC,LCM) = UNDEFINED

            BC_hw_T_s(LC, LCM) = UNDEFINED
            BC_Tw_s(LC, LCM) = UNDEFINED
            BC_C_T_s(LC, LCM) = UNDEFINED


            BC_hw_Theta_m(LC, LCM) = UNDEFINED
            BC_Thetaw_m(LC, LCM) = UNDEFINED
            BC_C_Theta_m(LC, LCM) = UNDEFINED
            BC_Theta_m(LC, LCM) = UNDEFINED

            DO 440 N = 1, DIMENSION_N_s
              BC_X_s(LC, LCM, N) = UNDEFINED
              BC_hw_X_s(LC, LCM, N) = UNDEFINED
              BC_Xw_s(LC, LCM, N) = UNDEFINED
              BC_C_X_s(LC, LCM, N) = UNDEFINED
440         CONTINUE
450      CONTINUE
500   CONTINUE
C
C INITIALIZE THE INTERNAL SURFACES
C
      DO 600 LC = 1,DIMENSION_IS
         IS_X_w(LC)   = UNDEFINED
         IS_X_e(LC)   = UNDEFINED
         IS_Y_s(LC)   = UNDEFINED
         IS_Y_n(LC)   = UNDEFINED
         IS_Z_b(LC)   = UNDEFINED
         IS_Z_t(LC)   = UNDEFINED
         IS_I_w(LC)   = UNDEFINED_I
         IS_I_e(LC)   = UNDEFINED_I
         IS_J_s(LC)   = UNDEFINED_I
         IS_J_n(LC)   = UNDEFINED_I
         IS_K_b(LC)   = UNDEFINED_I
         IS_K_t(LC)   = UNDEFINED_I
         IS_PC (LC, 1)   = LARGE_NUMBER
         IS_PC (LC, 2)   = ZERO
         IS_TYPE(LC)  = UNDEFINED_C
         DO 590 M = 1, DIMENSION_M
           IS_VEL_s(LC, M) = ZERO
590      CONTINUE
600   CONTINUE
C
C  Initialize keywords for defining reaction schemes and rates
C

       DO 610 LC = 1, DIMENSION_N_all
         SPECIES_NAME(LC) = UNDEFINED_C
610    CONTINUE

C
C  Initialize user-defined namelist items
C
      CALL USR_INIT_NAMELIST
C
      RETURN
      END
