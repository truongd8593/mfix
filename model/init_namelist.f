!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INIT_NAMELIST                                          C
!  Purpose: initialize the NAMELIST variables                          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 26-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 27-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Initialize Phi and Phi_w                                   C
!  Author: M. Syamlal                                 Date: 11-FEB-93  C
!  Revision Number: 2                                                  C
!  Purpose: Add L_scale0, L_scale                                      C
!  Author: W. Sams                                    Date: 04-MAY-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: RUN_NAME, DESCRIPTION, UNITS, RUN_TYPE, TIME    C
!                      TSTOP, DT, RES_DT, SPX_DT, OUT_DT, NLOG         C
!                      COORDINATES, IMAX, DX, XLENGTH, JMAX, DY        C
!                      YLENGTH, KMAX, DZ, ZLENGTH, MMAX, D_p, RO_s     C
!                      L_scale0, L_scale, EP_star, MU_g0, MW_AVG       C
!                      IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t  C
!                      IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t  C
!                      IC_EP_g, IC_P_g, IC_ROP_s, IC_T_g, IC_T_s      C
!                       IC_U_g, IC_U_s, IC_V_g, IC_V_s, IC_W_g C
!                      IC_W_s, BC_X_w, BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b  C
!                      BC_Z_t, BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b  C
!                      BC_K_t, BC_EP_g, BC_P_g, BC_RO_g, BC_ROP_g      C
!                      BC_ROP_s, BC_T_g, BC_T_s,  BC_U_g      C
!                      BC_U_s, BC_V_g,BC_V_s, BC_W_g, BC_W_s, BC_TYPE  C
!                      BC_VOLFLOW_g, BC_VOLFLOW_s, BC_MASSFLOW_g       C
!                      BC_MASSFLOW_s, BC_DT_0, BC_Jet_g0, BC_DT_h      C
!                      BC_Jet_gh, BC_DT_l, BC_Jet_gl, NO_I, NO_J, NO_K C
!                      RO_g0, MU_gmax                                  C
!                                                                      C
!  Local variables:  LC, LCM                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE INIT_NAMELIST 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------

      USE param 
      USE param1 
      USE run
      USE output
      USE physprop
      USE geometry
      USE ic
      USE bc
      USE fldvar
      USE constant
      USE indices
      USE is
      USE toleranc 
      USE scales 
      USE ur_facs 
      USE leqsol 
      USE residual
      USE rxns
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!              loop counters
      INTEGER :: LC, LCM, M, N 
!
!                      Coefficient of restitution (old symbol)
      DOUBLE PRECISION :: E 
!-----------------------------------------------
!
!
      INCLUDE 'namelist.inc'

!
!
! INITIALIZE THE RUN CONTROL SECTION
!
      RUN_NAME = UNDEFINED_C 
      DESCRIPTION = UNDEFINED_C 
      UNITS = UNDEFINED_C 
      RUN_TYPE = UNDEFINED_C 
      TIME = UNDEFINED 
      TSTOP = UNDEFINED 
      DT = UNDEFINED 
      DT_MAX = 1. 
      DT_MIN = 1.E-6 
      DT_FAC = 0.9 
      DETECT_STALL = .TRUE. 
      ENERGY_EQ = .TRUE.
      DEF_COR  =  .FALSE. 
      GRANULAR_ENERGY = .FALSE. 
! start peter 7/15
      V_sh=0d0
!
!
! start anuj 04/20
!
      FRICTION = .FALSE. 
      SAVAGE = 1 
! end anuj 04/20
!
 
! start loezos 
!
       SHEAR = .FALSE.
! end loezos  
	
     IF (DIM_M + 1 > 0) THEN 
         MOMENTUM_X_EQ(:DIM_M) = .TRUE. 
         MOMENTUM_Y_EQ(:DIM_M) = .TRUE. 
         MOMENTUM_Z_EQ(:DIM_M) = .TRUE. 
         SPECIES_EQ(:DIM_M) = .TRUE. 
      ENDIF 
      CALL_USR = .FALSE. 
      MODEL_B = .FALSE. 
      DISCRETIZE(1) = 0 
      DISCRETIZE(2) = 0 
      DISCRETIZE(3) = 0 
      DISCRETIZE(4) = 0 
      DISCRETIZE(5) = 0 
      DISCRETIZE(6) = 0 
      DISCRETIZE(7) = 0 
      DISCRETIZE(8) = 0 
      DISCRETIZE(9) = 0 
!
! INITIALIZE THE OUTPUT CONTROL SECTION
!
      NLOG = 25 
      FULL_LOG = .FALSE. 
      RES_DT = UNDEFINED 
      SPX_DT(:N_SPX) = UNDEFINED 
      LC = N_SPX + 1 
      OUT_DT = UNDEFINED 
      DO LC = 1, DIMENSION_USR 
         USR_DT(LC) = UNDEFINED 
         USR_TYPE(LC) = UNDEFINED_C 
         USR_VAR(LC) = UNDEFINED_C 
         USR_FORMAT(LC) = UNDEFINED_C 
         USR_EXT(LC) = UNDEFINED_C 
      END DO 
      DO LC = 1, 8 
         RESID_STRING(LC) = UNDEFINED_C 
      END DO 
      COORDINATES = UNDEFINED_C 
      NO_I = .FALSE. 
      NO_J = .FALSE. 
      NO_K = .FALSE. 
      IMAX = UNDEFINED_I 
      JMAX = UNDEFINED_I 
      KMAX = UNDEFINED_I 
      MMAX = 1 
      XMIN = ZERO 
      XLENGTH = UNDEFINED 
      YLENGTH = UNDEFINED 
      ZLENGTH = UNDEFINED 
      DX(:DIM_I) = UNDEFINED 
      DY(:DIM_J) = UNDEFINED 
      DZ(:DIM_K) = UNDEFINED 
      CYCLIC_X = .FALSE. 
      CYCLIC_Y = .FALSE. 
      CYCLIC_Z = .FALSE. 
      CYCLIC_X_PD = .FALSE. 
      CYCLIC_Y_PD = .FALSE. 
      CYCLIC_Z_PD = .FALSE. 
!
!  Constants
!
      DO LC = 1, DIMENSION_C 
         C(LC) = UNDEFINED 
         C_NAME(LC) = '....................' 
      END DO 
      GRAVITY = UNDEFINED 
      C_E = UNDEFINED 
      C_F = UNDEFINED 
      PHI = UNDEFINED 
      PHIP = 0.6D0 
      E_W = 1D0 
      EPS_F_MIN = 0.5D0 
      EPS_MAX = 0.65D0 
      PHI_W = ZERO 
      L_SCALE0 = ZERO 
      MU_GMAX = UNDEFINED 
      NORM_G = UNDEFINED 
      NORM_S = UNDEFINED 
      TOL_RESID = 1.0E-3 
      TOL_RESID_T = 1.0E-4 
      TOL_RESID_X = 1.0E-4 
      TOL_DIVERGE = 1.0E+4 
      V_EX = ZERO 
      P_REF = ZERO 
      P_SCALE = ONE 
      MAX_NIT = 500 
      LEQ_IT(1) = 20 
      LEQ_IT(2) = 20 
      LEQ_IT(3) = 5 
      LEQ_IT(4) = 5 
      LEQ_IT(5) = 5 
      LEQ_IT(6) = 15 
      LEQ_IT(7) = 15
      LEQ_IT(8) = 15
      LEQ_IT(9) = 15
      LEQ_METHOD(1) = 2 
      LEQ_METHOD(2) = 1 
      LEQ_METHOD(3) = 1 
      LEQ_METHOD(4) = 1 
      LEQ_METHOD(5) = 1 
      LEQ_METHOD(6) = 2 
      LEQ_METHOD(7) = 2 
      LEQ_METHOD(8) = 2 
      LEQ_METHOD(9) = 1 
      LEQ_SWEEP(1) = 'II' 
      LEQ_SWEEP(2) = 'II'  
      LEQ_SWEEP(3) = 'II' 
      LEQ_SWEEP(4) = 'II'  
      LEQ_SWEEP(5) = 'II'  
      LEQ_SWEEP(6) = 'II'  
      LEQ_SWEEP(7) = 'II' 
      LEQ_SWEEP(8) = 'II'  
      LEQ_SWEEP(9) = 'II'  
      LEQ_TOL(1) = 1.0D-4 
      LEQ_TOL(2) = 1.0D-4  
      LEQ_TOL(3) = 1.0D-4 
      LEQ_TOL(4) = 1.0D-4  
      LEQ_TOL(5) = 1.0D-4  
      LEQ_TOL(6) = 1.0D-4  
      LEQ_TOL(7) = 1.0D-4 
      LEQ_TOL(8) = 1.0D-4  
      LEQ_TOL(9) = 1.0D-4  
      UR_FAC(1) = 0.8                            !pressure 
      UR_FAC(2) = 0.5                            !rho, ep 
      UR_FAC(3) = 0.5                            !U 
      UR_FAC(4) = 0.5                            !V 
      UR_FAC(5) = 0.5                            !W 
      UR_FAC(6) = 0.8                            !T 
      UR_FAC(7) = 0.75                           !X 
      UR_FAC(8) = 0.5                            !Th 
      UR_FAC(8) = 0.9                            !Scalar
!
! INITIALIZE THE GAS PHASE SECTION
!
      RO_G0 = UNDEFINED 
      MU_G0 = UNDEFINED 
      K_G0 = UNDEFINED 
      DIF_G0 = UNDEFINED 
      C_PG0 = UNDEFINED 
      MW_AVG = UNDEFINED 
      NMAX(0) = UNDEFINED_I 
      MW_G(:DIM_N_G) = UNDEFINED 
      MU_S0 = UNDEFINED 
      K_S0 = UNDEFINED 
      DIF_S0 = UNDEFINED 
      C_PS0 = UNDEFINED 
      D_P(:DIM_M) = UNDEFINED 
      RO_S(:DIM_M) = UNDEFINED 
      NMAX(1:DIM_M) = UNDEFINED_I 
      CLOSE_PACKED(:DIM_M) = .TRUE. 
!      MW_S(1,1+:DIM_N_S+) = UNDEFINED 
      MW_S = UNDEFINED
      EP_STAR = UNDEFINED 
!
! INITIALIZE THE INITIAL CONDITIONS
!
      DO LC = 1, DIMENSION_IC 
         IC_X_W(LC) = UNDEFINED 
         IC_X_E(LC) = UNDEFINED 
         IC_Y_S(LC) = UNDEFINED 
         IC_Y_N(LC) = UNDEFINED 
         IC_Z_B(LC) = UNDEFINED 
         IC_Z_T(LC) = UNDEFINED 
         IC_I_W(LC) = UNDEFINED_I 
         IC_I_E(LC) = UNDEFINED_I 
         IC_J_S(LC) = UNDEFINED_I 
         IC_J_N(LC) = UNDEFINED_I 
         IC_K_B(LC) = UNDEFINED_I 
         IC_K_T(LC) = UNDEFINED_I 
         IC_TYPE(LC) = UNDEFINED_C 
         IC_EP_G(LC) = UNDEFINED 
         IC_P_G(LC) = UNDEFINED 
         IC_P_STAR(LC) = UNDEFINED 
         IC_L_SCALE(LC) = UNDEFINED 
         IC_T_G(LC) = UNDEFINED 
         IC_GAMA_RG(LC) = ZERO 
         IC_T_RG(LC) = UNDEFINED 
         IC_X_G(LC,:DIM_N_G) = UNDEFINED 
         IC_U_G(LC) = UNDEFINED 
         IC_V_G(LC) = UNDEFINED 
         IC_W_G(LC) = UNDEFINED 
         IC_ROP_S(LC,:DIM_M) = UNDEFINED 
         IC_U_S(LC,:DIM_M) = UNDEFINED 
         IC_V_S(LC,:DIM_M) = UNDEFINED 
         IC_W_S(LC,:DIM_M) = UNDEFINED 
         IC_T_S(LC,:DIM_M) = UNDEFINED 
         IC_THETA_M(LC,:DIM_M) = UNDEFINED 
         IC_GAMA_RS(LC,:DIM_M) = ZERO 
         IC_T_RS(LC,:DIM_M) = UNDEFINED 
!         IC_X_S(LC,1,1+:DIM_N_S+) = UNDEFINED 
      END DO 
      IC_X_S = UNDEFINED
      DELP_X = UNDEFINED 
      DELP_Y = UNDEFINED 
      DELP_Z = UNDEFINED 
      U_G0 = UNDEFINED 
      V_G0 = UNDEFINED 
      W_G0 = UNDEFINED 
      U_S0(:DIM_M) = UNDEFINED 
      V_S0(:DIM_M) = UNDEFINED 
      W_S0(:DIM_M) = UNDEFINED 
      DO LC = 1, DIMENSION_BC 
         BC_X_W(LC) = UNDEFINED 
         BC_X_E(LC) = UNDEFINED 
         BC_Y_S(LC) = UNDEFINED 
         BC_Y_N(LC) = UNDEFINED 
         BC_Z_B(LC) = UNDEFINED 
         BC_Z_T(LC) = UNDEFINED 
         BC_I_W(LC) = UNDEFINED_I 
         BC_I_E(LC) = UNDEFINED_I 
         BC_J_S(LC) = UNDEFINED_I 
         BC_J_N(LC) = UNDEFINED_I 
         BC_K_B(LC) = UNDEFINED_I 
         BC_K_T(LC) = UNDEFINED_I 
         BC_EP_G(LC) = UNDEFINED 
         BC_P_G(LC) = UNDEFINED 
         BC_ROP_G(LC) = UNDEFINED 
         BC_T_G(LC) = UNDEFINED 
         BC_X_G(LC,:DIM_N_G) = UNDEFINED 
         BC_HW_X_G(LC,:DIM_N_G) = UNDEFINED 
         BC_XW_G(LC,:DIM_N_G) = UNDEFINED 
         BC_C_X_G(LC,:DIM_N_G) = UNDEFINED 
         BC_U_G(LC) = UNDEFINED 
         BC_V_G(LC) = UNDEFINED 
         BC_W_G(LC) = UNDEFINED 
         BC_TYPE(LC) = UNDEFINED_C 
         BC_VOLFLOW_G(LC) = UNDEFINED 
         BC_MASSFLOW_G(LC) = UNDEFINED 
         BC_DT_0(LC) = UNDEFINED 
         BC_DT_H(LC) = UNDEFINED 
         BC_DT_L(LC) = UNDEFINED 
         BC_JET_G0(LC) = UNDEFINED 
         BC_JET_GH(LC) = UNDEFINED 
         BC_JET_GL(LC) = UNDEFINED 
!
         BC_HW_G(LC) = UNDEFINED 
         BC_UW_G(LC) = UNDEFINED 
         BC_VW_G(LC) = UNDEFINED 
         BC_WW_G(LC) = UNDEFINED 
!
         BC_HW_T_G(LC) = UNDEFINED 
         BC_TW_G(LC) = UNDEFINED 
         BC_C_T_G(LC) = UNDEFINED 
!
!start kapil&anuj 01/19/98
         BC_JJ_PS(LC) = UNDEFINED_I 
!end   kapil&anuj 01/19/98
!
         BC_ROP_S(LC,:DIM_M) = UNDEFINED 
         BC_U_S(LC,:DIM_M) = UNDEFINED 
         BC_V_S(LC,:DIM_M) = UNDEFINED 
         BC_W_S(LC,:DIM_M) = UNDEFINED 
         BC_T_S(LC,:DIM_M) = UNDEFINED 
         BC_VOLFLOW_S(LC,:DIM_M) = UNDEFINED 
         BC_MASSFLOW_S(LC,:DIM_M) = UNDEFINED 
!
         BC_HW_S(LC,:DIM_M) = UNDEFINED 
         BC_UW_S(LC,:DIM_M) = UNDEFINED 
         BC_VW_S(LC,:DIM_M) = UNDEFINED 
         BC_WW_S(LC,:DIM_M) = UNDEFINED 
!
         BC_HW_T_S(LC,:DIM_M) = UNDEFINED 
         BC_TW_S(LC,:DIM_M) = UNDEFINED 
         BC_C_T_S(LC,:DIM_M) = UNDEFINED 
!
!
         BC_HW_THETA_M(LC,:DIM_M) = UNDEFINED 
         BC_THETAW_M(LC,:DIM_M) = UNDEFINED 
         BC_C_THETA_M(LC,:DIM_M) = UNDEFINED 
!         BC_THETA_M(LC,:DIM_M) = UNDEFINED 
!
!         BC_X_S(LC,1,1+:DIM_N_S+) = UNDEFINED 
!         BC_HW_X_S(LC,1,1+:DIM_N_S+) = UNDEFINED 
!         BC_XW_S(LC,1,1+:DIM_N_S+) = UNDEFINED 
!         BC_C_X_S(LC,1,1+:DIM_N_S+) = UNDEFINED 
      END DO 
!   
      BC_THETA_M = UNDEFINED 
      BC_X_S = UNDEFINED 
      BC_HW_X_S = UNDEFINED 
      BC_XW_S = UNDEFINED 
      BC_C_X_S = UNDEFINED 
!
      DO LC = 1, DIMENSION_IS 
         IS_X_W(LC) = UNDEFINED 
         IS_X_E(LC) = UNDEFINED 
         IS_Y_S(LC) = UNDEFINED 
         IS_Y_N(LC) = UNDEFINED 
         IS_Z_B(LC) = UNDEFINED 
         IS_Z_T(LC) = UNDEFINED 
         IS_I_W(LC) = UNDEFINED_I 
         IS_I_E(LC) = UNDEFINED_I 
         IS_J_S(LC) = UNDEFINED_I 
         IS_J_N(LC) = UNDEFINED_I 
         IS_K_B(LC) = UNDEFINED_I 
         IS_K_T(LC) = UNDEFINED_I 
         IS_PC(LC,1) = LARGE_NUMBER 
         IS_PC(LC,2) = ZERO 
         IS_TYPE(LC) = UNDEFINED_C 
         IS_VEL_S(LC,:DIM_M) = ZERO 
      END DO 
      DO LC = 1, DIM_N_ALL
         SPECIES_NAME(LC) = UNDEFINED_C 
      END DO 

      NODESI = UNDEFINED_I
      NODESJ = UNDEFINED_I
      NODESK = UNDEFINED_I

      CALL USR_INIT_NAMELIST 
!
      RETURN  
      END SUBROUTINE INIT_NAMELIST 
