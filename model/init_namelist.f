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
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Add L_scale0, L_scale                                      C
!  Author: W. Sams                                    Date: 04-MAY-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: To call DES_Init_Namelist                                  C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: To call CARTESIAN_GRID_INIT_NAMELIST                       C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE INIT_NAMELIST 

!-----------------------------------------------
! Modules
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
      USE scalars
      USE compar
      USE parallel
      USE cdist
      USE stiff_chem          !added by Sebastien Dartevelle, LANL, May 2013
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop counters
      INTEGER :: LC, LCM, M, N 
! Coefficient of restitution (old symbol)
      DOUBLE PRECISION :: E       
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------      
      INCLUDE 'namelist.inc'
!-----------------------------------------------

! INITIALIZE THE RUN CONTROL SECTION
      RUN_NAME = UNDEFINED_C 
      DESCRIPTION = UNDEFINED_C 
      UNITS = UNDEFINED_C 
      RUN_TYPE = UNDEFINED_C 
      TIME = UNDEFINED 
      TSTOP = UNDEFINED 
! AEOLUS: STOP Trigger mechanism to terminate MFIX normally before 
!         batch queue terminates
      CHK_BATCHQ_END = .FALSE. 
      BATCH_WALLCLOCK = 9000.0    ! set to 2.5 hrs for jaguarcnl w/ nproc<=512
      TERM_BUFFER = 180.0         ! set to 3 minutes prior to end of job

      DT = UNDEFINED 
      DT_MAX = ONE 
      DT_MIN = 1.D-6 
      DT_FAC = 0.9D0
      DETECT_STALL = .TRUE. 
      ENERGY_EQ = .TRUE.
      DEF_COR  =  .FALSE. 
      C_FAC = UNDEFINED
      FPFOI = .FALSE.
      GRANULAR_ENERGY = .FALSE.
!     do not use revised JJ BC 
      BC_JJ_M = .false.
      PHIP_OUT_JJ=.false.
      PHIP_OUT_ITER=0

 !
 !SUBGRID model stuffs, Sebastien Dartevelel, LANL, May 2013
      SUBGRID_Igci = .FALSE.          !Igci & Sundar/Princeton's subgrid model
      SUBGRID_Milioli = .FALSE.       !Milioli & Sundar/Princeton's subgrid model
      SUBGRID_Wall = .FALSE.          !by default do not calculate the effects of the wall upon the filtered model
      filter_size_ratio = 2.0D0       !filter_size_ratio, which must be defined when SUBGRID = .T. is by default equal to TWO
 !
 !
      K_Epsilon = .FALSE.
      Added_Mass = .FALSE.
      M_AM = UNDEFINED_I
      SIMONIN = .FALSE.
      AHMADI = .FALSE.
      JENKINS = .FALSE.
      YU_STANDISH = .FALSE.
      FEDORS_LANDEL = .FALSE.
      AUTOMATIC_RESTART = .FALSE.
      AUTO_RESTART = .FALSE.
      ITER_RESTART = 1

! peter: 7/15
      V_sh=0d0

! jeg: 4/01/2005
      KT_TYPE = UNDEFINED_C
      RDF_TYPE = 'LEBOWITZ'

! anuj: 04/20
      FRICTION = .FALSE. 
      SAVAGE = 1 

! sof: 02/16/2005
      SCHAEFFER = .TRUE. 
 
! sp: 02/08/2006
      BLENDING_STRESS = .FALSE. 
      TANH_BLEND      = .TRUE.
      SIGM_BLEND      = .FALSE.

! sp: 06/15/2007
      DEBUG_RESID     = .TRUE.

! DISTIO distributed IO
      bDist_IO            = .false.
      bStart_with_one_RES = .false.

! loezos:
      SHEAR = .FALSE.

      DRAG_TYPE = 'SYAM_OBRIEN'
      drag_c1 = 0.8d0
      drag_d1 = 2.65d0
      LAM_HYS = UNDEFINED

! AE: 041601 Set the default to 1st order accurate time implementation
      CN_ON = .FALSE.
              
      IF (DIM_M + 1 > 0) THEN 
         MOMENTUM_X_EQ(:DIM_M) = .TRUE. 
         MOMENTUM_Y_EQ(:DIM_M) = .TRUE. 
         MOMENTUM_Z_EQ(:DIM_M) = .TRUE. 
         SPECIES_EQ(:DIM_M) = .TRUE. 
      ENDIF 
      CALL_USR = .FALSE. 
      MODEL_B = .FALSE. 
      DISCRETIZE(:) = 0 
      Chi_scheme = .FALSE.
      
      NScalar = 0
      Phase4Scalar(:) = UNDEFINED_I

      nRR = 0
      Call_DQMOM = .FALSE.      


! INITIALIZE THE OUTPUT CONTROL SECTION
      report_mass_balance_dt = UNDEFINED
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


! Constants
      DO LC = 1, DIMENSION_C 
         C(LC) = UNDEFINED 
         C_NAME(LC) = '....................' 
      ENDDO 
      GRAVITY = UNDEFINED 
      C_E = UNDEFINED 
      C_F = UNDEFINED 
      PHI = UNDEFINED 
      PHIP = 0.6D0 
! specularity coefficient as r->0
!
!!    IF ( (SUBGRID_Igci .OR. SUBGRID_Milioli) .AND. SUBGRID_Wall) PHIP = 0.0D0         !it must be ZEROED with any Filter models as we only take FREE_Slip wall for the time being
!                                                                                       !Sebastien Dartevelle, LANL, May 2013
!
      k4phi = undefined	
      phip0 = undefined            
      E_W = 1.D0 
      PHI_W = UNDEFINED 
      EPS_F_MIN = 0.5D0 
      L_SCALE0 = ZERO 
      MU_GMAX = UNDEFINED 
      V_EX = ZERO 
      P_REF = ZERO 
      P_SCALE = ONE 

! gera: 08/15/03
      SEGREGATION_SLOPE_COEFFICIENT=0.D0
      EP_S_MAX(:DIM_M) = UNDEFINED

! coefficient of restitution
      r_p(:DIM_M, :DIM_M) = UNDEFINED

! rong
      AGGREGATION_EFF=0.D0
      BREAKAGE_EFF=0.D0

! numerics     
      NORM_G = UNDEFINED 
      NORM_S = UNDEFINED 
      TOL_RESID = 1.0D-3 
      TOL_RESID_T = 1.0D-4 
      TOL_RESID_X = 1.0D-4 
      TOL_RESID_Scalar = 1.0D-4
      TOL_RESID_K_Epsilon = 1.0D-4
      TOL_RESID_Th = 1.0D-4
      TOL_DIVERGE = 1.0D+4 
      MAX_INLET_VEL_FAC = ONE
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
      LEQ_METHOD(2) = 2 
      LEQ_METHOD(3) = 2 
      LEQ_METHOD(4) = 2 
      LEQ_METHOD(5) = 2 
      LEQ_METHOD(6) = 2 
      LEQ_METHOD(7) = 2 
      LEQ_METHOD(8) = 2 
      LEQ_METHOD(9) = 2 
      LEQ_SWEEP(1) = 'RSRS' 
      LEQ_SWEEP(2) = 'RSRS'  
      LEQ_SWEEP(3) = 'RSRS' 
      LEQ_SWEEP(4) = 'RSRS'  
      LEQ_SWEEP(5) = 'RSRS'  
      LEQ_SWEEP(6) = 'RSRS'  
      LEQ_SWEEP(7) = 'RSRS' 
      LEQ_SWEEP(8) = 'RSRS'  
      LEQ_SWEEP(9) = 'RSRS'  
      LEQ_TOL(1) = 1.0D-4 
      LEQ_TOL(2) = 1.0D-4  
      LEQ_TOL(3) = 1.0D-4 
      LEQ_TOL(4) = 1.0D-4  
      LEQ_TOL(5) = 1.0D-4  
      LEQ_TOL(6) = 1.0D-4  
      LEQ_TOL(7) = 1.0D-4 
      LEQ_TOL(8) = 1.0D-4  
      LEQ_TOL(9) = 1.0D-4  
      LEQ_PC(1:9)  = 'LINE'
      DO_TRANSPOSE = .FALSE.
      icheck_bicgs = 1
      solver_statistics = .FALSE.
      opt_parallel = .FALSE.
! AEOLUS: set default value for new variable to debug print whole 
!         index layout 
      DBGPRN_LAYOUT = .FALSE.
! AEOLUS: set default value for enabling all processors write out
!         their *.LOG invidually 
      ENABLE_DMP_LOG = .FALSE.
      UR_FAC(1)  = 0.8D0             !pressure 
      UR_FAC(2)  = 0.5D0             !rho, ep 
      UR_FAC(3)  = 0.5D0             !U 
      UR_FAC(4)  = 0.5D0             !V 
      UR_FAC(5)  = 0.5D0             !W 
      UR_FAC(6)  = 1.0D0             !T 
      UR_FAC(7)  = 1.0D0             !X 
      UR_FAC(8)  = 0.5D0             !Th 
      UR_FAC(9)  = 0.8D0             !Scalar
      UR_F_gs    = 1.0D0             !drag coefficient update
      UR_Kth_sml = 1.0D0            ! conductivity term in IA theory

! INITIALIZE THE GAS PHASE SECTION
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
      D_P0(:DIM_M) = UNDEFINED 
      RO_S(:DIM_M) = UNDEFINED 
      NMAX(1:DIM_M) = UNDEFINED_I 
      CLOSE_PACKED(:DIM_M) = .TRUE. 
      MW_S = UNDEFINED
      EP_STAR = UNDEFINED 

!QX
!     density of each component
      RO_SS = UNDEFINED
      SOLID_RO_V = .FALSE.
!end

      NMAX_g = UNDEFINED_I
      SPECIES_g(:) = UNDEFINED_C
      SPECIES_ALIAS_g(:) = UNDEFINED_C

      NMAX_s(:) = UNDEFINED_I
      SPECIES_s(:,:) = UNDEFINED_C
      SPECIES_ALIAS_s(:,:) = UNDEFINED_C


      USE_RRATES = .FALSE.
! NO_OF_RXNS is not a keyword. However, it is initialized here so that
! if there are no reactions, this value is assigned.
      NO_OF_RXNS = UNDEFINED_I

! INITIALIZE THE INITIAL CONDITIONS
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
         IC_SCALAR(LC,:DIM_SCALAR) = UNDEFINED 

! sof: force users to set initial values for K and Epsilon. 
         IC_K_Turb_G(LC) = UNDEFINED  
         IC_E_Turb_G(LC) = UNDEFINED

         IC_GAMA_RS(LC,:DIM_M) = ZERO 
         IC_T_RS(LC,:DIM_M) = UNDEFINED 
!         IC_X_S(LC,1,1+:DIM_N_S+) = UNDEFINED 
      ENDDO 

      IC_X_S = UNDEFINED
      DELP_X = UNDEFINED 
      DELP_Y = UNDEFINED 
      DELP_Z = UNDEFINED 
      Flux_g = UNDEFINED 
      U_G0 = UNDEFINED 
      V_G0 = UNDEFINED 
      W_G0 = UNDEFINED 
      U_S0(:DIM_M) = UNDEFINED 
      V_S0(:DIM_M) = UNDEFINED 
      W_S0(:DIM_M) = UNDEFINED 

! Boundary Conditions      
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
         BC_VELMAG_G(LC) = UNDEFINED 
         BC_TYPE(LC) = UNDEFINED_C 
         BC_APPLY_TO_MPPIC(LC) = .true.
         BC_VOLFLOW_G(LC) = UNDEFINED 
         BC_MASSFLOW_G(LC) = UNDEFINED 
         BC_DT_0(LC) = UNDEFINED 
         BC_DT_H(LC) = UNDEFINED 
         BC_DT_L(LC) = UNDEFINED 
         BC_JET_G0(LC) = UNDEFINED 
         BC_JET_GH(LC) = UNDEFINED 
         BC_JET_GL(LC) = UNDEFINED 

         BC_HW_G(LC) = UNDEFINED 
         BC_UW_G(LC) = UNDEFINED 
         BC_VW_G(LC) = UNDEFINED 
         BC_WW_G(LC) = UNDEFINED 

         BC_HW_T_G(LC) = UNDEFINED 
         BC_TW_G(LC) = UNDEFINED 
         BC_C_T_G(LC) = UNDEFINED 

! Kapil & Anuj: 01/19/98
         BC_JJ_PS(LC) = UNDEFINED_I 

         BC_ROP_S(LC,:DIM_M) = UNDEFINED 
         BC_U_S(LC,:DIM_M) = UNDEFINED 
         BC_V_S(LC,:DIM_M) = UNDEFINED 
         BC_W_S(LC,:DIM_M) = UNDEFINED 
         BC_VELMAG_S(LC,:DIM_M) = UNDEFINED 
         BC_T_S(LC,:DIM_M) = UNDEFINED 
         BC_VOLFLOW_S(LC,:DIM_M) = UNDEFINED 
         BC_MASSFLOW_S(LC,:DIM_M) = UNDEFINED 

         BC_HW_S(LC,:DIM_M) = UNDEFINED 
         BC_UW_S(LC,:DIM_M) = UNDEFINED 
         BC_VW_S(LC,:DIM_M) = UNDEFINED 
         BC_WW_S(LC,:DIM_M) = UNDEFINED 

         BC_HW_T_S(LC,:DIM_M) = UNDEFINED 
         BC_TW_S(LC,:DIM_M) = UNDEFINED 
         BC_C_T_S(LC,:DIM_M) = UNDEFINED 

         BC_HW_THETA_M(LC,:DIM_M) = UNDEFINED 
         BC_THETAW_M(LC,:DIM_M) = UNDEFINED 
         BC_C_THETA_M(LC,:DIM_M) = UNDEFINED 

         BC_HW_Scalar(LC,:DIM_SCALAR) = UNDEFINED 
         BC_ScalarW(LC,:DIM_SCALAR) = UNDEFINED 
         BC_C_Scalar(LC,:DIM_SCALAR) = UNDEFINED 
      ENDDO 
   
      BC_THETA_M = UNDEFINED 
      BC_Scalar = UNDEFINED

! sof: force users to set inlet BC for K and Epsilon  
      BC_K_Turb_G = UNDEFINED  
      BC_E_Turb_G = UNDEFINED 

      BC_X_S = UNDEFINED 
      BC_HW_X_S = UNDEFINED 
      BC_XW_S = UNDEFINED 
      BC_C_X_S = UNDEFINED 

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
      ENDDO 

      DO LC = 1, DIM_N_ALL
         SPECIES_NAME(LC) = UNDEFINED_C 
      ENDDO 

      NODESI = UNDEFINED_I
      NODESJ = UNDEFINED_I
      NODESK = UNDEFINED_I

      IS_SERIAL = .TRUE.
      USE_DOLOOP = .FALSE.

! Stiff Chemistry Solver:
      STIFF_CHEMISTRY = .FALSE.
      CALL_DI   = .FALSE.      ! Legacy Keyword
      CALL_GROW = .FALSE.      ! Legacy Keyword
      CALL_ISAT = .FALSE.      ! Legacy Keyword
      ISATdt    = UNDEFINED    ! Legacy Keyword

      bWrite_netCDF(:) = .false.

      CALL DES_INIT_NAMELIST

! Initialize QMOMK namelist
      CALL QMOMK_INIT_NAMELIST

      CALL USR_INIT_NAMELIST 

! JFD: MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION      
      CALL CARTESIAN_GRID_INIT_NAMELIST


      RETURN  
      END SUBROUTINE INIT_NAMELIST 
