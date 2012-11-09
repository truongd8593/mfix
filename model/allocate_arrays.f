!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                      
!  Module name: ALLOCATE_ARRAYS                                     
!  Purpose: allocate arrays
!                                                                      
!  Author: M. Syamlal                                Date: 17-DEC-98 
!  Reviewer: 
!                                                                     
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE ALLOCATE_ARRAYS 

!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1
      Use ambm
      Use coeff
      Use cont
      Use drag
      Use energy
      Use fldvar
      Use geometry
      Use indices
      Use pgcor
      Use physprop
      Use pscor
      Use residual
      Use rxns
      Use run
      Use scalars
      Use turb
      Use tau_g
      Use tau_s
      Use tmp_array
      Use tmp_array1
      Use trace
      Use visc_g
      Use visc_s
      Use xsi_array
      Use vshear
      Use mflux
      Use mchem
      USE ghdtheory
      use kintheory
      Use cdist
      Use des_rxns

      IMPLICIT NONE

!-----------------------------------------------     
! Variables     
!-----------------------------------------------     
      INTEGER M

      integer :: dimension_3p   ! used during post_mfix to reduce allocations
!-----------------------------------------------     

!// Modified the DIMENSION_X based on the new domain decomposition variables
      DIMENSION_I   = IMAX3
      DIMENSION_J   = JMAX3
      DIMENSION_K   = KMAX3
      DIMENSION_3   = (kend3-kstart3+1)*(jend3-jstart3+1)*(iend3-istart3+1)
      DIMENSION_3G   = IJKMAX3            
      DIMENSION_3L  = ijksize3_all(myPE)      
      DIMENSION_M   = MAX(1, MMAX)
      DIMENSION_4   = (kend4-kstart4+1)*(jend4-jstart4+1)*(iend4-istart4+1)

      DIMENSION_N_g = 1
      IF(USE_RRATES) THEN
          IF(NMAX(0) .NE. UNDEFINED_I)DIMENSION_N_g = NMAX(0)
      ELSE
          IF(NMAX_g .NE. UNDEFINED_I)DIMENSION_N_g = NMAX_g
      ENDIF

      
! to reduce allocation space when doing post_mfix
      if (bDoing_postmfix) then
         dimension_3p = 1
      else
         dimension_3p = dimension_3
      endif

      DIMENSION_N_s = 1
      DO M = 1, MMAX
         IF(USE_RRATES) THEN
            IF(NMAX(M) .NE. UNDEFINED_I) &
               DIMENSION_N_s = MAX(DIMENSION_N_s, NMAX(M))
         ELSE
            IF(NMAX_s(M) .NE. UNDEFINED_I) &
               DIMENSION_N_s = MAX(DIMENSION_N_s, NMAX_s(M))
         ENDIF
      ENDDO
      DO M = 1, DIM_M
         IF(DES_NMAX_s(M) .NE. UNDEFINED_I) &
            DIMENSION_N_s = MAX(DIMENSION_N_s, DES_NMAX_s(M))
      ENDDO

      DIMENSION_LM    = (DIMENSION_M * (DIMENSION_M-1) / 2)+1
      DIMENSION_N_all = max(DIMENSION_N_g, DIMENSION_N_s)
      DIMENSION_Scalar = NScalar
!add by rong
      DIM_Scalar2 = 2*NScalar

!ambm
      Allocate( A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) )
      Allocate( B_m(DIMENSION_3, 0:DIMENSION_M) )

!coeff
      
      Allocate( DENSITY(0:DIMENSION_M) )
      Allocate( PSIZE(0:DIMENSION_M) )
      Allocate( SP_HEAT(0:DIMENSION_M) )
      Allocate( VISC(0:DIMENSION_M) )
      Allocate( COND(0:DIMENSION_M) )
      Allocate( DIFF(0:DIMENSION_M) )
      Allocate( DRAGCOEF(0:DIMENSION_M, 0:DIMENSION_M) )
      Allocate( HEAT_TR(0:DIMENSION_M, 0:DIMENSION_M))
      Allocate( GRAN_DISS(0:DIMENSION_M) )

!cont
      Allocate( DO_CONT(0:DIMENSION_M) )

!drag
      Allocate(  F_gs(DIMENSION_3, DIMENSION_M) )
      Allocate(  F_ss(DIMENSION_3, 0:DIMENSION_LM) )

!Off diagonal friction coefficient in HYS drag relation
      IF(TRIM(DRAG_TYPE).EQ.'HYS') &
         Allocate(  beta_ij(DIMENSION_3, 0:DIMENSION_M, 0:DIMENSION_M) )

!energy
      Allocate(  HOR_g (DIMENSION_3) )
      Allocate(  HOR_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  GAMA_gs (DIMENSION_3, DIMENSION_M) )
      Allocate(  GAMA_Rg (DIMENSION_3) )
      Allocate(  GAMA_Rs (DIMENSION_3, DIMENSION_M) )
      Allocate(  T_Rg (DIMENSION_3) )
      Allocate(  T_Rs (DIMENSION_3, DIMENSION_M) )

!fldvar
      Allocate(  EP_g (DIMENSION_3) )
      Allocate(  EP_go (DIMENSION_3p) )
      Allocate(  P_g (DIMENSION_3) )
      Allocate(  P_go (DIMENSION_3p) )
      Allocate(  RO_g (DIMENSION_3) )
      Allocate(  RO_go (DIMENSION_3p) )
      Allocate(  ROP_g (DIMENSION_3) )
      Allocate(  ROP_go (DIMENSION_3p) )
      Allocate(  ROP_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  ROP_so (DIMENSION_3p, DIMENSION_M) )
      Allocate(  T_g (DIMENSION_3) )
      Allocate(  T_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  T_go (DIMENSION_3p) )
      Allocate(  T_so (DIMENSION_3p, DIMENSION_M) )
      Allocate(  X_g (DIMENSION_3, DIMENSION_N_g) )
      Allocate(  X_s (DIMENSION_3, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  X_go (DIMENSION_3p, DIMENSION_N_g) )
      Allocate(  X_so (DIMENSION_3p, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  U_g (DIMENSION_3) )
      Allocate(  U_go (DIMENSION_3p) )
      Allocate(  U_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  U_so (DIMENSION_3p, DIMENSION_M) )
      Allocate(  V_g (DIMENSION_3) )
      Allocate(  V_go (DIMENSION_3p) )
      Allocate(  V_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  V_so (DIMENSION_3p, DIMENSION_M) )
      Allocate(  W_g (DIMENSION_3) )
      Allocate(  W_go (DIMENSION_3p) )
      Allocate(  W_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  W_so (DIMENSION_3p, DIMENSION_M) )
      Allocate(  P_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  P_s_c (DIMENSION_3, DIMENSION_M) )
      Allocate(  P_s_v (DIMENSION_3) )
      Allocate(  P_s_f (DIMENSION_3) )
      Allocate(  P_s_p (DIMENSION_3) )
      Allocate(  P_star (DIMENSION_3) )
      Allocate(  P_staro (DIMENSION_3p) )
      Allocate(  THETA_m (DIMENSION_3, DIMENSION_M) )
      Allocate(  THETA_mo (DIMENSION_3p, DIMENSION_M) )

! sof: MUST use k-epsilon model if using Simonin or Ahmadi model 
      IF(K_Epsilon .OR. SIMONIN .OR. AHMADI) THEN
        Allocate(  K_Turb_G (DIMENSION_3) )
        Allocate(  K_Turb_Go (DIMENSION_3p) )
        Allocate(  E_Turb_G (DIMENSION_3) )
        Allocate(  E_Turb_Go (DIMENSION_3p) )
      ENDIF
      
      IF(DIMENSION_Scalar /= 0) THEN
        Allocate(  Scalar (DIMENSION_3,  DIMENSION_Scalar) )
        Allocate(  Scalaro (DIMENSION_3p, DIMENSION_Scalar) )
      ENDIF


!geometry
      Allocate(  FLAG (DIMENSION_3) )
      Allocate(  FLAG_E (DIMENSION_3) )
      Allocate(  FLAG_N (DIMENSION_3) )
      Allocate(  FLAG_T (DIMENSION_3) )
      Allocate(  ICBC_FLAG (DIMENSION_3L) )
      Allocate(  oDX (0:DIMENSION_I) )
      Allocate(  oDY (0:DIMENSION_J) )
      Allocate(  oDZ (0:DIMENSION_K) )
      Allocate(  oDX_E (0:DIMENSION_I) )
      Allocate(  oDY_N (0:DIMENSION_J) )
      Allocate(  oDZ_T (0:DIMENSION_K) )
      Allocate(  X (0:DIMENSION_I) )
      Allocate(  X_E (0:DIMENSION_I) )
      Allocate(  oX (0:DIMENSION_I) )
      Allocate(  oX_E (0:DIMENSION_I) )
      Allocate(  Z (0:DIMENSION_K) )
      Allocate(  Z_T (0:DIMENSION_K) )
      Allocate(  FX (0:DIMENSION_I) )
      Allocate(  FX_bar (0:DIMENSION_I) )
      Allocate(  FX_E (0:DIMENSION_I) )
      Allocate(  FX_E_bar (0:DIMENSION_I) )
      Allocate(  FY_N (0:DIMENSION_J) )
      Allocate(  FY_N_bar (0:DIMENSION_J) )
      Allocate(  FZ_T (0:DIMENSION_K) )
      Allocate(  FZ_T_bar (0:DIMENSION_K) )
      Allocate(  AYZ (DIMENSION_3p) )
      Allocate(  AXZ (DIMENSION_3p) )
      Allocate(  AXY (DIMENSION_3p) )
      Allocate(  VOL (DIMENSION_3) )
      Allocate(  AYZ_U (DIMENSION_3p) )
      Allocate(  AXZ_U (DIMENSION_3p) )
      Allocate(  AXY_U (DIMENSION_3p) )
      Allocate(  VOL_U (DIMENSION_3) )
      Allocate(  AYZ_V (DIMENSION_3p) )
      Allocate(  AXZ_V (DIMENSION_3p) )
      Allocate(  AXY_V (DIMENSION_3p) )
      Allocate(  VOL_V (DIMENSION_3) )
      Allocate(  AYZ_W (DIMENSION_3p) )
      Allocate(  AXZ_W (DIMENSION_3p) )
      Allocate(  AXY_W (DIMENSION_3p) )
      Allocate(  VOL_W (DIMENSION_3) )

!indices
      Allocate(  STORE_LM (DIMENSION_M, DIMENSION_M) )
      Allocate(  CELL_CLASS (DIMENSION_3) )
      Allocate(  I_OF (DIMENSION_3) )
      Allocate(  J_OF (DIMENSION_3) )
      Allocate(  K_OF (DIMENSION_3) )
      Allocate(  Im1 (0:DIMENSION_I) )
      Allocate(  Ip1 (0:DIMENSION_I) )
      Allocate(  Jm1 (0:DIMENSION_J) )
      Allocate(  Jp1 (0:DIMENSION_J) )
      Allocate(  Km1 (0:DIMENSION_K) )
      Allocate(  Kp1 (0:DIMENSION_K) )
      
!pgcor
      Allocate(  d_e(DIMENSION_3p, 0:DIMENSION_M) )
      Allocate(  d_n(DIMENSION_3p, 0:DIMENSION_M) )
      Allocate(  d_t(DIMENSION_3p, 0:DIMENSION_M) )
      Allocate(  Pp_g(DIMENSION_3p) )
      Allocate(  PHASE_4_P_g(DIMENSION_3p) )

!physprop
      Allocate(  MU_g (DIMENSION_3) )
      Allocate(  C_pg (DIMENSION_3) )
      Allocate(  C_ps (DIMENSION_3, DIMENSION_M) )
      Allocate(  K_g (DIMENSION_3) )
      Allocate(  K_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  Kth_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  Kphi_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DIF_g (DIMENSION_3p, DIMENSION_N_g) )
      Allocate(  DIF_s (DIMENSION_3p, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  MW_MIX_g (DIMENSION_3) )

!pscor
      Allocate(  e_e(DIMENSION_3p) )
      Allocate(  e_n(DIMENSION_3p) )
      Allocate(  e_t(DIMENSION_3p) )
      Allocate(  K_cp(DIMENSION_3p) )
      Allocate(  EPp(DIMENSION_3p) )
      Allocate(  PHASE_4_P_s(DIMENSION_3p) )

!residual
      Allocate( RESID(NRESID, 0:DIMENSION_M) )
      Allocate( MAX_RESID(NRESID, 0:DIMENSION_M) )
      Allocate( IJK_RESID(NRESID, 0:DIMENSION_M) )
      Allocate( NUM_RESID(NRESID, 0:DIMENSION_M) )
      Allocate( DEN_RESID(NRESID, 0:DIMENSION_M) )
      Allocate( RESID_PACK(NRESID*2*(DIMENSION_M+1)))
 
!rxns
      if (nRR .gt. 0) Allocate( ReactionRates(DIMENSION_3,nRR) )
      Allocate(  R_gp (DIMENSION_3p, DIMENSION_N_g) )
      Allocate(  R_sp (DIMENSION_3p, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  RoX_gc (DIMENSION_3p, DIMENSION_N_g) )
      Allocate(  RoX_sc (DIMENSION_3p, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  SUM_R_g (DIMENSION_3p) )
      Allocate(  SUM_R_s (DIMENSION_3p, DIMENSION_M) )
      Allocate(  R_phase (DIMENSION_3, DIMENSION_LM+DIMENSION_M-1) )

! Undefined indicates that no reaction block was found in the deck file.
      IF(NO_OF_RXNS .NE. UNDEFINED_I) &
         Allocate( REACTION( NO_OF_RXNS ))
      
!scalars
      IF(DIMENSION_Scalar /= 0) then
        Allocate(  Scalar_c (DIMENSION_3p,  DIMENSION_Scalar) )
        Allocate(  Scalar_p (DIMENSION_3p,  DIMENSION_Scalar) )
        Allocate(  Dif_Scalar (DIMENSION_3p, DIMENSION_Scalar) )
      ENDIF

! isat
!  Insert user-defined code here
      Allocate( N_sh (DIMENSION_3, DIMENSION_M) )

! add by rong for dqmom
      Allocate(  D_p  (DIMENSION_3, DIMENSION_M) )
      Allocate(  D_po (DIMENSION_3, DIMENSION_M) )
!      Allocate(  ome  (DIMENSION_3, DIMENSION_M) )
!      Allocate(  ome_o (DIMENSION_3, DIMENSION_M) )
      Allocate(  Source_a(DIMENSION_3, DIMENSION_M) )
      Allocate(  S_bar( 0:DIM_Scalar2-1 ))
      Allocate(  Matrix_a(DIM_Scalar2,DIM_scalar2))
      Allocate(  Matrix_b(DIM_Scalar2,DIM_scalar2))
      Allocate(  Matrix_c(DIM_Scalar2,DIM_scalar2))
      Allocate(  Inv_a(DIM_Scalar2,DIM_scalar2))
      Allocate(  A( 1:DIMENSION_Scalar))
      Allocate(  omega(1:DIMENSION_m))
      ALLocate(  beta_a( DIM_Scalar,DIM_Scalar))
      ALLocate(  ystart( 1:DIM_Scalar2))
!     ALLocate(  g_a( 1:DIMENSION_Scalar))
 

! K-Epsilon Turbulence model
! sof (02/01/05): must use k-epsilon model with Simonin or Ahmadi models
      IF(K_Epsilon .OR. SIMONIN .OR. AHMADI)then
        Allocate(  K_Turb_G_c   (DIMENSION_3p) )
        Allocate(  K_Turb_G_p   (DIMENSION_3p) )
        Allocate(  Dif_K_Turb_G (DIMENSION_3p) )
        Allocate(  E_Turb_G_c   (DIMENSION_3p) )
        Allocate(  E_Turb_G_p   (DIMENSION_3p) )
        Allocate(  Dif_E_Turb_G (DIMENSION_3p) )
      ENDIF

! Simonin or Ahmadi model
      IF(SIMONIN .OR. AHMADI)then
        Allocate(  K_12 (DIMENSION_3) )
        Allocate(  Tau_12 (DIMENSION_3) )
        Allocate(  Tau_1 (DIMENSION_3) )
        Allocate(  Cos_theta (DIMENSION_3) )
      ENDIF

!tau_g
      Allocate(  TAU_U_g(DIMENSION_3p) )
      Allocate(  TAU_V_g(DIMENSION_3p) )
      Allocate(  TAU_W_g(DIMENSION_3p) )

!tau_s
      Allocate(  TAU_U_s(DIMENSION_3p, DIMENSION_M) )
      Allocate(  TAU_V_s(DIMENSION_3p, DIMENSION_M) )
      Allocate(  TAU_W_s(DIMENSION_3p, DIMENSION_M) )
      
!tmp_array
      Allocate(  Array1(DIMENSION_3) )
      Allocate(  Array2(DIMENSION_3) )
      Allocate(  Array3(DIMENSION_3) )
      Allocate(  Array4(DIMENSION_3) )
      Allocate(  Array1i(DIMENSION_3) ) 
      Allocate(  Array1c(DIMENSION_3) )
      Allocate(  TMP4(DIMENSION_4) )
      Allocate(  ArrayLM (DIMENSION_3, DIMENSION_LM) )  !S. Dartevelle, LANL, Feb. 2004


!tmp_array1
      Allocate(  Arraym1(DIMENSION_3, DIMENSION_M) )

!trace
      Allocate(  trD_s_C (DIMENSION_3, DIMENSION_M) )
      Allocate(  trD_s2 (DIMENSION_3, DIMENSION_M) )
      Allocate(  trD_s_Co (DIMENSION_3, DIMENSION_M) )

!visc_g
      Allocate(  trD_g(DIMENSION_3) )
      Allocate(  MU_gt (DIMENSION_3) )
      Allocate(  LAMBDA_gt (DIMENSION_3) )
      Allocate(  L_scale (DIMENSION_3) )

!visc_s
      Allocate(  trD_s(DIMENSION_3, DIMENSION_M) )
      Allocate(  MU_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  LAMBDA_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  ALPHA_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  MU_s_c (DIMENSION_3, DIMENSION_M) )
      Allocate(  LAMBDA_s_c (DIMENSION_3, DIMENSION_M) )
      Allocate(  LAMBDA_s_v (DIMENSION_3) )
      Allocate(  LAMBDA_s_f (DIMENSION_3) )
      Allocate(  LAMBDA_s_p (DIMENSION_3) )
      Allocate(  MU_s_v (DIMENSION_3) )
      Allocate(  MU_s_f (DIMENSION_3) )
      Allocate(  MU_s_p (DIMENSION_3) )
      Allocate(  MU_b_v (DIMENSION_3) )
      Allocate(  ALPHA_s_v (DIMENSION_3) )
      Allocate(  ALPHA_s_p (DIMENSION_3) )
      Allocate(  ALPHA_s_f (DIMENSION_3) )
      Allocate(  EP_star_array (DIMENSION_3) )
      Allocate(  EP_g_blend_start (DIMENSION_3) )
      Allocate(  EP_g_blend_end (DIMENSION_3) )
      Allocate(  VREL_array (DIMENSION_3) )
      Allocate(  I2_devD_s (DIMENSION_3) )
      Allocate(  TrM_s (DIMENSION_3) )
      Allocate(  TrDM_s (DIMENSION_3) )
      
!xsi_array
      Allocate(  Xsi_e(DIMENSION_3) )
      Allocate(  Xsi_n(DIMENSION_3) )
      Allocate(  Xsi_t(DIMENSION_3) )

!shear quantities
      Allocate(  VSH(DIMENSION_3) )
      Allocate(  VSHE(DIMENSION_3) )


! array allocation of add on packages, such as linear equation solvers


! array allocation for higher order implementation
      Allocate( FLAG3 (DIMENSION_4) )
      Allocate( CELL_CLASS3 (DIMENSION_4) )
      Allocate( I3_OF (DIMENSION_4) )
      Allocate( J3_OF (DIMENSION_4) )
      Allocate( K3_OF (DIMENSION_4) )
      Allocate( Im1_3 (-1:DIMENSION_I+1) )
      Allocate( Ip1_3 (-1:DIMENSION_I+1) )
      Allocate( Jm1_3 (-1:DIMENSION_J+1) )
      Allocate( Jp1_3 (-1:DIMENSION_J+1) )
      Allocate( Km1_3 (-1:DIMENSION_K+1) )
      Allocate( Kp1_3 (-1:DIMENSION_K+1) )
 
!mflux
      Allocate( Flux_gE(DIMENSION_3p) )
      Allocate( Flux_sE(DIMENSION_3p, DIMENSION_M) )
      Allocate( Flux_gN(DIMENSION_3p) )
      Allocate( Flux_sN(DIMENSION_3p, DIMENSION_M) )
      Allocate( Flux_gT(DIMENSION_3p) )
      Allocate( Flux_sT(DIMENSION_3p, DIMENSION_M) )
      IF(ADDED_MASS) THEN ! Fluxes calculated for just one 'bubble' species (M=M_AM)
         Allocate( Flux_gSE(DIMENSION_3p) )
         Allocate( Flux_sSE(DIMENSION_3p) )
         Allocate( Flux_gSN(DIMENSION_3p) )
         Allocate( Flux_sSN(DIMENSION_3p) )
         Allocate( Flux_gST(DIMENSION_3p) )
         Allocate( Flux_sST(DIMENSION_3p) )
      ENDIF
      Allocate( ROP_gE(DIMENSION_3p) )
      Allocate( ROP_sE(DIMENSION_3p, DIMENSION_M) )
      Allocate( ROP_gN(DIMENSION_3p) )
      Allocate( ROP_sN(DIMENSION_3p, DIMENSION_M) )
      Allocate( ROP_gT(DIMENSION_3p) )
      Allocate( ROP_sT(DIMENSION_3p, DIMENSION_M) )
 
! allocate variables for GHD Theory
      IF (TRIM(KT_TYPE) == 'GHD') THEN 
        Allocate(  Flux_nE(DIMENSION_3p) ) 
        Allocate(  Flux_nN(DIMENSION_3p) ) 
        Allocate(  Flux_nT(DIMENSION_3p) )  
        Allocate(  Zeta0(DIMENSION_3p) )   ! zeroth rate of cooling  
        Allocate(  ZetaU(DIMENSION_3p) )   ! 1st order cooling rate transport coefficient
        Allocate(  DiT(DIMENSION_3p, DIMENSION_M) )   ! thermal diffusivity
        Allocate(  DijF(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )   ! mass mobility
        Allocate(  Lij(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )   ! thermal mobility
        Allocate(  Dij(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )   ! ordinary diffusion
        Allocate(  DijQ(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )   ! Dufour coeff.
        Allocate(  JoiX(DIMENSION_3p, DIMENSION_M) )   ! X- species mass flux
        Allocate(  JoiY(DIMENSION_3p, DIMENSION_M) )   ! Y- species mass flux
        Allocate(  JoiZ(DIMENSION_3p, DIMENSION_M) )   ! Z- species mass flux
        Allocate(  FiX(DIMENSION_3p, DIMENSION_M) )   ! X- external force
        Allocate(  FiY(DIMENSION_3p, DIMENSION_M) )   ! Y- external force
        Allocate(  FiZ(DIMENSION_3p, DIMENSION_M) )   ! Z- external force
        Allocate(  FiXvel(DIMENSION_3p, DIMENSION_M) )   ! X- external force
        Allocate(  FiYvel(DIMENSION_3p, DIMENSION_M) )   ! Y- external force
        Allocate(  FiZvel(DIMENSION_3p, DIMENSION_M) )   ! Z- external force
        Allocate(  DELTAU(DIMENSION_3p, DIMENSION_M) )   
        Allocate(  DELTAV(DIMENSION_3p, DIMENSION_M) )   
        Allocate(  DELTAW(DIMENSION_3p, DIMENSION_M) )   
        Allocate(  dragFx(DIMENSION_3p, DIMENSION_M) )   ! X- drag force
        Allocate(  dragFy(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  dragFz(DIMENSION_3p, DIMENSION_M) )   ! Z- drag force
        Allocate(  dragFxflux(DIMENSION_3p, DIMENSION_M) )   ! X- drag force
        Allocate(  dragFyflux(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  dragFzflux(DIMENSION_3p, DIMENSION_M) )   ! Z- drag force
        Allocate(  FiMinusDragX(DIMENSION_3p, DIMENSION_M) )   ! X- drag force
        Allocate(  JoiMinusDragX(DIMENSION_3p, DIMENSION_M) )   ! X- drag force
        Allocate(  FiMinusDragY(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  JoiMinusDragY(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  FiMinusDragZ(DIMENSION_3p, DIMENSION_M) )   ! Z- drag force
        Allocate(  JoiMinusDragZ(DIMENSION_3p, DIMENSION_M) )   ! Z- drag force
        Allocate(  beta_cell_X(DIMENSION_3p, DIMENSION_M) )   ! X- drag force
        Allocate(  beta_cell_Y(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  beta_cell_Z(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  beta_ij_cell_X(DIMENSION_3p, DIMENSION_M,DIMENSION_M) )   ! X- drag force
        Allocate(  beta_ij_cell_Y(DIMENSION_3p, DIMENSION_M,DIMENSION_M) )   ! Y- drag force
        Allocate(  beta_ij_cell_Z(DIMENSION_3p, DIMENSION_M,DIMENSION_M) )   ! Y- drag force
        Allocate(  DEL_DOT_J(DIMENSION_3p, DIMENSION_M) )   
        Allocate(  DiT_HarmE(DIMENSION_3p) )   
        Allocate(  DiT_HarmN(DIMENSION_3p) )   
        Allocate(  DiT_HarmT(DIMENSION_3p) )   
        Allocate(  Dij_HarmE(DIMENSION_3p, DIMENSION_M) )   
        Allocate(  Dij_HarmN(DIMENSION_3p, DIMENSION_M) )   
        Allocate(  Dij_HarmT(DIMENSION_3p, DIMENSION_M) )   
        Allocate(  DijF_HarmE(DIMENSION_3p, DIMENSION_M) )   
        Allocate(  DijF_HarmN(DIMENSION_3p, DIMENSION_M) )   
        Allocate(  DijF_HarmT(DIMENSION_3p, DIMENSION_M) )   
      ENDIF


! We need to set this even when KT_TYPE is not set to IA_NONEP - at
! least in the current version of the code and needs to be revisited
      Allocate(  KTMOM_U_s(DIMENSION_3p, DIMENSION_M) )
      Allocate(  KTMOM_V_s(DIMENSION_3p, DIMENSION_M) )
      Allocate(  KTMOM_W_s(DIMENSION_3p, DIMENSION_M) )

! allocate variables for Iddir & Arastoopour (2005) kinetic theory
! EDvel_sM_ip & EDT_s_ip are also used for Garzy & Dufty (1999) kinetic theory
      IF (TRIM(KT_TYPE) == 'IA_NONEP') THEN      
         Allocate(  trD_s2_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  MU_sM_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  MU_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  XI_sM_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  XI_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  Fnu_s_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  FT_sM_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  FT_sL_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  Kth_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  Knu_sM_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  Knu_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  Kvel_s_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  EDvel_sL_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  ED_ss_ip(DIMENSION_3p, 0:DIMENSION_LM) )
      ENDIF
      IF (TRIM(KT_TYPE) == 'IA_NONEP' .OR. TRIM(KT_TYPE) == 'GD_99') THEN
         Allocate(  EDT_s_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  EDvel_sM_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
      ENDIF


      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS 
      

