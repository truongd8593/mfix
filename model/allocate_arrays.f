! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutinee: ALLOCATE_ARRAYS                                        C
!  Purpose: allocate arrays                                            C
!                                                                      C
!  Author: M. Syamlal                                Date: 17-DEC-98   C
!  Reviewer:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ALLOCATE_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use ambm
      use cdist
      use cont
      use des_rxns
      use drag
      use energy
      use fldvar
      use generate_particles, only: particle_count
      use geometry
      use ghdtheory
      use indices
      use kintheory
      use mflux
      use param
      use param1
      use pgcor
      use physprop
      use pscor
      use residual
      use run
      use rxns
      use scalars
      use iterate, only: errorpercent
      use tau_g
      use tau_s
      use trace
      use turb
      use visc_g
      use visc_s
      use vshear

      IMPLICIT NONE

!-----------------------------------------------
! Variables
!-----------------------------------------------

!ambm
      Allocate( A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) )
      Allocate( B_m(DIMENSION_3, 0:DIMENSION_M) )

!cont
      Allocate( DO_CONT(0:DIMENSION_M) )

!drag
      Allocate(  F_gs(DIMENSION_3, DIMENSION_M) )
      Allocate(  F_ss(DIMENSION_3, 0:DIMENSION_LM) )

!Off diagonal friction coefficient in HYS drag relation
      IF(DRAG_TYPE_ENUM.EQ.HYS) &
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
      Allocate(  epg_jfac (DIMENSION_3p) )
      Allocate(  epg_ifac (DIMENSION_3p) )
      Allocate(  eps_ifac (DIMENSION_3p, DIMENSION_M) )
      Allocate(  EP_go (DIMENSION_3p) )
      Allocate(  P_g (DIMENSION_3) )
      Allocate(  P_go (DIMENSION_3p) )
      Allocate(  RO_g (DIMENSION_3) )
      Allocate(  RO_go (DIMENSION_3p) )
      Allocate(  ROP_g (DIMENSION_3) )
      Allocate(  ROP_go (DIMENSION_3p) )
      Allocate(  RO_S (DIMENSION_3, DIMENSION_M) )
      Allocate(  RO_So (DIMENSION_3p, DIMENSION_M) )
      Allocate(  ROP_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  ROP_so (DIMENSION_3p, DIMENSION_M) )

      Allocate(  EP_SS(DIMENSION_3,DIMENSION_M,DIMENSION_N_S) )
      Allocate(  ERR_ARRAY(DIMENSION_3,DIMENSION_M) )

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

      IF(K_Epsilon)THEN
        Allocate(  K_Turb_G (DIMENSION_3) )
        Allocate(  K_Turb_Go (DIMENSION_3p) )
        Allocate(  E_Turb_G (DIMENSION_3) )
        Allocate(  E_Turb_Go (DIMENSION_3p) )
      ENDIF

      IF(DIMENSION_Scalar /= 0) THEN
        Allocate(  Scalar (DIMENSION_3,  DIMENSION_Scalar) )
        Allocate(  Scalaro (DIMENSION_3p, DIMENSION_Scalar) )
      ENDIF


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

!scalars
      IF(DIMENSION_Scalar /= 0) then
        Allocate(  Scalar_c (DIMENSION_3p,  DIMENSION_Scalar) )
        Allocate(  Scalar_p (DIMENSION_3p,  DIMENSION_Scalar) )
        Allocate(  Dif_Scalar (DIMENSION_3p, DIMENSION_Scalar) )
      ENDIF

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
      IF(K_Epsilon) THEN
        Allocate(  K_Turb_G_c   (DIMENSION_3p) )
        Allocate(  K_Turb_G_p   (DIMENSION_3p) )
        Allocate(  Dif_K_Turb_G (DIMENSION_3p) )
        Allocate(  E_Turb_G_c   (DIMENSION_3p) )
        Allocate(  E_Turb_G_p   (DIMENSION_3p) )
        Allocate(  Dif_E_Turb_G (DIMENSION_3p) )
      ENDIF

! Simonin or Ahmadi model
      IF(KT_TYPE_ENUM==SIMONIN_1996 .OR.&
         KT_TYPE_ENUM==AHMADI_1995) THEN
        Allocate(  K_12 (DIMENSION_3) )
        Allocate(  Tau_12 (DIMENSION_3) )
        Allocate(  Tau_1 (DIMENSION_3) )
      ENDIF

!tau_g
      Allocate(  TAU_U_g(DIMENSION_3p) )
      Allocate(  TAU_V_g(DIMENSION_3p) )
      Allocate(  TAU_W_g(DIMENSION_3p) )
      Allocate(  DF_gu(DIMENSION_3p, -3:3) )
      Allocate(  DF_gv(DIMENSION_3p, -3:3) )
      Allocate(  DF_gw(DIMENSION_3p, -3:3) )
      Allocate(  CTAU_U_G(DIMENSION_3P))
      Allocate(  CTAU_V_G(DIMENSION_3P))
      Allocate(  CTAU_W_G(DIMENSION_3P))

!tau_s
      Allocate(  TAU_U_s(DIMENSION_3p, DIMENSION_M) )
      Allocate(  TAU_V_s(DIMENSION_3p, DIMENSION_M) )
      Allocate(  TAU_W_s(DIMENSION_3p, DIMENSION_M) )

! generate_particles / particle_count
      Allocate(  PARTICLE_COUNT(DIMENSION_3) )

!trace
      Allocate(  trD_s_C (DIMENSION_3, DIMENSION_M) )
      Allocate(  trD_s2 (DIMENSION_3, DIMENSION_M) )
      Allocate(  trD_s_Co (DIMENSION_3, DIMENSION_M) )
      Allocate(  trD_s_Co2 (DIMENSION_3, DIMENSION_M) )
!visc_g
      Allocate(  trD_g(DIMENSION_3) )
      Allocate(  MU_gt (DIMENSION_3) )
      Allocate(  EPMU_gt (DIMENSION_3p) )
      Allocate(  LAMBDA_gt (DIMENSION_3p) )
      Allocate(  EPLAMBDA_gt (DIMENSION_3) )
      Allocate(  L_scale (DIMENSION_3) )

!visc_s
      Allocate(  MU_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  EPMU_s (DIMENSION_3p, DIMENSION_M) )
      Allocate(  LAMBDA_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  EPLAMBDA_s (DIMENSION_3p, DIMENSION_M) )
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
      Allocate(  EP_star_array (DIMENSION_3) )
      Allocate(  EP_g_blend_start (DIMENSION_3) )
      Allocate(  EP_g_blend_end (DIMENSION_3) )
      Allocate(  trD_s(DIMENSION_3, DIMENSION_M) )
      Allocate(  I2_devD_s (DIMENSION_3) )
      Allocate(  TrM_s (DIMENSION_3) )
      Allocate(  TrDM_s (DIMENSION_3) )

!shear quantities
      Allocate(  VSH(DIMENSION_3) )
      Allocate(  VSHE(DIMENSION_3) )

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
      IF (KT_TYPE_ENUM == GHD_2007) THEN
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
      IF (KT_TYPE_ENUM == IA_2005) THEN
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
      IF (KT_TYPE_ENUM == GTSH_2012) THEN
         Allocate(  A2_gtsh(DIMENSION_3) )
         Allocate(  xsi_gtsh(DIMENSION_3) )
      ENDIF
      IF (KT_TYPE_ENUM == IA_2005 .OR. &
          KT_TYPE_ENUM == GD_1999 .OR. &
          KT_TYPE_ENUM == GTSH_2012) THEN
         Allocate(  EDT_s_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  EDvel_sM_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
      ENDIF

      Allocate(errorpercent(0:MMAX))

      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_ARRAYS_GEOMETRY                               !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Calculate X, X_E,  oX, oX_E                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_ARRAYS_GEOMETRY

! Global Variables:
!---------------------------------------------------------------------//
! Domain decomposition and dimensions
      use geometry, only: oDX, oDX_E
      use geometry, only: oDZ, oDZ_T
      use geometry, only: oDY, oDY_N
      use geometry, only: X, X_E, oX, oX_E, cyl_X, cyl_X_E
      use geometry, only: Z, Z_T
! Averaging factors.
      use geometry, only: FX_E, FX_E_bar, FX, FX_bar
      use geometry, only: FY_N, FY_N_bar
      use geometry, only: FZ_T, FZ_T_bar
! Domain flags.
      use geometry, only: ICBC_FLAG
      use geometry, only: FLAG, FLAG3
      use geometry, only: FLAG_E, FLAG_N, FLAG_T
! Domain volumes and areas.
      use geometry, only: VOL, VOL_SURR, AYZ, AXZ, AXY! Scalar grid
      use geometry, only: VOL_U, AYZ_U, AXZ_U, AXY_U  ! X-Momentum
      use geometry, only: VOL_V, AYZ_V, AXZ_V, AXY_V  ! Y-Momentum
      use geometry, only: VOL_W, AYZ_W, AXZ_W, AXY_W  ! Z-Momentum
! Axis decomposition
      USE param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
      USE param, only: DIMENSION_3, DIMENSION_4
      USE param, only: DIMENSION_3L, DIMENSION_3P
! Flag for POST_MFIX
      use cdist, only: bDoing_postmfix

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager
      use compar, only:ADJUST_PARTITION

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Error Flag
      INTEGER :: IER
! Flag indicating that the arrays were previously allocated.
      INTEGER, SAVE :: CALLED = -1
!......................................................................!

      IF(ADJUST_PARTITION) CALLED = -1

      CALLED = CALLED + 1

      IF(CALLED > 0) THEN
         IF(.NOT.bDoing_postmfix) THEN
            RETURN
         ELSEIF(mod(CALLED,2) /= 0) THEN
            RETURN
         ENDIF
      ENDIF

! Initialize the error manager.
      CALL INIT_ERR_MSG("ALLOCATE_ARRAYS_GEOMETRY")

! Allocate geometry components related to the mesh. Check the
! allocation error status and abort if any failure is detected.
      ALLOCATE( X     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( cyl_X     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( X_E   (0:DIMENSION_I), STAT=IER)
      ALLOCATE( cyl_X_E   (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oX    (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oX_E  (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oDX   (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oDX_E (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( oDY   (0:DIMENSION_J), STAT=IER )
      ALLOCATE( oDY_N (0:DIMENSION_J), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( Z     (0:DIMENSION_K), STAT=IER )
      ALLOCATE( Z_T   (0:DIMENSION_K), STAT=IER )
      ALLOCATE( oDZ   (0:DIMENSION_K), STAT=IER )
      ALLOCATE( oDZ_T (0:DIMENSION_K), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( FX     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( FX_bar (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( FX_E     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( FX_E_bar (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( FY_N     (0:DIMENSION_J), STAT=IER )
      ALLOCATE( FY_N_bar (0:DIMENSION_J), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( FZ_T     (0:DIMENSION_K), STAT=IER )
      ALLOCATE( FZ_T_bar (0:DIMENSION_K), STAT=IER )
      IF(IER /= 0) goto 500

! Flags for the scalar grid.
      Allocate( FLAG  (DIMENSION_3), STAT=IER )
      Allocate( FLAG3 (DIMENSION_4), STAT=IER )
      IF(IER /= 0) goto 500

! Flags for the momentum grids.
      Allocate( FLAG_E (DIMENSION_3), STAT=IER )
      Allocate( FLAG_N (DIMENSION_3), STAT=IER )
      Allocate( FLAG_T (DIMENSION_3), STAT=IER )
      IF(IER /= 0) goto 500

! Text flags for scalar grid.
      Allocate( ICBC_FLAG (DIMENSION_3L), STAT=IER )
      IF(IER /= 0) goto 500

! Volume and face-areas of scalar grid.
      Allocate( VOL (DIMENSION_3),  STAT=IER )
      Allocate( AYZ (DIMENSION_3P), STAT=IER )
      Allocate( AXZ (DIMENSION_3P), STAT=IER )
      Allocate( AXY (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

      ! total volume of each cell's surrounding stencil cells
      Allocate( VOL_SURR (DIMENSION_3), STAT=IER )

! Volume and face-areas of X-Momentumn grid.
      Allocate( VOL_U (DIMENSION_3),  STAT=IER )
      Allocate( AYZ_U (DIMENSION_3P), STAT=IER )
      Allocate( AXZ_U (DIMENSION_3P), STAT=IER )
      Allocate( AXY_U (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

! Volume and face-areas of Y-Momentum grid.
      Allocate( VOL_V (DIMENSION_3),  STAT=IER )
      Allocate( AYZ_V (DIMENSION_3P), STAT=IER )
      Allocate( AXZ_V (DIMENSION_3P), STAT=IER )
      Allocate( AXY_V (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

! Volume and face-areas of Z-Momentum grid.
      Allocate( VOL_W (DIMENSION_3),  STAT=IER )
      Allocate( AYZ_W (DIMENSION_3P), STAT=IER )
      Allocate( AXZ_W (DIMENSION_3P), STAT=IER )
      Allocate( AXY_W (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

! Collect the error flags from all ranks. If all allocaitons were
! successfull, do nothing. Otherwise, flag the error and abort.
! Note that the allocation status is checked in groups. This can
! be increase if tracking the source of an allocation failure.
  500 CALL GLOBAL_ALL_SUM(IER)

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Failure during array allocation.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_ARRAYS_INCREMENTS                             !
!  Author: M. Syamlal, W. Rogers                      Date: 10-DEC-91  !
!                                                                      !
!  Purpose: The purpose of this module is to create increments to be   !
!           stored in the array STORE_INCREMENT which will be added    !
!           to cell index ijk to find the effective indices of its     !
!           neighbors. These increments are found using the 'class'    !
!           of cell ijk. The class is determined based on the          !
!           neighboring cell type, i.e. wall or fluid.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_ARRAYS_INCREMENTS

      USE param
      USE param1
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE fldvar
      USE funits

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Error flag.
      INTEGER :: IER
! Flag indicating that the arrays were previously allocated.
      LOGICAL, SAVE :: ALREADY_ALLOCATED = .FALSE.
!......................................................................!

      IF(ADJUST_PARTITION) THEN
         CALL DEALLOCATE_ARRAYS_INCREMENTS
         ALREADY_ALLOCATED = .FALSE.
         INCREMENT_ARRAYS_ALLOCATED = .FALSE.
      ENDIF

      IF(ALREADY_ALLOCATED) RETURN

! Initialize the error manager.
      CALL INIT_ERR_MSG("ALLOCATE_ARRAYS_INCREMENTS")

! Allocate increment arrays and report an allocation errors.
      Allocate( I_OF (DIMENSION_3), STAT=IER)
      Allocate( J_OF (DIMENSION_3), STAT=IER)
      Allocate( K_OF (DIMENSION_3), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Im1 (0:DIMENSION_I), STAT=IER)
      Allocate( Ip1 (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Jm1 (0:DIMENSION_J), STAT=IER)
      Allocate( Jp1 (0:DIMENSION_J), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Km1 (0:DIMENSION_K), STAT=IER)
      Allocate( Kp1 (0:DIMENSION_K), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( STORE_LM (DIMENSION_M, DIMENSION_M), STAT=IER)
      Allocate( CELL_CLASS (DIMENSION_3), STAT=IER)
      IF(IER /= 0) goto 500


! Allocate increment arrays and report an allocation errors.
      Allocate( I3_OF (DIMENSION_4), STAT=IER)
      Allocate( J3_OF (DIMENSION_4), STAT=IER)
      Allocate( K3_OF (DIMENSION_4), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Im1_3 (-1:DIMENSION_I+1), STAT=IER)
      Allocate( Ip1_3 (-1:DIMENSION_I+1), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Jm1_3 (-1:DIMENSION_J+1), STAT=IER)
      Allocate( Jp1_3 (-1:DIMENSION_J+1), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Km1_3 (-1:DIMENSION_K+1), STAT=IER)
      Allocate( Kp1_3 (-1:DIMENSION_K+1), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( CELL_CLASS3 (DIMENSION_4), STAT=IER)
      IF(IER /= 0) goto 500

! Collect the error flags from all ranks. If all allocaitons were
! successfull, do nothing. Otherwise, flag the error and abort.
! Note that the allocation status is checked in groups. This can
! be increase if tracking the source of an allocation failure.
  500 CALL GLOBAL_ALL_SUM(IER)

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Failure during array allocation.')

      ALREADY_ALLOCATED = .TRUE.

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS_INCREMENTS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutinee: DEALLOCATE_ARRAYS_                                     C
!  Purpose: Deallocate arrays                                          C
!                                                                      C
!  Author: Jeff Dietiker                             Date: 03-MAR-2016 C
!  Reviewer:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DEALLOCATE_ARRAYS  !@@
      
      CALL DEALLOCATE_ARRAYS_MAIN
      CALL DEALLOCATE_ARRAYS_GEOMETRY      
      CALL DEALLOCATE_ARRAYS_INCREMENTS
      CALL DEALLOCATE_ARRAYS_PARALLEL
      CALL DEALLOCATE_CUT_CELL_ARRAYS
      CALL DEALLOCATE_DEM_MI
      CALL DEALLOCATE_PIC_MIO
      CALL DES_DEALLOCATE_ARRAYS
      
      
      END SUBROUTINE DEALLOCATE_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutinee: DEALLOCATE_ARRAYS_MAIN                                 C
!  Purpose: Deallocate arrays                                          C
!                                                                      C
!  Author: Jeff Dietiker                             Date: 03-MAR-2016 C
!  Reviewer:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DEALLOCATE_ARRAYS_MAIN

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use ambm
      use cdist
      use cont
      use des_rxns
      use drag
      use energy
      use fldvar
      use generate_particles, only: particle_count
      use geometry
      use ghdtheory
      use indices
      use kintheory
      use mflux
      use param
      use param1
      use pgcor
      use physprop
      use pscor
      use residual
      use run
      use rxns
      use scalars
      use iterate, only: errorpercent
      use tau_g
      use tau_s
      use trace
      use turb
      use visc_g
      use visc_s
      use vshear

      IMPLICIT NONE

!-----------------------------------------------
! Variables
!-----------------------------------------------

!ambm
      if(allocated(  A_m )) deallocate(  A_m )
      if(allocated(  B_m )) deallocate(  B_m )

!cont
      if(allocated(  DO_CONT )) deallocate(  DO_CONT )

!drag
      if(allocated(   F_gs )) deallocate(   F_gs )
      if(allocated(   F_ss )) deallocate(   F_ss )


!Off diagonal friction coefficient in HYS drag relation
      if(allocated(   beta_ij )) deallocate(   beta_ij )


!energy
      if(allocated(   HOR_g  )) deallocate(   HOR_g  )
      if(allocated(   HOR_s  )) deallocate(   HOR_s  )
      if(allocated(   GAMA_gs  )) deallocate(   GAMA_gs  )
      if(allocated(   GAMA_Rg  )) deallocate(   GAMA_Rg  )
      if(allocated(   GAMA_Rs  )) deallocate(   GAMA_Rs  )
      if(allocated(   T_Rg  )) deallocate(   T_Rg  )
      if(allocated(   T_Rs  )) deallocate(   T_Rs  )


!fldvar

      if(allocated(   EP_g  )) deallocate(   EP_g  )
      if(allocated(   epg_jfac  )) deallocate(   epg_jfac  )
      if(allocated(   epg_ifac  )) deallocate(   epg_ifac  )
      if(allocated(   eps_ifac  )) deallocate(   eps_ifac  )
      if(allocated(   EP_go  )) deallocate(   EP_go  )
      if(allocated(   P_g  )) deallocate(   P_g  )
      if(allocated( P_go )) deallocate( P_go )
      if(allocated( RO_g )) deallocate( RO_g )
      if(allocated( RO_go )) deallocate( RO_go )
      if(allocated( ROP_g )) deallocate( ROP_g )
      if(allocated( ROP_go )) deallocate( ROP_go )
      if(allocated( RO_S  )) deallocate( RO_S  )
      if(allocated( RO_So  )) deallocate( RO_So  )
      if(allocated( ROP_s  )) deallocate( ROP_s  )
      if(allocated( ROP_so  )) deallocate( ROP_so  )

      if(allocated( EP_SS )) deallocate( EP_SS )
      if(allocated( ERR_ARRAY )) deallocate( ERR_ARRAY )

      if(allocated( T_g  )) deallocate( T_g  )
      if(allocated( T_s  )) deallocate( T_s  )
      if(allocated( T_go  )) deallocate( T_go  )
      if(allocated( T_so  )) deallocate( T_so  )
      if(allocated( X_g  )) deallocate( X_g  )
      if(allocated( X_s  )) deallocate( X_s  )
      if(allocated( X_go  )) deallocate( X_go  )
      if(allocated( X_so  )) deallocate( X_so  )
      if(allocated( U_g  )) deallocate( U_g  )
      if(allocated( U_go  )) deallocate( U_go  )
      if(allocated( U_s  )) deallocate( U_s  )
      if(allocated( U_so  )) deallocate( U_so  )
      if(allocated( V_g  )) deallocate( V_g  )
      if(allocated( V_go  )) deallocate( V_go  )
      if(allocated( V_s  )) deallocate( V_s  )
      if(allocated( V_so  )) deallocate( V_so  )
      if(allocated( W_g  )) deallocate( W_g  )
      if(allocated( W_go  )) deallocate( W_go  )
      if(allocated( W_s  )) deallocate( W_s  )
      if(allocated( W_so  )) deallocate( W_so  )
      if(allocated( P_s  )) deallocate( P_s  )
      if(allocated( P_s_c  )) deallocate( P_s_c  )
      if(allocated( P_s_v  )) deallocate( P_s_v  )
      if(allocated( P_s_f  )) deallocate( P_s_f  )
      if(allocated( P_s_p  )) deallocate( P_s_p  )
      if(allocated( P_star  )) deallocate( P_star  )
      if(allocated( P_staro  )) deallocate( P_staro  )
      if(allocated( THETA_m  )) deallocate( THETA_m  )
      if(allocated( THETA_mo  )) deallocate( THETA_mo  )


      IF(K_Epsilon)THEN
        if(allocated( K_Turb_G  )) deallocate( K_Turb_G  )
        if(allocated( K_Turb_Go  )) deallocate( K_Turb_Go  )
        if(allocated( E_Turb_G  )) deallocate( E_Turb_G  )
        if(allocated( E_Turb_Go  )) deallocate( E_Turb_Go  )
      ENDIF

      if(allocated( Scalar  )) deallocate( Scalar  )
      if(allocated( Scalaro  )) deallocate( Scalaro  )



!pgcor
      if(allocated( d_e )) deallocate( d_e )
      if(allocated( d_n )) deallocate( d_n )
      if(allocated( d_t )) deallocate( d_t )
      if(allocated( Pp_g )) deallocate( Pp_g )
      if(allocated( PHASE_4_P_g )) deallocate( PHASE_4_P_g )

!physprop
      if(allocated( MU_g  )) deallocate( MU_g  )
      if(allocated( C_pg  )) deallocate( C_pg  )
      if(allocated( C_ps  )) deallocate( C_ps  )
      if(allocated( K_g  )) deallocate( K_g  )
      if(allocated( K_s  )) deallocate( K_s  )
      if(allocated( Kth_s  )) deallocate( Kth_s  )
      if(allocated( Kphi_s  )) deallocate( Kphi_s  )
      if(allocated( DIF_g  )) deallocate( DIF_g  )
      if(allocated( DIF_s  )) deallocate( DIF_s  )
      if(allocated( MW_MIX_g  )) deallocate( MW_MIX_g  )

!pscor
      if(allocated( e_e )) deallocate( e_e )
      if(allocated( e_n )) deallocate( e_n )
      if(allocated( e_t )) deallocate( e_t )
      if(allocated( K_cp )) deallocate( K_cp )
      if(allocated( EPp )) deallocate( EPp )
      if(allocated( PHASE_4_P_s )) deallocate( PHASE_4_P_s )


!residual
      if(allocated(  RESID )) deallocate(  RESID )
      if(allocated(  MAX_RESID )) deallocate(  MAX_RESID )
      if(allocated(  IJK_RESID )) deallocate(  IJK_RESID )
      if(allocated(  NUM_RESID )) deallocate(  NUM_RESID )
      if(allocated(  DEN_RESID )) deallocate(  DEN_RESID )
      if(allocated(  RESID_PACK )) deallocate(  RESID_PACK )

!rxns
      if(allocated( R_gp  )) deallocate( R_gp  )
      if(allocated( R_sp  )) deallocate( R_sp  )
      if(allocated( RoX_gc  )) deallocate( RoX_gc  )
      if(allocated( RoX_sc  )) deallocate( RoX_sc  )
      if(allocated( SUM_R_g  )) deallocate( SUM_R_g  )
      if(allocated( SUM_R_s  )) deallocate( SUM_R_s  )
      if(allocated( R_phase  )) deallocate( R_phase  )

!scalars
      if(allocated( Scalar_c  )) deallocate( Scalar_c  )
      if(allocated( Scalar_p  )) deallocate( Scalar_p  )
      if(allocated( Dif_Scalar  )) deallocate( Dif_Scalar  )

! add by rong for dqmom
      if(allocated( D_p   )) deallocate( D_p   )
      if(allocated( D_po  )) deallocate( D_po  )
      if(allocated( Source_a )) deallocate( Source_a )
      if(allocated( S_bar )) deallocate( S_bar )
      if(allocated( Matrix_a )) deallocate( Matrix_a )
      if(allocated( Matrix_b )) deallocate( Matrix_b )
      if(allocated( Matrix_c )) deallocate( Matrix_c )
      if(allocated( Inv_a )) deallocate( Inv_a )
      if(allocated( A )) deallocate( A )
      if(allocated( omega )) deallocate( omega )
      if(allocated( beta_a )) deallocate( beta_a )
      if(allocated( ystart )) deallocate( ystart )

! K-Epsilon Turbulence model
      if(allocated( K_Turb_G_c    )) deallocate( K_Turb_G_c    )
      if(allocated( K_Turb_G_p    )) deallocate( K_Turb_G_p    )
      if(allocated( Dif_K_Turb_G  )) deallocate( Dif_K_Turb_G  )
      if(allocated( E_Turb_G_c    )) deallocate( E_Turb_G_c    )
      if(allocated( E_Turb_G_p    )) deallocate( E_Turb_G_p    )
      if(allocated( Dif_E_Turb_G  )) deallocate( Dif_E_Turb_G  )

! Simonin or Ahmadi model
      if(allocated( K_12  )) deallocate( K_12  )
      if(allocated( Tau_12  )) deallocate( Tau_12  )
      if(allocated( Tau_1  )) deallocate( Tau_1  )


!tau_g
      if(allocated( TAU_U_g )) deallocate( TAU_U_g )
      if(allocated( TAU_V_g )) deallocate( TAU_V_g )
      if(allocated( TAU_W_g )) deallocate( TAU_W_g )
      if(allocated( DF_gu )) deallocate( DF_gu )
      if(allocated( DF_gv )) deallocate( DF_gv )
      if(allocated( DF_gw )) deallocate( DF_gw )
      if(allocated( CTAU_U_G )) deallocate( CTAU_U_G )
      if(allocated( CTAU_V_G )) deallocate( CTAU_V_G )
      if(allocated( CTAU_W_G )) deallocate( CTAU_W_G )

!tau_s
      if(allocated( TAU_U_s )) deallocate( TAU_U_s )
      if(allocated( TAU_V_s )) deallocate( TAU_V_s )
      if(allocated( TAU_W_s )) deallocate( TAU_W_s )


! generate_particles / particle_count
      if(allocated( PARTICLE_COUNT )) deallocate( PARTICLE_COUNT )


!trace
      if(allocated( trD_s_C  )) deallocate( trD_s_C  )
      if(allocated( trD_s2  )) deallocate( trD_s2  )
      if(allocated( trD_s_Co  )) deallocate( trD_s_Co  )
      if(allocated( trD_s_Co2  )) deallocate( trD_s_Co2  )

!visc_g
      if(allocated( trD_g )) deallocate( trD_g )
      if(allocated( MU_gt  )) deallocate( MU_gt  )
      if(allocated( EPMU_gt  )) deallocate( EPMU_gt  )
      if(allocated( LAMBDA_gt  )) deallocate( LAMBDA_gt  )
      if(allocated( EPLAMBDA_gt  )) deallocate( EPLAMBDA_gt  )
      if(allocated( L_scale  )) deallocate( L_scale  )


!visc_s
      if(allocated( MU_s  )) deallocate( MU_s  )
      if(allocated( EPMU_s  )) deallocate( EPMU_s  )
      if(allocated( LAMBDA_s  )) deallocate( LAMBDA_s  )
      if(allocated( EPLAMBDA_s  )) deallocate( EPLAMBDA_s  )
      if(allocated( ALPHA_s  )) deallocate( ALPHA_s  )
      if(allocated( MU_s_c  )) deallocate( MU_s_c  )
      if(allocated( LAMBDA_s_c  )) deallocate( LAMBDA_s_c  )
      if(allocated( LAMBDA_s_v  )) deallocate( LAMBDA_s_v  )
      if(allocated( LAMBDA_s_f  )) deallocate( LAMBDA_s_f  )
      if(allocated( LAMBDA_s_p  )) deallocate( LAMBDA_s_p  )
      if(allocated( MU_s_v  )) deallocate( MU_s_v  )
      if(allocated( MU_s_f  )) deallocate( MU_s_f  )
      if(allocated( MU_s_p  )) deallocate( MU_s_p  )
      if(allocated( MU_b_v  )) deallocate( MU_b_v  )
      if(allocated( EP_star_array  )) deallocate( EP_star_array  )
      if(allocated( EP_g_blend_start  )) deallocate( EP_g_blend_start  )
      if(allocated( EP_g_blend_end  )) deallocate( EP_g_blend_end  )
      if(allocated( trD_s )) deallocate( trD_s )
      if(allocated( I2_devD_s  )) deallocate( I2_devD_s  )
      if(allocated( TrM_s  )) deallocate( TrM_s  )
      if(allocated( TrDM_s  )) deallocate( TrDM_s  )


!shear quantities
      if(allocated( VSH )) deallocate( VSH )
      if(allocated( VSHE )) deallocate( VSHE )

!mflux
      if(allocated(  Flux_gE )) deallocate(  Flux_gE )
      if(allocated(  Flux_sE )) deallocate(  Flux_sE )
      if(allocated(  Flux_gN )) deallocate(  Flux_gN )
      if(allocated(  Flux_sN )) deallocate(  Flux_sN )
      if(allocated(  Flux_gT )) deallocate(  Flux_gT )
      if(allocated(  Flux_sT )) deallocate(  Flux_sT )

      if(allocated(  Flux_gSE )) deallocate(  Flux_gSE )
      if(allocated(  Flux_sSE )) deallocate(  Flux_sSE )
      if(allocated(  Flux_gSN )) deallocate(  Flux_gSN )
      if(allocated(  Flux_sSN )) deallocate(  Flux_sSN )
      if(allocated(  Flux_gST )) deallocate(  Flux_gST )
      if(allocated(  Flux_sST )) deallocate(  Flux_sST )

      if(allocated(  ROP_gE )) deallocate(  ROP_gE )
      if(allocated(  ROP_sE )) deallocate(  ROP_sE )
      if(allocated(  ROP_gN )) deallocate(  ROP_gN )
      if(allocated(  ROP_sN )) deallocate(  ROP_sN )
      if(allocated(  ROP_gT )) deallocate(  ROP_gT )
      if(allocated(  ROP_sT )) deallocate(  ROP_sT )

! allocate variables for GHD Theory
      if(allocated( Flux_nE )) deallocate( Flux_nE )
      if(allocated( Flux_nN )) deallocate( Flux_nN )
      if(allocated( Flux_nT )) deallocate( Flux_nT )
      if(allocated( Zeta0 )) deallocate( Zeta0 )
      if(allocated( ZetaU )) deallocate( ZetaU )
      if(allocated( DiT )) deallocate( DiT )
      if(allocated( DijF )) deallocate( DijF )
      if(allocated( Lij )) deallocate( Lij )
      if(allocated( Dij )) deallocate( Dij )
      if(allocated( DijQ )) deallocate( DijQ )
      if(allocated( JoiX )) deallocate( JoiX )
      if(allocated( JoiY )) deallocate( JoiY )
      if(allocated( JoiZ )) deallocate( JoiZ )
      if(allocated( FiX )) deallocate( FiX )
      if(allocated( FiY )) deallocate( FiY )
      if(allocated( FiZ )) deallocate( FiZ )
      if(allocated( FiXvel )) deallocate( FiXvel )
      if(allocated( FiYvel )) deallocate( FiYvel )
      if(allocated( FiZvel )) deallocate( FiZvel )
      if(allocated( DELTAU )) deallocate( DELTAU )
      if(allocated( DELTAV )) deallocate( DELTAV )
      if(allocated( DELTAW )) deallocate( DELTAW )
      if(allocated( dragFx )) deallocate( dragFx )
      if(allocated( dragFy )) deallocate( dragFy )
      if(allocated( dragFz )) deallocate( dragFz )
      if(allocated( dragFxflux )) deallocate( dragFxflux )
      if(allocated( dragFyflux )) deallocate( dragFyflux )
      if(allocated( dragFzflux )) deallocate( dragFzflux )
      if(allocated( FiMinusDragX )) deallocate( FiMinusDragX )
      if(allocated( JoiMinusDragX )) deallocate( JoiMinusDragX )
      if(allocated( FiMinusDragY )) deallocate( FiMinusDragY )
      if(allocated( JoiMinusDragY )) deallocate( JoiMinusDragY )
      if(allocated( FiMinusDragZ )) deallocate( FiMinusDragZ )
      if(allocated( JoiMinusDragZ )) deallocate( JoiMinusDragZ )
      if(allocated( beta_cell_X )) deallocate( beta_cell_X )
      if(allocated( beta_cell_Y )) deallocate( beta_cell_Y )
      if(allocated( beta_cell_Z )) deallocate( beta_cell_Z )
      if(allocated( beta_ij_cell_X )) deallocate( beta_ij_cell_X )
      if(allocated( beta_ij_cell_Y )) deallocate( beta_ij_cell_Y )
      if(allocated( beta_ij_cell_Z )) deallocate( beta_ij_cell_Z )
      if(allocated( DEL_DOT_J )) deallocate( DEL_DOT_J )
      if(allocated( DiT_HarmE )) deallocate( DiT_HarmE )
      if(allocated( DiT_HarmN )) deallocate( DiT_HarmN )
      if(allocated( DiT_HarmT )) deallocate( DiT_HarmT )
      if(allocated( Dij_HarmE )) deallocate( Dij_HarmE )
      if(allocated( Dij_HarmN )) deallocate( Dij_HarmN )
      if(allocated( Dij_HarmT )) deallocate( Dij_HarmT )
      if(allocated( DijF_HarmE )) deallocate( DijF_HarmE )
      if(allocated( DijF_HarmN )) deallocate( DijF_HarmN )
      if(allocated( DijF_HarmT )) deallocate( DijF_HarmT )


! We need to set this even when KT_TYPE is not set to IA_NONEP - at
! least in the current version of the code and needs to be revisited
      if(allocated( KTMOM_U_s )) deallocate( KTMOM_U_s )
      if(allocated( KTMOM_V_s )) deallocate( KTMOM_V_s )
      if(allocated( KTMOM_W_s )) deallocate( KTMOM_W_s )

! allocate variables for Iddir & Arastoopour (2005) kinetic theory
! EDvel_sM_ip & EDT_s_ip are also used for Garzy & Dufty (1999) kinetic theory
      if(allocated( trD_s2_ip )) deallocate( trD_s2_ip )
      if(allocated( MU_sM_ip )) deallocate( MU_sM_ip )
      if(allocated( MU_sL_ip )) deallocate( MU_sL_ip )
      if(allocated( XI_sM_ip )) deallocate( XI_sM_ip )
      if(allocated( XI_sL_ip )) deallocate( XI_sL_ip )
      if(allocated( Fnu_s_ip )) deallocate( Fnu_s_ip )
      if(allocated( FT_sM_ip )) deallocate( FT_sM_ip )
      if(allocated( FT_sL_ip )) deallocate( FT_sL_ip )
      if(allocated( Kth_sL_ip )) deallocate( Kth_sL_ip )
      if(allocated( Knu_sM_ip )) deallocate( Knu_sM_ip )
      if(allocated( Knu_sL_ip )) deallocate( Knu_sL_ip )
      if(allocated( Kvel_s_ip )) deallocate( Kvel_s_ip )
      if(allocated( EDvel_sL_ip )) deallocate( EDvel_sL_ip )
      if(allocated( ED_ss_ip )) deallocate( ED_ss_ip )

      if(allocated( A2_gtsh )) deallocate( A2_gtsh )
      if(allocated( xsi_gtsh )) deallocate( xsi_gtsh )

      if(allocated( EDT_s_ip )) deallocate( EDT_s_ip )
      if(allocated( EDvel_sM_ip )) deallocate( EDvel_sM_ip )

      if(allocated( errorpercent )) deallocate( errorpercent )

      RETURN
      END SUBROUTINE DEALLOCATE_ARRAYS_MAIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DEALLOCATE_ARRAYS_GEOMETRY                             !
!  Author: Jeff DIetiker                              Date: 16-MAR-2016!
!                                                                      !
!  Purpose: Calculate X, X_E,  oX, oX_E                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEALLOCATE_ARRAYS_GEOMETRY

! Global Variables:
!---------------------------------------------------------------------//
! Domain decomposition and dimensions
      use geometry, only: oDX, oDX_E
      use geometry, only: oDZ, oDZ_T
      use geometry, only: oDY, oDY_N
      use geometry, only: X, X_E, oX, oX_E, cyl_X, cyl_X_E
      use geometry, only: Z, Z_T
! Averaging factors.
      use geometry, only: FX_E, FX_E_bar, FX, FX_bar
      use geometry, only: FY_N, FY_N_bar
      use geometry, only: FZ_T, FZ_T_bar
! Domain flags.
      use geometry, only: ICBC_FLAG
      use geometry, only: FLAG, FLAG3
      use geometry, only: FLAG_E, FLAG_N, FLAG_T
! Domain volumes and areas.
      use geometry, only: VOL, VOL_SURR, AYZ, AXZ, AXY! Scalar grid
      use geometry, only: VOL_U, AYZ_U, AXZ_U, AXY_U  ! X-Momentum
      use geometry, only: VOL_V, AYZ_V, AXZ_V, AXY_V  ! Y-Momentum
      use geometry, only: VOL_W, AYZ_W, AXZ_W, AXY_W  ! Z-Momentum
! Axis decomposition
      USE param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
      USE param, only: DIMENSION_3, DIMENSION_4
      USE param, only: DIMENSION_3L, DIMENSION_3P
! Flag for POST_MFIX
      use cdist, only: bDoing_postmfix

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Allocate geometry components related to the mesh. Check the
! allocation error status and abort if any failure is detected.
      if(allocated(  X      )) deallocate(  X      )
      if(allocated(  cyl_X      )) deallocate(  cyl_X      )
      if(allocated(  X_E    )) deallocate(  X_E    )
      if(allocated(  cyl_X_E    )) deallocate(  cyl_X_E    )
      if(allocated(  oX     )) deallocate(  oX     )
      if(allocated(  oX_E   )) deallocate(  oX_E   )
      if(allocated(  oDX    )) deallocate(  oDX    )
      if(allocated(  oDX_E  )) deallocate(  oDX_E  )
      
     
      if(allocated(  oDY    )) deallocate(  oDY    )
      if(allocated(  oDY_N  )) deallocate(  oDY_N  )
    
   
      if(allocated(  Z      )) deallocate(  Z      )
      if(allocated(  Z_T    )) deallocate(  Z_T    )
      if(allocated(  oDZ    )) deallocate(  oDZ    )
      if(allocated(  oDZ_T  )) deallocate(  oDZ_T  )
  
 
      if(allocated(  FX      )) deallocate(  FX      )
      if(allocated(  FX_bar  )) deallocate(  FX_bar  )
      
     
      if(allocated(  FX_E      )) deallocate(  FX_E      )
      if(allocated(  FX_E_bar  )) deallocate(  FX_E_bar  )
    
   
      if(allocated(  FY_N      )) deallocate(  FY_N      )
      if(allocated(  FY_N_bar  )) deallocate(  FY_N_bar  )
 
      if(allocated(  FZ_T      )) deallocate(  FZ_T      )
      if(allocated(  FZ_T_bar  )) deallocate(  FZ_T_bar  )


! Flags for the scalar grid.      
      if(allocated(  FLAG   )) deallocate(  FLAG   )
      if(allocated(  FLAG3  )) deallocate(  FLAG3  )

! Flags for the momentum grids.      
      if(allocated(  FLAG_E  )) deallocate(  FLAG_E  )
      if(allocated(  FLAG_N  )) deallocate(  FLAG_N  )
      if(allocated(  FLAG_T  )) deallocate(  FLAG_T  )

! Text flags for scalar grid.
!      if(allocated(  ICBC_FLAG  )) deallocate(  ICBC_FLAG  )

! Volume and face-areas of scalar grid.
      if(allocated(  VOL  )) deallocate(  VOL  )
      if(allocated(  AYZ  )) deallocate(  AYZ  )
      if(allocated(  AXZ  )) deallocate(  AXZ  )
      if(allocated(  AXY  )) deallocate(  AXY  )

! total volume of each cell's surrounding stencil cells
      if(allocated(  VOL_SURR  )) deallocate(  VOL_SURR  )

! Volume and face-areas of X-Momentumn grid.
      if(allocated(  VOL_U  )) deallocate(  VOL_U  )
      if(allocated(  AYZ_U  )) deallocate(  AYZ_U  )
      if(allocated(  AXZ_U  )) deallocate(  AXZ_U  )
      if(allocated(  AXY_U  )) deallocate(  AXY_U  )

! Volume and face-areas of Y-Momentum grid.
      if(allocated(  VOL_V  )) deallocate(  VOL_V  )
      if(allocated(  AYZ_V  )) deallocate(  AYZ_V  )
      if(allocated(  AXZ_V  )) deallocate(  AXZ_V  )
      if(allocated(  AXY_V  )) deallocate(  AXY_V  )

! Volume and face-areas of Z-Momentum grid.
      if(allocated(  VOL_W  )) deallocate(  VOL_W  )
      if(allocated(  AYZ_W  )) deallocate(  AYZ_W  )
      if(allocated(  AXZ_W  )) deallocate(  AXZ_W  )
      if(allocated(  AXY_W  )) deallocate(  AXY_W  )


      RETURN
      END SUBROUTINE DEALLOCATE_ARRAYS_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DEALLOCATE_ARRAYS_INCREMENTS                           !
!  Author: Jeff Dietiker                              Date: 16-MAR-2016!
!                                                                      !
!  Purpose: The purpose of this module is to create increments to be   !
!           stored in the array STORE_INCREMENT which will be added    !
!           to cell index ijk to find the effective indices of its     !
!           neighbors. These increments are found using the 'class'    !
!           of cell ijk. The class is determined based on the          !
!           neighboring cell type, i.e. wall or fluid.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEALLOCATE_ARRAYS_INCREMENTS

      USE param
      USE param1
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE fldvar
      USE funits

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Error flag.
      INTEGER :: IER


! Allocate increment arrays and report an allocation errors.
      if(allocated(  I_OF  )) deallocate(  I_OF  )
      if(allocated(  J_OF  )) deallocate(  J_OF  )
      if(allocated(  K_OF  )) deallocate(  K_OF  )

      if(allocated(  Im1  )) deallocate(  Im1  )
      if(allocated(  Ip1  )) deallocate(  Ip1  )

      if(allocated(  Jm1  )) deallocate(  Jm1  )
      if(allocated(  Jp1  )) deallocate(  Jp1  )

      if(allocated(  Km1  )) deallocate(  Km1  )
      if(allocated(  Kp1  )) deallocate(  Kp1  )

      if(allocated(  STORE_LM  )) deallocate(  STORE_LM  )
      if(allocated(  CELL_CLASS  )) deallocate(  CELL_CLASS  )

! Allocate increment arrays and report an allocation errors.
      if(allocated(  I3_OF  )) deallocate(  I3_OF  )
      if(allocated(  J3_OF  )) deallocate(  J3_OF  )
      if(allocated(  K3_OF  )) deallocate(  K3_OF  )

      if(allocated(  Im1_3  )) deallocate(  Im1_3  )
      if(allocated(  Ip1_3  )) deallocate(  Ip1_3  )

      if(allocated(  Jm1_3  )) deallocate(  Jm1_3  )
      if(allocated(  Jp1_3  )) deallocate(  Jp1_3  )

      if(allocated(  Km1_3  )) deallocate(  Km1_3  )
      if(allocated(  Kp1_3  )) deallocate(  Kp1_3  )

      if(allocated(  CELL_CLASS3  )) deallocate(  CELL_CLASS3  )

      if(allocated(WEST_ARRAY_OF)) deallocate(WEST_ARRAY_OF)
      if(allocated(EAST_ARRAY_OF)) deallocate(EAST_ARRAY_OF)
      if(allocated(SOUTH_ARRAY_OF)) deallocate(SOUTH_ARRAY_OF)
      if(allocated(NORTH_ARRAY_OF)) deallocate(NORTH_ARRAY_OF)
      if(allocated(BOTTOM_ARRAY_OF)) deallocate(BOTTOM_ARRAY_OF)
      if(allocated(TOP_ARRAY_OF)) deallocate(TOP_ARRAY_OF)

      if(allocated(IM_ARRAY_OF)) deallocate(IM_ARRAY_OF)
      if(allocated(IP_ARRAY_OF)) deallocate(IP_ARRAY_OF)
      if(allocated(JM_ARRAY_OF)) deallocate(JM_ARRAY_OF)
      if(allocated(JP_ARRAY_OF)) deallocate(JP_ARRAY_OF)
      if(allocated(KM_ARRAY_OF)) deallocate(KM_ARRAY_OF)
      if(allocated(KP_ARRAY_OF)) deallocate(KP_ARRAY_OF)


      RETURN
      END SUBROUTINE DEALLOCATE_ARRAYS_INCREMENTS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DEALLOCATE_ARRAYS_PARALLEL                             !
!  Author: Jeff Dietiker                              Date: 16-MAR-2016!
!                                                                      !
!  Purpose: The purpose of this module is to create increments to be   !
!           stored in the array STORE_INCREMENT which will be added    !
!           to cell index ijk to find the effective indices of its     !
!           neighbors. These increments are found using the 'class'    !
!           of cell ijk. The class is determined based on the          !
!           neighboring cell type, i.e. wall or fluid.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEALLOCATE_ARRAYS_PARALLEL

      USE param
      USE param1
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE fldvar
      USE funits

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager
      use desgrid
      use derived_types, only: dg_pic
      use desmpi


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Error flag.
      INTEGER :: IER


! gridmap_mod

      if(allocated(ISIZE_ALL)) deallocate(ISIZE_ALL)
      if(allocated(JSIZE_ALL)) deallocate(JSIZE_ALL)
      if(allocated(KSIZE_ALL)) deallocate(KSIZE_ALL)
 
      if(allocated(imap)) deallocate(imap)
      if(allocated(jmap)) deallocate(jmap)
      if(allocated(kmap)) deallocate(kmap)
 
      if(allocated(imap_c)) deallocate(imap_c)
      if(allocated(jmap_c)) deallocate(jmap_c)
      if(allocated(kmap_c)) deallocate(kmap_c)

! desgrid_mod
      if(allocated(dg_istart1_all)) deallocate(dg_istart1_all)
      if(allocated(dg_iend1_all)) deallocate(dg_iend1_all)
      if(allocated(dg_istart2_all)) deallocate(dg_istart2_all)
      if(allocated(dg_iend2_all)) deallocate(dg_iend2_all)
      if(allocated(dg_isize_all)) deallocate(dg_isize_all)

      if(allocated(dg_jstart1_all)) deallocate(dg_jstart1_all)
      if(allocated(dg_jend1_all)) deallocate(dg_jend1_all)
      if(allocated(dg_jstart2_all)) deallocate(dg_jstart2_all)
      if(allocated(dg_jend2_all)) deallocate(dg_jend2_all)
      if(allocated(dg_jsize_all)) deallocate(dg_jsize_all)

      if(allocated(dg_kstart1_all)) deallocate(dg_kstart1_all)
      if(allocated(dg_kend1_all)) deallocate(dg_kend1_all)
      if(allocated(dg_kstart2_all)) deallocate(dg_kstart2_all)
      if(allocated(dg_kend2_all)) deallocate(dg_kend2_all)
      if(allocated(dg_ksize_all)) deallocate(dg_ksize_all)

!      if(allocated(dg_dx_all)) deallocate(dg_dx_all)
!      if(allocated(dg_dy_all)) deallocate(dg_dy_all)
!      if(allocated(dg_dz_all)) deallocate(dg_dz_all)

      if(allocated(dg_cycoffset)) deallocate(dg_cycoffset)
      if(allocated(icycoffset)) deallocate(icycoffset)

      if(allocated(dg_c1_all)) deallocate(dg_c1_all)
      if(allocated(dg_c2_all)) deallocate(dg_c2_all)
      if(allocated(dg_c3_all)) deallocate(dg_c3_all)            

      if(allocated(dg_pic)) deallocate(dg_pic)
      
      
! des/mpi_init_des_mod.f
      if(allocated(dsendbuf)) deallocate(dsendbuf)
      if(allocated(drecvbuf)) deallocate(drecvbuf)
      if(allocated(isendindices)) deallocate(isendindices)
      if(allocated(irecvindices)) deallocate(irecvindices)                  
      if(allocated(isendreq)) deallocate(isendreq)
      if(allocated(irecvreq)) deallocate(irecvreq)
      if(allocated(isendcnt)) deallocate(isendcnt)
      if(allocated(dcycl_offset)) deallocate(dcycl_offset)                  
      if(allocated(ineighproc)) deallocate(ineighproc)
      if(allocated(iexchflag)) deallocate(iexchflag)
      if(allocated(iscattercnts)) deallocate(iscattercnts)
      if(allocated(igathercnts)) deallocate(igathercnts)                  
      if(allocated(idispls)) deallocate(idispls)   
      

      RETURN
      END SUBROUTINE DEALLOCATE_ARRAYS_PARALLEL
      
      SUBROUTINE DEALLOCATE_CUT_CELL_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: ALLOCATE_ARRAYS
!  Purpose: allocate arrays
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      Use indices

      USE cutcell
      USE stl
      USE discretelement

      IMPLICIT NONE

      if(allocated(INTERIOR_CELL_AT)) Deallocate(  INTERIOR_CELL_AT   )

      if(allocated(XG_E)) Deallocate( XG_E )
      if(allocated(	 YG_N 	)) deallocate(	 YG_N 	)
      if(allocated(	 ZG_T 	)) deallocate(	 ZG_T 	)

      if(allocated(	  X_U  	)) deallocate(	  X_U  	)
      if(allocated(	  Y_U  	)) deallocate(	  Y_U  	)
      if(allocated(	  Z_U  	)) deallocate(	  Z_U  	)

      if(allocated(	  X_V  	)) deallocate(	  X_V  	)
      if(allocated(	  Y_V  	)) deallocate(	  Y_V  	)
      if(allocated(	  Z_V  	)) deallocate(	  Z_V  	)

      if(allocated(	  X_W  	)) deallocate(	  X_W  	)
      if(allocated(	  Y_W  	)) deallocate(	  Y_W  	)
      if(allocated(	  Z_W  	)) deallocate(	  Z_W  	)

      if(allocated(	  INTERSECT_X   	)) deallocate(	  INTERSECT_X   	)
      if(allocated(	  INTERSECT_Y   	)) deallocate(	  INTERSECT_Y   	)
      if(allocated(	  INTERSECT_Z   	)) deallocate(	  INTERSECT_Z   	)

      if(allocated(	  X_int  	)) deallocate(	  X_int  	)
      if(allocated(	  Y_int  	)) deallocate(	  Y_int  	)
      if(allocated(	  Z_int  	)) deallocate(	  Z_int  	)

      if(allocated(	  X_NEW_POINT   	)) deallocate(	  X_NEW_POINT   	)
      if(allocated(	  Y_NEW_POINT   	)) deallocate(	  Y_NEW_POINT   	)
      if(allocated(	  Z_NEW_POINT   	)) deallocate(	  Z_NEW_POINT   	)

      if(allocated(	  X_NEW_U_POINT   	)) deallocate(	  X_NEW_U_POINT   	)
      if(allocated(	  Y_NEW_U_POINT   	)) deallocate(	  Y_NEW_U_POINT   	)
      if(allocated(	  Z_NEW_U_POINT   	)) deallocate(	  Z_NEW_U_POINT   	)

      if(allocated(	  X_NEW_V_POINT   	)) deallocate(	  X_NEW_V_POINT   	)
      if(allocated(	  Y_NEW_V_POINT   	)) deallocate(	  Y_NEW_V_POINT   	)
      if(allocated(	  Z_NEW_V_POINT   	)) deallocate(	  Z_NEW_V_POINT   	)

      if(allocated(	  X_NEW_W_POINT   	)) deallocate(	  X_NEW_W_POINT   	)
      if(allocated(	  Y_NEW_W_POINT   	)) deallocate(	  Y_NEW_W_POINT   	)
      if(allocated(	  Z_NEW_W_POINT   	)) deallocate(	  Z_NEW_W_POINT   	)

      if(allocated(	  NUMBER_OF_NODES   	)) deallocate(	  NUMBER_OF_NODES   	)
      if(allocated(	  NUMBER_OF_U_NODES   	)) deallocate(	  NUMBER_OF_U_NODES   	)
      if(allocated(	  NUMBER_OF_V_NODES   	)) deallocate(	  NUMBER_OF_V_NODES   	)
      if(allocated(	  NUMBER_OF_W_NODES   	)) deallocate(	  NUMBER_OF_W_NODES   	)

      if(allocated(	  CONNECTIVITY  	)) deallocate(	  CONNECTIVITY  	)
      if(allocated(	  CONNECTIVITY_U  	)) deallocate(	  CONNECTIVITY_U  	)
      if(allocated(	  CONNECTIVITY_V  	)) deallocate(	  CONNECTIVITY_V  	)
      if(allocated(	  CONNECTIVITY_W  	)) deallocate(	  CONNECTIVITY_W  	)

      if(allocated(	  PARTITION   	)) deallocate(	  PARTITION   	)

      if(allocated(	  WALL_U_AT  	)) deallocate(	  WALL_U_AT  	)
      if(allocated(	  WALL_V_AT  	)) deallocate(	  WALL_V_AT  	)
      if(allocated(	  WALL_W_AT  	)) deallocate(	  WALL_W_AT  	)

      if(allocated(	 Area_CUT   	)) deallocate(	 Area_CUT   	)
      if(allocated(	 Area_U_CUT   	)) deallocate(	 Area_U_CUT   	)
      if(allocated(	 Area_V_CUT   	)) deallocate(	 Area_V_CUT   	)
      if(allocated(	 Area_W_CUT   	)) deallocate(	 Area_W_CUT   	)


      if(allocated(	 DELX_Ue   	)) deallocate(	 DELX_Ue   	)
      if(allocated(	 DELX_Uw   	)) deallocate(	 DELX_Uw   	)
      if(allocated(	 DELY_Un   	)) deallocate(	 DELY_Un   	)
      if(allocated(	 DELY_Us   	)) deallocate(	 DELY_Us   	)
      if(allocated(	 DELZ_Ut   	)) deallocate(	 DELZ_Ut   	)
      if(allocated(	 DELZ_Ub   	)) deallocate(	 DELZ_Ub   	)

      if(allocated(	 DELX_Ve   	)) deallocate(	 DELX_Ve   	)
      if(allocated(	 DELX_Vw   	)) deallocate(	 DELX_Vw   	)
      if(allocated(	 DELY_Vn   	)) deallocate(	 DELY_Vn   	)
      if(allocated(	 DELY_Vs   	)) deallocate(	 DELY_Vs   	)
      if(allocated(	 DELZ_Vt   	)) deallocate(	 DELZ_Vt   	)
      if(allocated(	 DELZ_Vb   	)) deallocate(	 DELZ_Vb   	)

      if(allocated(	 DELX_We   	)) deallocate(	 DELX_We   	)
      if(allocated(	 DELX_Ww   	)) deallocate(	 DELX_Ww   	)
      if(allocated(	 DELY_Wn   	)) deallocate(	 DELY_Wn   	)
      if(allocated(	 DELY_Ws   	)) deallocate(	 DELY_Ws   	)
      if(allocated(	 DELZ_Wt   	)) deallocate(	 DELZ_Wt   	)
      if(allocated(	 DELZ_Wb   	)) deallocate(	 DELZ_Wb   	)

      if(allocated(	 X_U_ec   	)) deallocate(	 X_U_ec   	)
      if(allocated(	 Y_U_ec   	)) deallocate(	 Y_U_ec   	)
      if(allocated(	 Z_U_ec   	)) deallocate(	 Z_U_ec   	)
      if(allocated(	 X_U_nc   	)) deallocate(	 X_U_nc   	)
      if(allocated(	 Y_U_nc   	)) deallocate(	 Y_U_nc   	)
      if(allocated(	 Z_U_nc   	)) deallocate(	 Z_U_nc   	)
      if(allocated(	 X_U_tc   	)) deallocate(	 X_U_tc   	)
      if(allocated(	 Y_U_tc   	)) deallocate(	 Y_U_tc   	)
      if(allocated(	 Z_U_tc   	)) deallocate(	 Z_U_tc   	)

      if(allocated(	 X_V_ec   	)) deallocate(	 X_V_ec   	)
      if(allocated(	 Y_V_ec   	)) deallocate(	 Y_V_ec   	)
      if(allocated(	 Z_V_ec   	)) deallocate(	 Z_V_ec   	)
      if(allocated(	 X_V_nc   	)) deallocate(	 X_V_nc   	)
      if(allocated(	 Y_V_nc   	)) deallocate(	 Y_V_nc   	)
      if(allocated(	 Z_V_nc   	)) deallocate(	 Z_V_nc   	)
      if(allocated(	 X_V_tc   	)) deallocate(	 X_V_tc   	)
      if(allocated(	 Y_V_tc   	)) deallocate(	 Y_V_tc   	)
      if(allocated(	 Z_V_tc   	)) deallocate(	 Z_V_tc   	)

      if(allocated(	 X_W_ec   	)) deallocate(	 X_W_ec   	)
      if(allocated(	 Y_W_ec   	)) deallocate(	 Y_W_ec   	)
      if(allocated(	 Z_W_ec   	)) deallocate(	 Z_W_ec   	)
      if(allocated(	 X_W_nc   	)) deallocate(	 X_W_nc   	)
      if(allocated(	 Y_W_nc   	)) deallocate(	 Y_W_nc   	)
      if(allocated(	 Z_W_nc   	)) deallocate(	 Z_W_nc   	)
      if(allocated(	 X_W_tc   	)) deallocate(	 X_W_tc   	)
      if(allocated(	 Y_W_tc   	)) deallocate(	 Y_W_tc   	)
      if(allocated(	 Z_W_tc   	)) deallocate(	 Z_W_tc   	)
       
      if(allocated(	 DELH_Scalar   	)) deallocate(	 DELH_Scalar   	)
       
      if(allocated(	 DELH_U   	)) deallocate(	 DELH_U   	)
      if(allocated(	 Theta_Ue   	)) deallocate(	 Theta_Ue   	)
      if(allocated(	 Theta_Ue_bar  	)) deallocate(	 Theta_Ue_bar  	)
      if(allocated(	 Theta_U_ne   	)) deallocate(	 Theta_U_ne   	)
      if(allocated(	 Theta_U_nw   	)) deallocate(	 Theta_U_nw   	)
      if(allocated(	 Theta_U_te   	)) deallocate(	 Theta_U_te   	)
      if(allocated(	 Theta_U_tw   	)) deallocate(	 Theta_U_tw   	)
      if(allocated(	 ALPHA_Ue_c   	)) deallocate(	 ALPHA_Ue_c   	)
      if(allocated(	 NOC_U_E   	)) deallocate(	 NOC_U_E   	)
      if(allocated(	 Theta_Un   	)) deallocate(	 Theta_Un   	)
      if(allocated(	 Theta_Un_bar  	)) deallocate(	 Theta_Un_bar  	)
      if(allocated(	 ALPHA_Un_c   	)) deallocate(	 ALPHA_Un_c   	)
      if(allocated(	 NOC_U_N   	)) deallocate(	 NOC_U_N   	)
      if(allocated(	 Theta_Ut   	)) deallocate(	 Theta_Ut   	)
      if(allocated(	 Theta_Ut_bar  	)) deallocate(	 Theta_Ut_bar  	)
      if(allocated(	 ALPHA_Ut_c   	)) deallocate(	 ALPHA_Ut_c   	)
      if(allocated(	 NOC_U_T   	)) deallocate(	 NOC_U_T   	)
      if(allocated(	 A_UPG_E  	)) deallocate(	 A_UPG_E  	)
      if(allocated(	 A_UPG_W  	)) deallocate(	 A_UPG_W  	)

      if(allocated(	 DELH_V   	)) deallocate(	 DELH_V   	)
      if(allocated(	 Theta_V_ne   	)) deallocate(	 Theta_V_ne   	)
      if(allocated(	 Theta_V_se   	)) deallocate(	 Theta_V_se   	)
      if(allocated(	 Theta_Vn   	)) deallocate(	 Theta_Vn   	)
      if(allocated(	 Theta_Vn_bar  	)) deallocate(	 Theta_Vn_bar  	)
      if(allocated(	 Theta_V_nt   	)) deallocate(	 Theta_V_nt   	)
      if(allocated(	 Theta_V_st  	)) deallocate(	 Theta_V_st  	)
      if(allocated(	 Theta_Ve   	)) deallocate(	 Theta_Ve   	)
      if(allocated(	 Theta_Ve_bar  	)) deallocate(	 Theta_Ve_bar  	)
      if(allocated(	 ALPHA_Ve_c   	)) deallocate(	 ALPHA_Ve_c   	)
      if(allocated(	 NOC_V_E   	)) deallocate(	 NOC_V_E   	)
      if(allocated(	 ALPHA_Vn_c   	)) deallocate(	 ALPHA_Vn_c   	)
      if(allocated(	 NOC_V_N   	)) deallocate(	 NOC_V_N   	)
      if(allocated(	 Theta_Vt   	)) deallocate(	 Theta_Vt   	)
      if(allocated(	 Theta_Vt_bar  	)) deallocate(	 Theta_Vt_bar  	)
      if(allocated(	 ALPHA_Vt_c   	)) deallocate(	 ALPHA_Vt_c   	)
      if(allocated(	 NOC_V_T   	)) deallocate(	 NOC_V_T   	)
      if(allocated(	 A_VPG_N  	)) deallocate(	 A_VPG_N  	)
      if(allocated(	 A_VPG_S  	)) deallocate(	 A_VPG_S  	)

      if(allocated(	 DELH_W  	)) deallocate(	 DELH_W  	)
      if(allocated(	 Theta_W_te  	)) deallocate(	 Theta_W_te  	)
      if(allocated(	 Theta_W_be  	)) deallocate(	 Theta_W_be  	)
      if(allocated(	 Theta_W_tn  	)) deallocate(	 Theta_W_tn  	)
      if(allocated(	 Theta_W_bn  	)) deallocate(	 Theta_W_bn  	)
      if(allocated(	 Theta_Wt  	)) deallocate(	 Theta_Wt  	)
      if(allocated(	 Theta_Wt_bar  	)) deallocate(	 Theta_Wt_bar  	)
      if(allocated(	 Theta_We  	)) deallocate(	 Theta_We  	)
      if(allocated(	 Theta_We_bar  	)) deallocate(	 Theta_We_bar  	)
      if(allocated(	 ALPHA_We_c  	)) deallocate(	 ALPHA_We_c  	)
      if(allocated(	 NOC_W_E  	)) deallocate(	 NOC_W_E  	)
      if(allocated(	 Theta_Wn  	)) deallocate(	 Theta_Wn  	)
      if(allocated(	 Theta_Wn_bar  	)) deallocate(	 Theta_Wn_bar  	)
      if(allocated(	 ALPHA_Wn_c  	)) deallocate(	 ALPHA_Wn_c  	)
      if(allocated(	 NOC_W_N  	)) deallocate(	 NOC_W_N  	)
      if(allocated(	 ALPHA_Wt_c  	)) deallocate(	 ALPHA_Wt_c  	)
      if(allocated(	 NOC_W_T  	)) deallocate(	 NOC_W_T  	)
      if(allocated(	 A_WPG_T  	)) deallocate(	 A_WPG_T  	)
      if(allocated(	 A_WPG_B  	)) deallocate(	 A_WPG_B  	)

      if(allocated(	 NORMAL_S  	)) deallocate(	 NORMAL_S  	)
      if(allocated(	 NORMAL_U  	)) deallocate(	 NORMAL_U  	)
      if(allocated(	 NORMAL_V  	)) deallocate(	 NORMAL_V  	)
      if(allocated(	 NORMAL_W  	)) deallocate(	 NORMAL_W  	)

      if(allocated(	 REFP_S  	)) deallocate(	 REFP_S  	)
      if(allocated(	 REFP_U  	)) deallocate(	 REFP_U  	)
      if(allocated(	 REFP_V  	)) deallocate(	 REFP_V  	)
      if(allocated(	 REFP_W  	)) deallocate(	 REFP_W  	)
      
      
      

      if(allocated(ONEoDX_E_U)) Deallocate(  ONEoDX_E_U  )
      if(allocated(ONEoDY_N_U)) Deallocate(  ONEoDY_N_U  )
      if(allocated(ONEoDZ_T_U)) Deallocate(  ONEoDZ_T_U  )

      if(allocated(ONEoDX_E_V )) Deallocate(  ONEoDX_E_V  )
      if(allocated(ONEoDY_N_V)) Deallocate(  ONEoDY_N_V  )
      if(allocated(ONEoDZ_T_V)) Deallocate(  ONEoDZ_T_V  )

      if(allocated(ONEoDX_E_W)) Deallocate(  ONEoDX_E_W  )
      if(allocated(ONEoDY_N_W)) Deallocate(  ONEoDY_N_W  )
      if(allocated(ONEoDZ_T_W)) Deallocate(  ONEoDZ_T_W  )

      if(allocated(Xn_int)) Deallocate(  Xn_int  )
      if(allocated(Xn_U_int)) Deallocate(  Xn_U_int  )
      if(allocated(Xn_V_int)) Deallocate(  Xn_V_int  )
      if(allocated(Xn_W_int)) Deallocate(  Xn_W_int  )

      if(allocated(Ye_int)) Deallocate(  Ye_int  )
      if(allocated(Ye_U_int)) Deallocate(  Ye_U_int  )
      if(allocated(Ye_V_int)) Deallocate(  Ye_V_int  )
      if(allocated(Ye_W_int)) Deallocate(  Ye_W_int  )

      if(allocated(Zt_int)) Deallocate(  Zt_int  )
      if(allocated(Zt_U_int)) Deallocate(  Zt_U_int  )
      if(allocated(Zt_V_int)) Deallocate(  Zt_V_int  )
      if(allocated(Zt_W_int)) Deallocate(  Zt_W_int  )

      if(allocated(SNAP)) Deallocate(  SNAP  )

      if(allocated(CUT_TREATMENT_AT)) Deallocate(  CUT_TREATMENT_AT  )
      if(allocated(CUT_U_TREATMENT_AT)) Deallocate(  CUT_U_TREATMENT_AT  )
      if(allocated(CUT_V_TREATMENT_AT)) Deallocate(  CUT_V_TREATMENT_AT  )
      if(allocated(CUT_W_TREATMENT_AT)) Deallocate(  CUT_W_TREATMENT_AT  )

      if(allocated(CUT_CELL_AT)) Deallocate(  CUT_CELL_AT  )
      if(allocated(CUT_U_CELL_AT)) Deallocate(  CUT_U_CELL_AT  )
      if(allocated(CUT_V_CELL_AT)) Deallocate(  CUT_V_CELL_AT  )
      if(allocated(CUT_W_CELL_AT)) Deallocate(  CUT_W_CELL_AT  )

      if(allocated(SMALL_CELL_AT)) Deallocate( SMALL_CELL_AT   )

      if(allocated(SMALL_CELL_FLAG)) Deallocate( SMALL_CELL_FLAG   )

      if(allocated(BLOCKED_CELL_AT)) Deallocate(  BLOCKED_CELL_AT  )
      if(allocated(BLOCKED_U_CELL_AT)) Deallocate(  BLOCKED_U_CELL_AT  )
      if(allocated(BLOCKED_V_CELL_AT)) Deallocate(  BLOCKED_V_CELL_AT  )
      if(allocated(BLOCKED_W_CELL_AT)) Deallocate(  BLOCKED_W_CELL_AT  )

      if(allocated(STANDARD_CELL_AT)) Deallocate(  STANDARD_CELL_AT  )
      if(allocated(STANDARD_U_CELL_AT)) Deallocate(  STANDARD_U_CELL_AT  )
      if(allocated(STANDARD_V_CELL_AT)) Deallocate(  STANDARD_V_CELL_AT  )
      if(allocated(STANDARD_W_CELL_AT)) Deallocate(  STANDARD_W_CELL_AT  )


      if(allocated(VORTICITY )) Deallocate(  VORTICITY  )
      if(allocated(LAMBDA2)) Deallocate(  LAMBDA2  )

      if(allocated(TRD_G_OUT )) Deallocate(  TRD_G_OUT  )
      if(allocated(PP_G_OUT)) Deallocate(  PP_G_OUT  )
      if(allocated(EPP_OUT)) Deallocate(  EPP_OUT  )

      if(allocated(dudx_OUT)) Deallocate(  dudx_OUT  )
      if(allocated(dvdy_OUT)) Deallocate(  dvdy_OUT  )
      if(allocated(delv_OUT)) Deallocate(  delv_OUT  )

      if(allocated(U_MASTER_OF)) Deallocate(  U_MASTER_OF  )
      if(allocated(V_MASTER_OF)) Deallocate(  V_MASTER_OF  )
      if(allocated(W_MASTER_OF)) Deallocate(  W_MASTER_OF  )

      if(allocated(BC_ID)) Deallocate(  BC_ID ) 
      if(allocated(BC_U_ID)) Deallocate(  BC_U_ID  )
      if(allocated(BC_V_ID)) Deallocate(  BC_V_ID  )
      if(allocated(BC_W_ID)) Deallocate(  BC_W_ID  )

      if(allocated(DEBUG_CG)) Deallocate(  DEBUG_CG  )

      if(allocated(U_g_CC)) Deallocate(  U_g_CC  )
      if(allocated(V_g_CC)) Deallocate(  V_g_CC  )
      if(allocated(W_g_CC)) Deallocate(  W_g_CC  )

      if(allocated(U_s_CC)) Deallocate(  U_s_CC  )
      if(allocated(V_s_CC)) Deallocate(  V_s_CC  )
      if(allocated(W_s_CC )) Deallocate(  W_s_CC  )

      if(allocated(N_FACET_AT)) Deallocate(N_FACET_AT)
   

      if(allocated(LIST_FACET_AT)) Deallocate(LIST_FACET_AT)

      if(allocated(POTENTIAL_CUT_CELL_AT)) Deallocate(POTENTIAL_CUT_CELL_AT)


      if(allocated(F_AT)) Deallocate(  F_AT  )

      if(allocated(DWALL)) Deallocate(  DWALL )


      RETURN
      END SUBROUTINE DEALLOCATE_CUT_CELL_ARRAYS    


      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEALLOCATE_DEM_MIO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DEALLOCATE_DEM_MI

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: undefined
      USE des_bc, only: dem_bcmi
      USE des_bc, only: pi_factor, pi_count
      use des_bc, only: numfrac_limit
      use des_bc, only: dem_mi_time, dem_bc_poly_layout
      use des_bc, only: dem_mi
      use des_bc, only: dem_bcmi_ijkstart, dem_bcmi_ijkend
      IMPLICIT NONE
!-----------------------------------------------

! Particle injection factor
      if(allocated(	 PI_FACTOR 	)) deallocate(	 PI_FACTOR 	)
! Particle injection count (injection number)
      if(allocated(	 PI_COUNT 	)) deallocate(	 PI_COUNT 	)
! Particle injection time scale
      if(allocated(	 DEM_MI_TIME 	)) deallocate(	 DEM_MI_TIME 	)
! Array used for polydisperse inlets: stores the particle number
! distribution of an inlet scaled with numfrac_limit
      if(allocated(	 DEM_BC_POLY_LAYOUT	)) deallocate(	 DEM_BC_POLY_LAYOUT	)
! Data structure for storing BC data.
      if(allocated(	 DEM_MI	)) deallocate(	 DEM_MI	)

      if(allocated(	 DEM_BCMI_IJKSTART	)) deallocate(	 DEM_BCMI_IJKSTART	)
      if(allocated(	 DEM_BCMI_IJKEND	)) deallocate(	 DEM_BCMI_IJKEND	)

 

! Boundary classification
!         Allocate( PARTICLE_PLCMNT (DES_BCMI) )
! Character precision arrays
!         PARTICLE_PLCMNT(:) = UNDEFINED_C

      RETURN
      END SUBROUTINE DEALLOCATE_DEM_MI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEALLOCATE_PIC_MIO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: R. Garg                                    Date: 11-Jun-14  !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEALLOCATE_PIC_MIO

! Modules
!-----------------------------------------------
      USE pic_bc, only: pic_bcmi, pic_bcmo
      USE pic_bc, only: pic_bcmi_ijkstart, pic_bcmi_ijkend
      USE pic_bc, only: pic_bcmo_ijkstart, pic_bcmo_ijkend
      USE pic_bc, only: pic_bcmi_normdir
      USE pic_bc, only: pic_bcmi_offset
      USE pic_bc, only: pic_bcmi_incl_cutcell
      IMPLICIT NONE
!-----------------------------------------------

! Allocate/Initialize for inlets

      if(allocated(	 PIC_BCMI_IJKSTART	)) deallocate(	 PIC_BCMI_IJKSTART	)
      if(allocated(	 PIC_BCMI_IJKEND  	)) deallocate(	 PIC_BCMI_IJKEND  	)
      if(allocated(	 PIC_BCMI_NORMDIR 	)) deallocate(	 PIC_BCMI_NORMDIR 	)
      if(allocated(	 PIC_BCMI_OFFSET  	)) deallocate(	 PIC_BCMI_OFFSET  	)
      if(allocated(	 PIC_BCMI_INCL_CUTCELL	)) deallocate(	 PIC_BCMI_INCL_CUTCELL	)

      if(allocated(	 PIC_BCMO_IJKSTART	)) deallocate(	 PIC_BCMO_IJKSTART	)
      if(allocated(	 PIC_BCMO_IJKEND	)) deallocate(	 PIC_BCMO_IJKEND	)

      RETURN
      END SUBROUTINE DEALLOCATE_PIC_MIO


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DEALLOCATE_ARRAYS                                   C
!  Purpose: Deallocate arrays subroutines for DES                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DES_DEALLOCATE_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE cutcell
      USE derived_types, only: boxhandle, pic
      USE des_bc
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE functions
      USE funits
      USE geometry
      USE indices
      USE mfix_pic
      USE param
      USE param1
      USE physprop
      USE pic_bc, only: pic_bcmo, pic_bcmi

      USE run, only: ENERGY_EQ
      USE run, only: ANY_SPECIES_EQ

      USE particle_filter, only: DES_INTERP_SCHEME_ENUM
      USE particle_filter, only: DES_INTERP_GARG
      USE particle_filter, only: DES_INTERP_DPVM
      USE particle_filter, only: DES_INTERP_GAUSS
      USE particle_filter, only: DES_INTERP_LHAT
      USE particle_filter, only: FILTER_SIZE
      USE particle_filter, only: FILTER_CELL
      USE particle_filter, only: FILTER_WEIGHT
      use des_bc, only: DEM_BCMO_IJKSTART, DEM_BCMO_IJKEND
      use stl, only: FACETS_AT_DG
! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      USE error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: IJK
!-----------------------------------------------

      CALL INIT_ERR_MSG("DES_DEALLOCATE_ARRAYS")


! DES Allocatable arrays
!-----------------------------------------------
! Dynamic particle info including another index for parallel
! processing for ghost
      if(allocated(	 PARTICLE_STATE 	)) deallocate(	 PARTICLE_STATE 	)
      if(allocated(	iglobal_id	)) deallocate(	iglobal_id	)

! R.Garg: Allocate necessary arrays for PIC mass inlet/outlet BCs
            if(PIC_BCMI /= 0 .OR. PIC_BCMO /=0) CALL DEALLOCATE_PIC_MIO

! Particle attributes
! Radius, density, mass, moment of inertia
      if(allocated(	  DES_RADIUS 	)) deallocate(	  DES_RADIUS 	)
      if(allocated(	  RO_Sol 	)) deallocate(	  RO_Sol 	)
      if(allocated(	  PVOL 	)) deallocate(	  PVOL 	)
      if(allocated(	  PMASS 	)) deallocate(	  PMASS 	)
      if(allocated(	  OMOI 	)) deallocate(	  OMOI 	)

! Old and new particle positions, velocities (translational and
! rotational)
      if(allocated(	  DES_POS_NEW 	)) deallocate(	  DES_POS_NEW 	)
      if(allocated(	  DES_VEL_NEW 	)) deallocate(	  DES_VEL_NEW 	)
      if(allocated(	  OMEGA_NEW 	)) deallocate(	  OMEGA_NEW 	)

      if(allocated(	  ORIENTATION 	)) deallocate(	  ORIENTATION 	)

      if(allocated(	  DES_POS_OLD 	)) deallocate(	  DES_POS_OLD 	)
      if(allocated(	  DES_VEL_OLD 	)) deallocate(	  DES_VEL_OLD 	)
      if(allocated(	  DES_ACC_OLD 	)) deallocate(	  DES_ACC_OLD 	)
      if(allocated(	  OMEGA_OLD 	)) deallocate(	  OMEGA_OLD 	)
      if(allocated(	  ROT_ACC_OLD 	)) deallocate(	  ROT_ACC_OLD 	)


! Allocating user defined array
      if(allocated(	 DES_USR_VAR	)) deallocate(	 DES_USR_VAR	)

! Particle positions at the last call neighbor search algorithm call
      if(allocated(	  PPOS 	)) deallocate(	  PPOS 	)

! Total, normal and tangetial forces
      if(allocated(	  FC 	)) deallocate(	  FC 	)

! Torque
      if(allocated(	  TOW 	)) deallocate(	  TOW 	)


! allocate variable for des grid binning
      if(allocated(	dg_pijk	)) deallocate(	dg_pijk	)
      if(allocated(	dg_pijkprv	)) deallocate(	dg_pijkprv	)

! allocate variables related to ghost particles
      if(allocated(	ighost_updated	)) deallocate(	ighost_updated	)


      if(allocated(	  wall_collision_facet_id 	)) deallocate(	  wall_collision_facet_id 	)

      if(allocated(	  wall_collision_PFT 	)) deallocate(	  wall_collision_PFT 	)

! Temporary variables to store wall position, velocity and normal vector
      if(allocated(	  WALL_NORMAL  	)) deallocate(	  WALL_NORMAL  	)

      if(allocated(	  NEIGHBOR_INDEX 	)) deallocate(	  NEIGHBOR_INDEX 	)
      if(allocated(	  NEIGHBOR_INDEX_OLD 	)) deallocate(	  NEIGHBOR_INDEX_OLD 	)
      if(allocated(	  NEIGHBORS 	)) deallocate(	  NEIGHBORS 	)

      if(allocated(	  NEIGHBORS_OLD 	)) deallocate(	  NEIGHBORS_OLD 	)
      if(allocated(	  PFT_NEIGHBOR 	)) deallocate(	  PFT_NEIGHBOR 	)
      if(allocated(	  PFT_NEIGHBOR_OLD 	)) deallocate(	  PFT_NEIGHBOR_OLD 	)

#ifdef do_sap
      if(allocated(	  boxhandle	)) deallocate(	  boxhandle	)
#endif

! Variable that stores the particle in cell information (ID) on the
! computational fluid grid defined by imax, jmax and kmax in mfix.dat
      if(allocated(	PIC	)) deallocate(	PIC	)

! Particles in a computational fluid cell (for volume fraction)
      if(allocated(	  PINC 	)) deallocate(	  PINC 	)

! For each particle track its i,j,k location on computational fluid grid
! defined by imax, jmax and kmax in mfix.dat and phase no.
      if(allocated(	  PIJK 	)) deallocate(	  PIJK 	)

      if(allocated(	DRAG_AM	)) deallocate(	DRAG_AM	)
      if(allocated(	DRAG_BM	)) deallocate(	DRAG_BM	)
      if(allocated(	F_gp	)) deallocate(	F_gp	)


! Explicit drag force acting on a particle.
      if(allocated(	DRAG_FC 	)) deallocate(	DRAG_FC 	)

! force due to gas-pressure gradient
      if(allocated(	P_FORCE	)) deallocate(	P_FORCE	)

! Volume of nodes
      if(allocated(	DES_VOL_NODE	)) deallocate(	DES_VOL_NODE	)

      if(allocated(	F_GDS	)) deallocate(	F_GDS	)
      if(allocated(	VXF_GDS	)) deallocate(	VXF_GDS	)

      if(allocated(	FILTER_CELL	)) deallocate(	FILTER_CELL	)
      if(allocated(	FILTER_WEIGHT	)) deallocate(	FILTER_WEIGHT	)
      if(allocated(	DES_ROPS_NODE	)) deallocate(	DES_ROPS_NODE	)
      if(allocated(	DES_VEL_NODE	)) deallocate(	DES_VEL_NODE	)

! Variables for hybrid model

      if(allocated(	SDRAG_AM	)) deallocate(	SDRAG_AM	)
      if(allocated(	SDRAG_BM	)) deallocate(	SDRAG_BM	)

      if(allocated(	F_SDS	)) deallocate(	F_SDS	)
      if(allocated(	VXF_SDS	)) deallocate(	VXF_SDS	)


! MP-PIC related
      if(allocated(	PS_FORCE_PIC	)) deallocate(	PS_FORCE_PIC	)
      if(allocated(	DES_STAT_WT	)) deallocate(	DES_STAT_WT	)
      if(allocated(	DES_VEL_MAX	)) deallocate(	DES_VEL_MAX	)
      if(allocated(	PS_GRAD	)) deallocate(	PS_GRAD	)
      if(allocated(	AVGSOLVEL_P	)) deallocate(	AVGSOLVEL_P	)
      if(allocated(	EPG_P	)) deallocate(	EPG_P	)

      if(allocated(	PIC_U_s 	)) deallocate(	PIC_U_s 	)
      if(allocated(	PIC_V_s 	)) deallocate(	PIC_V_s 	)
      if(allocated(	PIC_W_s 	)) deallocate(	PIC_W_s 	)
      if(allocated(	PIC_P_s 	)) deallocate(	PIC_P_s 	)

! Averaged velocity obtained by averaging over all the particles
      if(allocated(	DES_VEL_AVG	)) deallocate(	DES_VEL_AVG	)

! Global Granular Energy
      if(allocated(	GLOBAL_GRAN_ENERGY	)) deallocate(	GLOBAL_GRAN_ENERGY	)
      if(allocated(	GLOBAL_GRAN_TEMP	)) deallocate(	GLOBAL_GRAN_TEMP	)

! variable for bed height of solids phase M
      if(allocated(	BED_HEIGHT	)) deallocate(	BED_HEIGHT	)

! ---------------------------------------------------------------->>>
! BEGIN COHESION
! Matrix location of particle  (should be allocated in case user wishes
! to invoke routines in /cohesion subdirectory
      if(allocated(	  PostCohesive 	)) deallocate(	  PostCohesive 	)
! END COHESION
! ----------------------------------------------------------------<<<

! ---------------------------------------------------------------->>>
! BEGIN Thermodynamic Allocation
! Particle temperature
      if(allocated(	 DES_T_s	)) deallocate(	 DES_T_s	)
! Spec      ific heat
      if(allocated(	 DES_C_PS	)) deallocate(	 DES_C_PS	)
! Species mass fractions comprising a particle. This array may not be
! needed for all thermo problems.
      if(allocated(	 DES_X_s	)) deallocate(	 DES_X_s	)
! Total rate of heat transfer to individual particles.
      if(allocated(	 Q_Source	)) deallocate(	 Q_Source	)
! Average solids temperature in fluid cell
      if(allocated(	avgDES_T_s	)) deallocate(	avgDES_T_s	)
! Gas/Solids convective heat transfer coupling

! Fluid phase energy equation source terms
      if(allocated(	CONV_Sc	)) deallocate(	CONV_Sc	)
      if(allocated(	CONV_Sp	)) deallocate(	CONV_Sp	)
! Particle convection source term (explicit coupled)
      if(allocated(	CONV_Qs	)) deallocate(	CONV_Qs	)
! Gas-particle heat transfer coefficient TIMES surface area
      if(allocated(	GAMMAxSA	)) deallocate(	GAMMAxSA	)

! Allocate the history variables for Adams-Bashforth integration
      if(allocated(	 Q_Source0	)) deallocate(	 Q_Source0	)

! End Thermodynamic Allocation
! ----------------------------------------------------------------<<<


! ---------------------------------------------------------------->>>
! BEGIN Species Allocation
! Rate of solids phase production/consumption for each species
      if(allocated(	 DES_R_s	)) deallocate(	 DES_R_s	)

      if(allocated(	 DES_R_gp	)) deallocate(	 DES_R_gp	)
      if(allocated(	 DES_R_gc	)) deallocate(	 DES_R_gc	)
      if(allocated(	 DES_SUM_R_g	)) deallocate(	 DES_SUM_R_g	)
      if(allocated(	 DES_R_PHASE	)) deallocate(	 DES_R_PHASE	)
      if(allocated(	 DES_HOR_g	)) deallocate(	 DES_HOR_g	)


! Allocate the history variables for Adams-Bashforth integration
 
! Rate of change of particle mass
      if(allocated(	 dMdt_OLD	)) deallocate(	 dMdt_OLD	)
! Rate of change of particle mass percent species
      if(allocated(	 dXdt_OLD	)) deallocate(	 dXdt_OLD	)
  

! Energy generation from reaction (cal/sec)
      if(allocated(	 RXNS_Qs	)) deallocate(	 RXNS_Qs	)
  
! End Species Allocation
! ----------------------------------------------------------------<<<

      if(allocated(PARTICLE_STATE)) deallocate(PARTICLE_STATE)
      if(allocated(iglobal_id)) deallocate(iglobal_id)

      if(allocated(DEM_BCMO_IJKSTART)) deallocate(DEM_BCMO_IJKSTART)
      if(allocated(DEM_BCMO_IJKEND)) deallocate(DEM_BCMO_IJKEND)

! stl
      if(allocated(FACETS_AT_DG)) deallocate(FACETS_AT_DG)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE DES_DEALLOCATE_ARRAYS
