
      SUBROUTINE ALLOCATE_ARRAYS 
      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                      
!  Module name: ALLOCATE_ARRAYS                                     
!  Purpose: allocate arrays
!                                                                      C
!  Author: M. Syamlal                                Date: 17-DEC-98 
!  Reviewer: 
!                                                                     
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
      IMPLICIT NONE
      
      INTEGER M

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
      IF(NMAX(0) .NE. UNDEFINED_I)DIMENSION_N_g = NMAX(0)
      
      
      DIMENSION_N_s = 1
      DO M = 1, MMAX
        IF(NMAX(M) .NE. UNDEFINED_I)DIMENSION_N_s = MAX(DIMENSION_N_s, NMAX(M))
        
      END DO
      
      DIMENSION_LM    = (DIMENSION_M * (DIMENSION_M-1) / 2)+1
      DIMENSION_N_all = DIMENSION_N_g + DIMENSION_M * DIMENSION_N_s
      
      DIMENSION_Scalar = NScalar
!     add by rong
      DIM_Scalar2 = 2*NScalar
!     add by rong
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
      
!cont
      Allocate( DO_CONT(0:DIMENSION_M) )


!drag
      Allocate(  F_gs(DIMENSION_3, DIMENSION_M) )
      Allocate(  F_ss(DIMENSION_3, DIMENSION_LM) )


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
      Allocate(  EP_go (DIMENSION_3) )
      Allocate(  P_g (DIMENSION_3) )
      Allocate(  P_go (DIMENSION_3) )
      Allocate(  RO_g (DIMENSION_3) )
      Allocate(  RO_go (DIMENSION_3) )
      Allocate(  ROP_g (DIMENSION_3) )
      Allocate(  ROP_go (DIMENSION_3) )
      Allocate(  ROP_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  ROP_so (DIMENSION_3, DIMENSION_M) )
      Allocate(  T_g (DIMENSION_3) )
      Allocate(  T_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  T_go (DIMENSION_3) )
      Allocate(  T_so (DIMENSION_3, DIMENSION_M) )
      Allocate(  X_g (DIMENSION_3, DIMENSION_N_g) )
      Allocate(  X_s (DIMENSION_3, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  X_go (DIMENSION_3, DIMENSION_N_g) )
      Allocate(  X_so (DIMENSION_3, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  U_g (DIMENSION_3) )
      Allocate(  U_go (DIMENSION_3) )
      Allocate(  U_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  U_so (DIMENSION_3, DIMENSION_M) )
      Allocate(  V_g (DIMENSION_3) )
      Allocate(  V_go (DIMENSION_3) )
      Allocate(  V_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  V_so (DIMENSION_3, DIMENSION_M) )
      Allocate(  W_g (DIMENSION_3) )
      Allocate(  W_go (DIMENSION_3) )
      Allocate(  W_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  W_so (DIMENSION_3, DIMENSION_M) )
      Allocate(  P_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  P_s_c (DIMENSION_3, DIMENSION_M) )
      Allocate(  P_star (DIMENSION_3) )
      Allocate(  P_staro (DIMENSION_3) )
      Allocate(  THETA_m (DIMENSION_3, DIMENSION_M) )
      Allocate(  THETA_mo (DIMENSION_3, DIMENSION_M) )



      IF(K_Epsilon)then
        Allocate(  K_Turb_G (DIMENSION_3) )
        Allocate(  K_Turb_Go (DIMENSION_3) )
        Allocate(  E_Turb_G (DIMENSION_3) )
        Allocate(  E_Turb_Go (DIMENSION_3) )
      ENDIF
      
      IF(DIMENSION_Scalar /= 0)then
        Allocate(  Scalar (DIMENSION_3,  DIMENSION_Scalar) )
        Allocate(  Scalaro (DIMENSION_3, DIMENSION_Scalar) )
      
      ENDIF


!geometry
      Allocate(           FLAG (DIMENSION_3) )
      Allocate(           FLAG_E (DIMENSION_3) )
      Allocate(           FLAG_N (DIMENSION_3) )
      Allocate(           FLAG_T (DIMENSION_3) )
      Allocate(           ICBC_FLAG (DIMENSION_3L) )
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
      Allocate(  AYZ (DIMENSION_3) )
      Allocate(  AXZ (DIMENSION_3) )
      Allocate(  AXY (DIMENSION_3) )
      Allocate(  VOL (DIMENSION_3) )
      Allocate(  AYZ_U (DIMENSION_3) )
      Allocate(  AXZ_U (DIMENSION_3) )
      Allocate(  AXY_U (DIMENSION_3) )
      Allocate(  VOL_U (DIMENSION_3) )
      Allocate(  AYZ_V (DIMENSION_3) )
      Allocate(  AXZ_V (DIMENSION_3) )
      Allocate(  AXY_V (DIMENSION_3) )
      Allocate(  VOL_V (DIMENSION_3) )
      Allocate(  AYZ_W (DIMENSION_3) )
      Allocate(  AXZ_W (DIMENSION_3) )
      Allocate(  AXY_W (DIMENSION_3) )
      Allocate(  VOL_W (DIMENSION_3) )

!indices
      Allocate(           STORE_LM (DIMENSION_M, DIMENSION_M) )
      Allocate(           CELL_CLASS (DIMENSION_3) )
      Allocate(           I_OF (DIMENSION_3) )
      Allocate(           J_OF (DIMENSION_3) )
      Allocate(           K_OF (DIMENSION_3) )
      Allocate(           Im1 (0:DIMENSION_I) )
      Allocate(           Ip1 (0:DIMENSION_I) )
      Allocate(           Jm1 (0:DIMENSION_J) )
      Allocate(           Jp1 (0:DIMENSION_J) )
      Allocate(           Km1 (0:DIMENSION_K) )
      Allocate(           Kp1 (0:DIMENSION_K) )
      


!pgcor
      Allocate(  d_e(DIMENSION_3, 0:DIMENSION_M) )
      Allocate(  d_n(DIMENSION_3, 0:DIMENSION_M) )
      Allocate(  d_t(DIMENSION_3, 0:DIMENSION_M) )
      Allocate(  Pp_g(DIMENSION_3) )
      Allocate(           PHASE_4_P_g(DIMENSION_3) )

!physprop
      Allocate(  MU_g (DIMENSION_3) )
      Allocate(  C_pg (DIMENSION_3) )
      Allocate(  C_ps (DIMENSION_3, DIMENSION_M) )
      Allocate(  K_g (DIMENSION_3) )
      Allocate(  K_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  Kth_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  Kphi_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DIF_g (DIMENSION_3, DIMENSION_N_g) )
      Allocate(  DIF_s (DIMENSION_3, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  MW_MIX_g (DIMENSION_3) )

!pscor
      Allocate(  e_e(DIMENSION_3) )
      Allocate(  e_n(DIMENSION_3) )
      Allocate(  e_t(DIMENSION_3) )
      Allocate(  K_cp(DIMENSION_3) )
      Allocate(  EPp(DIMENSION_3) )
      Allocate(           PHASE_4_P_s(DIMENSION_3) )

!residual
      Allocate( RESID(NRESID, 0:DIMENSION_M) )
      Allocate( MAX_RESID(NRESID, 0:DIMENSION_M) )
      Allocate( IJK_RESID(NRESID, 0:DIMENSION_M) )
 
!rxns
      if (nRR .gt. 0) Allocate( ReactionRates(DIMENSION_3,nRR) )
      Allocate(  R_gp (DIMENSION_3, DIMENSION_N_g) )
      Allocate(  R_sp (DIMENSION_3, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  RoX_gc (DIMENSION_3, DIMENSION_N_g) )
      Allocate(  RoX_sc (DIMENSION_3, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  SUM_R_g (DIMENSION_3) )
      Allocate(  SUM_R_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  R_phase (DIMENSION_3, DIMENSION_LM+DIMENSION_M-1) )
      Allocate(  MW_all (DIMENSION_N_all) )
      Allocate(  SPECIES_ID2N(DIMENSION_N_all, 2) )
      Allocate(  SPECIES_N2IDg(DIMENSION_N_g) )
      Allocate(  SPECIES_N2IDs(DIMENSION_M, DIMENSION_N_s) )
      
!scalars
      
      IF(DIMENSION_Scalar /= 0)then
        Allocate(  Scalar_c (DIMENSION_3,  DIMENSION_Scalar) )
        Allocate(  Scalar_p (DIMENSION_3,  DIMENSION_Scalar) )
        Allocate(  Dif_Scalar (DIMENSION_3, DIMENSION_Scalar) )
      
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
 
! add by rong for dqmom


!K-Epsilon Turbulence model
      
      IF(K_Epsilon)then
        Allocate(  K_Turb_G_c (DIMENSION_3) )
        Allocate(  K_Turb_G_p (DIMENSION_3) )
        Allocate(  Dif_K_Turb_G (DIMENSION_3) )
        Allocate(  E_Turb_G_c (DIMENSION_3) )
        Allocate(  E_Turb_G_p (DIMENSION_3) )
        Allocate(  Dif_E_Turb_G (DIMENSION_3) )
      
      ENDIF
! Simonin or Ahmadi model

      IF(SIMONIN .OR. AHMADI)then
        Allocate(  K_12 (DIMENSION_3) )
	Allocate(  Tau_12 (DIMENSION_3) )
        Allocate(  Tau_1 (DIMENSION_3) )
      ENDIF

!tau_g
      Allocate(  TAU_U_g(DIMENSION_3) )
      Allocate(  TAU_V_g(DIMENSION_3) )
      Allocate(  TAU_W_g(DIMENSION_3) )

!tau_s
      Allocate(  TAU_U_s(DIMENSION_3, DIMENSION_M) )
      Allocate(  TAU_V_s(DIMENSION_3, DIMENSION_M) )
      Allocate(  TAU_W_s(DIMENSION_3, DIMENSION_M) )
      
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
      
!xsi_array
      Allocate(  Xsi_e(DIMENSION_3) )
      Allocate(  Xsi_n(DIMENSION_3) )
      Allocate(  Xsi_t(DIMENSION_3) )
      

!VSH
      Allocate(  VSH(DIMENSION_3) )

!VSHE
      Allocate(  VSHE(DIMENSION_3) )

!
! array allocation of add on packages, such as linear equation solvers
!

! array allocation for higher order implementation
      Allocate(           FLAG3 (DIMENSION_4) )
      Allocate(           CELL_CLASS3 (DIMENSION_4) )
      Allocate(           I3_OF (DIMENSION_4) )
      Allocate(           J3_OF (DIMENSION_4) )
      Allocate(           K3_OF (DIMENSION_4) )
      Allocate(           Im1_3 (-1:DIMENSION_I+1) )
      Allocate(           Ip1_3 (-1:DIMENSION_I+1) )
      Allocate(           Jm1_3 (-1:DIMENSION_J+1) )
      Allocate(           Jp1_3 (-1:DIMENSION_J+1) )
      Allocate(           Km1_3 (-1:DIMENSION_K+1) )
      Allocate(           Kp1_3 (-1:DIMENSION_K+1) )
 
!mflux
      Allocate(    Flux_gE(DIMENSION_3) ) 
      Allocate(    Flux_sE(DIMENSION_3, DIMENSION_M) ) 
      Allocate(    Flux_gN(DIMENSION_3) ) 
      Allocate(    Flux_sN(DIMENSION_3, DIMENSION_M) ) 
      Allocate(    Flux_gT(DIMENSION_3) ) 
      Allocate(    Flux_sT(DIMENSION_3, DIMENSION_M) ) 
      Allocate(    ROP_gE(DIMENSION_3) ) 
      Allocate(    ROP_sE(DIMENSION_3, DIMENSION_M) ) 
      Allocate(    ROP_gN(DIMENSION_3) ) 
      Allocate(    ROP_sN(DIMENSION_3, DIMENSION_M) ) 
      Allocate(    ROP_gT(DIMENSION_3) ) 
      Allocate(    ROP_sT(DIMENSION_3, DIMENSION_M) ) 

     
      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS 
      

