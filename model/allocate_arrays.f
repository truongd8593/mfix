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
      Use dsdgmr_a
      Use energy
      Use fldvar
      Use geometry
      Use indices
      Use pgcor
      Use physprop
      Use pscor
      Use residual
      Use rxns
      Use tau_g
      Use tau_s
      Use tmp_array
      Use tmp_array1
      Use trace
      Use visc_g
      Use visc_s
      Use xsi_array


!//EFD 0913
      use compar
      use debug

      IMPLICIT NONE
      
      INTEGER M
!// 300 0912 Relocated set_max2 call to get_data before calling gridmap_init
!//TD 0912      CALL SET_MAX2


!//EFD 0913 many changes to dimension_i,dimension_j,dimension_j
!//         to conform to new allocation method
!//         change dimension_i to istart3:iend3
!//         change dimension_j to jstart3:jend3
!//         change dimension_k to kstart3:kend3
!//
!//      DIMENSION_I   = IMAX2
!//      DIMENSION_J   = JMAX2
!//      DIMENSION_K   = KMAX2
!//      DIMENSION_3   = IJKMAX2

      DIMENSION_I = (IEND3 - ISTART3 + 1)
      DIMENSION_J = (JEND3 - JSTART3 + 1)
      DIMENSION_K = (KEND3 - KSTART3 + 1)
      DIMENSION_3 = DIMENSION_I*DIMENSION_J*DIMENSION_K
     
!//EFD 0913 extra checks
      
      call assert( ijkstart3 .eq. 1, &
        '** allocate_arrays: ijkstart3 .ne. 1, ijkstart3 =  ', ijkstart3 )
      call assert( dimension_3 .eq. (ijkend3 - ijkstart3 + 1), &
	'** allocate_arrays: invalid dimension_3 ' // & 
        ' dimension_3,  ijkend3 - ijkstart3+1 ', &
          dimension_3,  ijkend3 - ijkstart3+1 )
		
      
      DIMENSION_M   = MAX(1, MMAX)
      
      DIMENSION_N_g = 1
      IF(NMAX(0) .NE. UNDEFINED_I)DIMENSION_N_g = NMAX(0)
      
      
      DIMENSION_N_s = 1
      DO M = 1, MMAX
        IF(NMAX(M) .NE. UNDEFINED_I)DIMENSION_N_s = MAX(DIMENSION_N_s, NMAX(M))
        
      END DO
      
      DIMENSION_LM    = (DIMENSION_M * (DIMENSION_M-1) / 2)+1
      DIMENSION_N_all = DIMENSION_N_g + DIMENSION_M * DIMENSION_N_s

!ambm
      Allocate( A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) )
      Allocate( B_m(DIMENSION_3, 0:DIMENSION_M) )


!coeff
      
      Allocate( DENSITY(0:DIMENSION_M) )
      Allocate( SIZE(0:DIMENSION_M) )
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

!geometry
      Allocate(           FLAG (DIMENSION_3) )
      Allocate(           FLAG_E (DIMENSION_3) )
      Allocate(           FLAG_N (DIMENSION_3) )
      Allocate(           FLAG_T (DIMENSION_3) )
      Allocate(           ICBC_FLAG (DIMENSION_3) )
      Allocate(  oDX (istart3:iend3) )
      Allocate(  oDY (jstart3:jend3) )
      Allocate(  oDZ (kstart3:kend3) )
      Allocate(  oDX_E (istart3:iend3) )
      Allocate(  oDY_N (jstart3:jend3) )
      Allocate(  oDZ_T (kstart3:kend3) )
      Allocate(  X (istart3:iend3) )
      Allocate(  X_E (istart3:iend3) )
      Allocate(  oX (istart3:iend3) )
      Allocate(  oX_E (istart3:iend3) )
      Allocate(  Z (kstart3:kend3) )
      Allocate(  Z_T (kstart3:kend3) )
      Allocate(  FX (istart3:iend3) )
      Allocate(  FX_bar (istart3:iend3) )
      Allocate(  FX_E (istart3:iend3) )
      Allocate(  FX_E_bar (istart3:iend3) )
      Allocate(  FY_N (jstart3:jend3) )
      Allocate(  FY_N_bar (jstart3:jend3) )
      Allocate(  FZ_T (kstart3:kend3) )
      Allocate(  FZ_T_bar (kstart3:kend3) )
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
      Allocate(           Im1 (istart3:iend3) )
      Allocate(           Ip1 (istart3:iend3) )
      Allocate(           Jm1 (jstart3:jend3) )
      Allocate(           Jp1 (jstart3:jend3) )
      Allocate(           Km1 (kstart3:kend3) )
      Allocate(           Kp1 (kstart3:kend3) )
      


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
      
!
! array allocation of add on packages, such as linear equation solvers
!

!//EFD 0913  avoid unnecessary allocation for SLAP and other linear solver
!//      CALL Allocate_dsdgmr
!//      CALL Allocate_dslucs
!//      CALL Allocate_dslugm
!//      CALL Allocate_igcg  
     
      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS 
      
      
      
      SUBROUTINE Allocate_dsdgmr
        USE param 
        USE param1
        Use dsdgmr_a
        IMPLICIT NONE
      

        LENW = 2+DIMENSION_3*(NSAVE+7)+NSAVE*(NSAVE+3) 
        NELTMAX = 7*DIMENSION_3 
        Allocate( IA(NELTMAX) )
        Allocate( JA(NELTMAX) )
        Allocate( IWORK(LENIW) )
        Allocate( A(NELTMAX) )
        Allocate( RWORK(LENW) )
      
      END SUBROUTINE Allocate_dsdgmr

      SUBROUTINE Allocate_dslucs
        USE param 
        USE param1
        Use dslucs_a
        IMPLICIT NONE
      

        LENW = 4*DIMENSION_3 + 4*DIMENSION_3 + 8*DIMENSION_3 
         
        LENIW = 4*DIMENSION_3 + 4*DIMENSION_3 + 4*DIMENSION_3 + 12 
         
        NELTMAX = 7*DIMENSION_3 
	
        Allocate( IA(NELTMAX) )
        Allocate( JA(NELTMAX) )
        Allocate( IWORK(LENIW) )
        Allocate( A(NELTMAX) )
        Allocate( RWORK(LENW) )
      
      END SUBROUTINE Allocate_dslucs

      SUBROUTINE Allocate_dslugm
      
        USE param 
        USE param1
        Use dslugm_a
        IMPLICIT NONE

        LENW = 2 + DIMENSION_3*(NSAVE + 7) + NSAVE*(NSAVE&
          + 3) + 4*DIMENSION_3 + 4*DIMENSION_3 
        LENIW = 4*DIMENSION_3 + 4*DIMENSION_3 + 4*&
         DIMENSION_3 + 32 
        NELTMAX = 7*DIMENSION_3 
	
        Allocate( IA(NELTMAX) )
        Allocate( JA(NELTMAX) )
        Allocate( IWORK(LENIW) )
        Allocate( A(NELTMAX) )
        Allocate( RWORK(LENW) )
      
      END SUBROUTINE Allocate_dslugm
      
      SUBROUTINE Allocate_igcg
        USE param 
        USE param1
        Use igcg_i
        Use igcg_a
        IMPLICIT NONE
      

        DIM_I1 = 1 
        DIM_I3 = DIMENSION_I 
        DIM_I9 = DIMENSION_I*DIMENSION_J 
	
        Allocate( BD00N(DIMENSION_3) )
        Allocate( BL09N(DIMENSION_3 + DIM_I9) ) 
        Allocate( BL03N(DIMENSION_3 + DIM_I3) )  
        Allocate( BL01N(DIMENSION_3 + DIM_I1) )  
        Allocate( BU09N(1 - DIM_I9:DIMENSION_3) ) 
        Allocate( BU03N (1 - DIM_I3:DIMENSION_3) )
        Allocate( BU01N (1 - DIM_I1:DIMENSION_3))

      
      END SUBROUTINE Allocate_igcg

