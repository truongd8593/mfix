
      SUBROUTINE Deallocate_ARRAYS 
      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                      
!  Module name: DEDeallocate_ARRAYS                                     
!  Purpose: deDeallocate arrays
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
      Use scalars
      Use tau_g
      Use tau_s
      Use tmp_array
      Use tmp_array1
      Use trace
      Use visc_g
      Use visc_s
      Use xsi_array
      Use vshear
      IMPLICIT NONE
      
      INTEGER M
      
      
!ambm
      Deallocate( A_m )
      Deallocate( B_m )


!coeff
      
      Deallocate( DENSITY )
      Deallocate( SIZE )
      Deallocate( SP_HEAT )
      Deallocate( VISC )
      Deallocate( COND )
      Deallocate( DIFF )
      Deallocate( DRAGCOEF )
      Deallocate( HEAT_TR)
      
!cont
      Deallocate( DO_CONT )


!drag
      Deallocate(  F_gs )
      Deallocate(  F_ss )


!energy
      Deallocate(  HOR_g  )
      Deallocate(  HOR_s  )
      Deallocate(  GAMA_gs  )
      Deallocate(  GAMA_Rg  )
      Deallocate(  GAMA_Rs  )
      Deallocate(  T_Rg  )
      Deallocate(  T_Rs  )

!fldvar
      Deallocate(  EP_g  )
      Deallocate(  EP_go  )
      Deallocate(  P_g  )
      Deallocate(  P_go  )
      Deallocate(  RO_g  )
      Deallocate(  RO_go  )
      Deallocate(  ROP_g  )
      Deallocate(  ROP_go  )
      Deallocate(  ROP_s  )
      Deallocate(  ROP_so  )
      Deallocate(  T_g  )
      Deallocate(  T_s  )
      Deallocate(  T_go  )
      Deallocate(  T_so  )
      Deallocate(  X_g  )
      Deallocate(  X_s  )
      Deallocate(  X_go  )
      Deallocate(  X_so  )
      Deallocate(  U_g  )
      Deallocate(  U_go  )
      Deallocate(  U_s  )
      Deallocate(  U_so  )
      Deallocate(  V_g  )
      Deallocate(  V_go  )
      Deallocate(  V_s  )
      Deallocate(  V_so  )
      Deallocate(  W_g  )
      Deallocate(  W_go  )
      Deallocate(  W_s  )
      Deallocate(  W_so  )
      Deallocate(  P_s  )
      Deallocate(  P_s_c  )
      Deallocate(  P_star  )
      Deallocate(  P_staro  )
      Deallocate(  THETA_m  )
      Deallocate(  THETA_mo  )
      
      IF(DIMENSION_Scalar /= 0)then
        Deallocate(  Scalar  )
        Deallocate(  Scalaro  )
      
      ENDIF


!geometry
      Deallocate(           FLAG  )
      Deallocate(           FLAG_E  )
      Deallocate(           FLAG_N  )
      Deallocate(           FLAG_T  )
      Deallocate(           ICBC_FLAG  )
      Deallocate(  oDX  )
      Deallocate(  oDY  )
      Deallocate(  oDZ  )
      Deallocate(  oDX_E  )
      Deallocate(  oDY_N  )
      Deallocate(  oDZ_T  )
      Deallocate(  X  )
      Deallocate(  X_E  )
      Deallocate(  oX  )
      Deallocate(  oX_E  )
      Deallocate(  Z  )
      Deallocate(  Z_T  )
      Deallocate(  FX  )
      Deallocate(  FX_bar  )
      Deallocate(  FX_E  )
      Deallocate(  FX_E_bar  )
      Deallocate(  FY_N  )
      Deallocate(  FY_N_bar  )
      Deallocate(  FZ_T  )
      Deallocate(  FZ_T_bar  )
      Deallocate(  AYZ  )
      Deallocate(  AXZ  )
      Deallocate(  AXY  )
      Deallocate(  VOL  )
      Deallocate(  AYZ_U  )
      Deallocate(  AXZ_U  )
      Deallocate(  AXY_U  )
      Deallocate(  VOL_U  )
      Deallocate(  AYZ_V  )
      Deallocate(  AXZ_V  )
      Deallocate(  AXY_V  )
      Deallocate(  VOL_V  )
      Deallocate(  AYZ_W  )
      Deallocate(  AXZ_W  )
      Deallocate(  AXY_W  )
      Deallocate(  VOL_W  )

!indices
      Deallocate(           STORE_LM  )
      Deallocate(           CELL_CLASS  )
      Deallocate(           I_OF  )
      Deallocate(           J_OF  )
      Deallocate(           K_OF  )
      Deallocate(           Im1  )
      Deallocate(           Ip1  )
      Deallocate(           Jm1  )
      Deallocate(           Jp1  )
      Deallocate(           Km1  )
      Deallocate(           Kp1  )
      


!pgcor
      Deallocate(  d_e )
      Deallocate(  d_n )
      Deallocate(  d_t )
      Deallocate(  Pp_g )
      Deallocate(  PHASE_4_P_g )

!physprop
      Deallocate(  MU_g  )
      Deallocate(  C_pg  )
      Deallocate(  C_ps  )
      Deallocate(  K_g  )
      Deallocate(  K_s  )
      Deallocate(  Kth_s  )
      Deallocate(  Kphi_s  )
      Deallocate(  DIF_g  )
      Deallocate(  DIF_s  )
      Deallocate(  MW_MIX_g  )

!pscor
      Deallocate(  e_e )
      Deallocate(  e_n )
      Deallocate(  e_t )
      Deallocate(  K_cp )
      Deallocate(  EPp )
      Deallocate(  PHASE_4_P_s )

!residual
      Deallocate( RESID )
      Deallocate( MAX_RESID )
      Deallocate( IJK_RESID )
 
!rxns
      Deallocate(  R_gp  )
      Deallocate(  R_sp  )
      Deallocate(  RoX_gc  )
      Deallocate(  RoX_sc  )
      Deallocate(  SUM_R_g  )
      Deallocate(  SUM_R_s  )
      Deallocate(  R_phase  )
      Deallocate(  MW_all  )
      Deallocate(  SPECIES_ID2N )
      Deallocate(  SPECIES_N2IDg )
      Deallocate(  SPECIES_N2IDs )
      
!scalars
      
      IF(DIMENSION_Scalar /= 0)then
        Deallocate(  Scalar_c  )
        Deallocate(  Scalar_p  )
        Deallocate(  Dif_Scalar  )
      
      ENDIF

!tau_g
      Deallocate(  TAU_U_g )
      Deallocate(  TAU_V_g )
      Deallocate(  TAU_W_g )

!tau_s
      Deallocate(  TAU_U_s )
      Deallocate(  TAU_V_s )
      Deallocate(  TAU_W_s )
      
!tmp_array
      Deallocate(  Array1 )
      Deallocate(  Array2 )
      Deallocate(  Array3 )
      Deallocate(  Array4 )
      Deallocate(  Array1i )
      Deallocate(  Array1c )

!tmp_array1
      Deallocate(  Arraym1 )

!trace
      Deallocate(  trD_s_C  )
      Deallocate(  trD_s2  )
      Deallocate(  trD_s_Co  )

!visc_g
      Deallocate(  trD_g )
      Deallocate(  MU_gt  )
      Deallocate(  LAMBDA_gt  )
      Deallocate(  L_scale  )

!visc_s
      Deallocate(  trD_s )
      Deallocate(  MU_s  )
      Deallocate(  LAMBDA_s  )
      Deallocate(  ALPHA_s  )
      Deallocate(  MU_s_c  )
      Deallocate(  LAMBDA_s_c  )
      
!xsi_array
      Deallocate(  Xsi_e )
      Deallocate(  Xsi_n )
      Deallocate(  Xsi_t )
      

!VSH
      Deallocate(  VSH )

!VSHE
      Deallocate(  VSHE )

!
! array allocation of add on packages, such as linear equation solvers
!

     
      RETURN
      END SUBROUTINE Deallocate_ARRAYS 
      
