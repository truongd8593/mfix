!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_GRANULAR_ENERGY(sourcelhs,sourcerhs,IJK,M,IER)  C
!  Purpose: Calculate the source terms in the granular energy equation C
!                                                                      C
!  Author: Kapil Agrawal, Princeton University        Date: 04-FEB-98  C
!  Reviewer: M. Syamlal                               Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!     Local variables: sourcelhs, sourcerhs                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_GRANULAR_ENERGY(SOURCELHS, SOURCERHS, IJK, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE run
      USE drag
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE trace
      USE turb
      USE indices
      USE constant
      USE toleranc
      USE compar        !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!     
!                      Error index
      INTEGER          IER
!     
!                      Indices
      INTEGER          IJK, I, J, K, M, MM
!
!                      Source terms to be kept on rhs
      DOUBLE PRECISION sourcerhs
!
!                      Source terms to be kept on lhs
      DOUBLE PRECISION sourcelhs
!
!                      Particle relaxation time
      DOUBLE PRECISION Tau_12_st
!
!                      Sum of eps*G_0
      DOUBLE PRECISION SUM_EpsGo
!
!                      Slip velocity
      DOUBLE PRECISION VSLIP
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: G_0 
!-----------------------------------------------
      INCLUDE 's_pr1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 's_pr2.inc'
!
      I = I_OF(IJK) 
      J = J_OF(IJK) 
      K = K_OF(IJK) 
!
      SOURCERHS = (ZMAX(LAMBDA_S_C(IJK,M))*TRD_S_C(IJK,M)**2d0+2d0*MU_S_C(IJK,M)*&
         TRD_S2(IJK,M)+P_S_C(IJK,M)*ZMAX((-TRD_S_C(IJK,M))))*VOL(IJK)  
!
      IF(SIMONIN) THEN
        SOURCERHS = SOURCERHS +  F_GS(IJK,M)*K_12(IJK)*VOL(IJK)
!
      ELSE IF(AHMADI) THEN
	IF(Ep_s(IJK,M) > DIL_EP_S .AND. F_GS(IJK,1) > small_number) THEN
          Tau_12_st = Ep_s(IJK,M)*RO_s(M)/F_GS(IJK,1)
	  SOURCERHS = SOURCERHS + 2.D+0*F_GS(IJK,M)* (ONE/(ONE+Tau_12_st/  &
                      (Tau_1(ijk)+small_number)))*K_Turb_G(IJK)*VOL(IJK)
	ELSE
	  SOURCERHS = SOURCERHS
	ENDIF
!
      ELSE IF(SWITCH > ZERO .AND. RO_g0 /= ZERO) THEN ! no modifications done to the original KTGF
!
        VSLIP = (U_S(IJK,M)-U_G(IJK))**2 + (V_S(IJK,M)-V_G(IJK))**2 + (W_S(IJK,M)&
           -W_G(IJK))**2 
        VSLIP = DSQRT(VSLIP) 
!
        SOURCERHS = SOURCERHS + (SWITCH*81D0*EP_S(IJK,M)*(MU_G(IJK)*VSLIP)**2D0/(&
           G_0(IJK,M,M)*D_P(IJK,M)**3D0*RO_S(M)*(PI*THETA_M(IJK,M)+SMALL_NUMBER)**&
           0.5D0))*VOL(IJK) 
      ENDIF
!
!---------------------------------------------------------------------
!  The following lines are commented out since Kphi_s has been set to zero
!  in subroutine calc_mu_s.  To activate the feature uncomment the following
!  lines and the lines in calc_mu_s.
!
!
!      sourcerhs=sourcerhs+
!     &    AVG_X_h(Kphi_s(IJK,M),Kphi_s(EAST_OF(IJK),M),I)
!     &    *(EP_s(EAST_OF(IJK),M)-EP_s(IJK,M))*oDX_E(I)
!     &    *AYZ(IJK)
!
!      sourcerhs=sourcerhs-
!     &    AVG_X_h(Kphi_s(WEST_OF(IJK),M),Kphi_s(IJK,M)
!     &    ,I_OF(WEST_OF(IJK)))*(EP_s(IJK,M)-EP_s(WEST_OF(IJK),M))
!     &    *oDX_E(I_OF(WEST_OF(IJK)))*AYZ(WEST_OF(IJK))
!
!      sourcerhs=sourcerhs+
!     &    AVG_Y_h(Kphi_s(IJK,M),Kphi_s(NORTH_OF(IJK),M),J)
!     &    *(EP_s(NORTH_OF(IJK),M)-EP_s(IJK,M))*oDY_N(J)
!     &    *AXZ(IJK)
!
!      sourcerhs=sourcerhs-
!     &    AVG_X_h(Kphi_s(SOUTH_OF(IJK),M),Kphi_s(IJK,M)
!     &    ,J_OF(SOUTH_OF(IJK)))*(EP_s(IJK,M)-EP_s(SOUTH_OF(IJK),M))
!     &    *oDY_N(J_OF(SOUTH_OF(IJK)))
!     &    *AXZ(JM_OF(IJK))
!
!      sourcerhs=sourcerhs+
!     &    AVG_Z_h(Kphi_s(IJK,M),Kphi_s(TOP_OF(IJK),M),K)
!     &    *(EP_s(TOP_OF(IJK),M)-EP_s(IJK,M))
!     &    *oX(I)*oDZ_T(K)
!     &    *AXY(IJK)
!
!      sourcerhs=sourcerhs-
!     &    AVG_Z_h(Kphi_s(BOTTOM_OF(IJK),M),Kphi_s(IJK,M)
!     &    ,K_OF(BOTTOM_OF(IJK)))*(EP_s(IJK,M)-
!     &    EP_s(BOTTOM_OF(IJK),M))
!     &    *oX(I)*oDZ_T(K_OF(BOTTOM_OF(IJK)))
!     &    *AXY(KM_OF(IJK))
!---------------------------------------------------------------------
!
! Changes needed for multitype particles, sof June 16 2005
! Sum of eps*G_0 is used instead of Eps*G_0
!
      SUM_EpsGo = ZERO
      DO MM = 1, MMAX
        SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,MM)*G_0(IJK,MM,MM)
      ENDDO
      SOURCELHS = ((48d0/DSQRT(PI))*ETA*(ONE-ETA)*ROP_S(IJK&
         ,M)*SUM_EpsGo*DSQRT(THETA_M(IJK,M))/D_P(IJK,M)+ &
	 P_S_C(IJK,M)*ZMAX((TRD_S_C(IJK,M)))/(THETA_M(IJK,M)+SMALL_NUMBER) &
	 +ZMAX((-LAMBDA_S_C(IJK,M)))*TRD_S_C(IJK,M)**2d0/(THETA_M(IJK,M)+&
         SMALL_NUMBER))*VOL(IJK)  
!
      IF(SIMONIN .OR. AHMADI) THEN
        SOURCELHS = SOURCELHS +  3d0 *F_GS(IJK,M)*VOL(IJK)
!
      ELSE IF(SWITCH > ZERO .AND. RO_g0 /= ZERO) THEN ! no modifications done to the original KTGF
!
        SOURCELHS = SOURCELHS + SWITCH *3d0 *F_GS(IJK,M)*VOL(IJK)
      ENDIF
!
      RETURN  
      END SUBROUTINE SOURCE_GRANULAR_ENERGY 
!-----------------------------------------------  
!
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                      
!  Module name: SOURCE_IA_NONEP_GRANULAR_ENERGY(SOURCELHS, SOURCERHS, IJK, M, IER)
!  Purpose: This module computes the source term of the species
!           granular energy equation
!
!  Literature/Document References:    
!	Iddir, Y.H., "Modeling of the multiphase mixture of particles 
!	  using the kinetic theory approach," PhD Thesis, Illinois
!	  Institute of Technology, Chicago, Illinois, 2004
!    Iddir, Y.H., & H. Arastoopour, "Modeling of multitype particle
!      flow using the kinetic theory approach," AIChE J., Vol 51,
!      No 6, June 2005
!                                                         
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      SUBROUTINE SOURCE_IA_NONEP_GRANULAR_ENERGY(SOURCELHS, SOURCERHS, IJK, M, IER) 
!
!-----------------------------------------------
!     Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE run
      USE drag
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE trace
      USE turb
      USE indices
      USE constant
      USE toleranc
      USE residual
      use kintheory
      USE compar        !//d
      IMPLICIT NONE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------  
!                      Error index
      INTEGER          IER
!     
!                      Indices
      INTEGER          IJK, I, J, K, IM, JM, KM, IMJK, IJMK, IJKM,&
                       IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
!
!                      phase index 
      INTEGER          M, L, LM
!
!                      velocities
      DOUBLE PRECISION UsM_e, UsM_w, VsM_n, VsM_s, WsM_t, WsM_b,&
                       UsL_e, UsL_w, VsL_n, VsL_s, WsL_t, WsL_b,&
                       UsM_p, VsM_p, WsM_P, UsL_p, VsL_p, WsL_p
!
!                      number densities
      DOUBLE PRECISION NU_PM_E, NU_PM_W, NU_PM_N, NU_PM_S, NU_PM_T,&
                       NU_PM_B, NU_PM_p,&
                       NU_PL_E, NU_PL_W, NU_PL_N, NU_PL_S, NU_PL_T,&
                       NU_PL_B, NU_PL_p
!
!                      temperature of species L
      DOUBLE PRECISION T_PL_E, T_PL_W, T_PL_N, T_PL_S, T_PL_T,&
                       T_PL_B, T_PL_p
!
!                      particle characteristics
      DOUBLE PRECISION M_PM, M_PL, D_PM, D_PL
!
!                      coefficients in heat flux term
      DOUBLE PRECISION Knu_sL_e, Knu_sL_w, Knu_sL_n, Knu_sL_s, Knu_sL_t,&
                       Knu_sL_b, Knu_sM_e, Knu_sM_w, Knu_sM_n, Knu_sM_s,&
                       Knu_sM_t, Knu_sM_b,&
                       Kvel_s_e, Kvel_s_w, Kvel_s_n, Kvel_s_s, Kvel_s_t,&
                       Kvel_s_b,&
                       Kth_sL_e, Kth_sL_w, Kth_sL_n, Kth_sL_s, Kth_sL_t,&
                       Kth_sL_b
!
!                      Source terms to be kept on rhs
      DOUBLE PRECISION sourcerhs, S10_rhs, S15_rhs, S16_rhs,&
                       S11a_sum_rhs, S11b_sum_rhs, S11c_sum_rhs, S12a_sum_rhs,&
                       S12b_sum_rhs, S12c_sum_rhs, &
                       S13_sum_lhs, S13_sum_rhs, &
                       S14a_sum_rhs, S14b_sum_rhs, S14c_sum_rhs,&
                       S17_sum_rhs, S18_sum_rhs, &
                       S9_sum_rhs, s21a_sum_rhs, S21b_sum_rhs, S21c_sum_rhs
!
!                      Source terms to be kept on lhs
      DOUBLE PRECISION sourcelhs, S10_lhs, S16_lhs,&
                       S11a_sum_lhs, S11b_sum_lhs, S11c_sum_lhs, S12a_sum_lhs,&
                       S12b_sum_lhs, S12c_sum_lhs, &
                       S17_sum_lhs, S18_sum_lhs, S20_sum_lhs

!
!----------------------------------------------- 
!     Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------
!
      I = I_OF(IJK) 
      J = J_OF(IJK) 
      K = K_OF(IJK) 
      IM = Im1(I)
      JM = Jm1(J)
      KM = Km1(K)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK) 
      IJKE = EAST_OF(IJK) 
      IJKW = WEST_OF(IJK)
      IJKN = NORTH_OF(IJK) 
      IJKS = SOUTH_OF(IJK)
      IJKT = TOP_OF(IJK) 
      IJKB = BOTTOM_OF(IJK) 
!
!     initialize summation variables
      S9_sum_rhs = ZERO
      S11a_sum_rhs = ZERO
      S12a_sum_rhs = ZERO 
      S13_sum_rhs = ZERO  
      S14a_sum_rhs = ZERO
      S14b_sum_rhs = ZERO
      S14c_sum_rhs = ZERO   
      S17_sum_rhs = ZERO
      S18_sum_rhs = ZERO
      S21a_sum_rhs = ZERO
      S21b_sum_rhs = ZERO
      S21c_sum_rhs = ZERO

      S11a_sum_lhs = ZERO
      S11b_sum_lhs = ZERO
      S11c_sum_lhs = ZERO
      S12a_sum_lhs = ZERO
      S12b_sum_lhs = ZERO
      S12c_sum_lhs = ZERO
      S13_sum_lhs = ZERO  
      S17_sum_lhs = ZERO
      S18_sum_lhs = ZERO
      S20_sum_lhs = ZERO
!
      UsM_e = U_S(IJK,M)
      UsM_w = U_S(IMJK,M)
      VsM_n = V_S(IJK,M)
      VsM_s = V_S(IJMK,M)
      WsM_t = W_S(IJK,M)
      WsM_b = W_S(IJKM,M)
      UsM_p = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
      VsM_p = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M) )
      WsM_p = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M) )
! 
      D_PM = D_P(IJK,M) 
      M_PM = (Pi/6.d0)*D_PM**3 * RO_S(M)
      NU_PM_p = ROP_S(IJK,M)/M_PM
      NU_PM_E = ROP_S(IJKE,M)/M_PM
      NU_PM_W = ROP_S(IJKW,M)/M_PM
      NU_PM_N = ROP_S(IJKN,M)/M_PM
      NU_PM_S = ROP_S(IJKS,M)/M_PM
      NU_PM_T = ROP_S(IJKT,M)/M_PM
      NU_PM_B = ROP_S(IJKB,M)/M_PM
!
!
!     Production by shear: (S:grad(vi))
!         Pi_s*tr(Di)

      S10_lhs = P_S_C(IJK,M) * ZMAX(TRD_S_C(IJK,M)) 
      S10_rhs = P_S_C(IJK,M) * ZMAX(-TRD_S_C(IJK,M))
!
!
!     Production by shear: (S:grad(vi))  
!         Mu_s*tr(Di^2)
      S15_rhs = 2.d0*Mu_s_c(IJK,M)*TRD_S2(IJK,M)
!
!
!     Production by shear: (S:grad(vi))  
!         Lambda_s*tr(Di)^2
      S16_lhs = (TRD_S_C(IJK,M)**2)*ZMAX( -LAMBDA_s_c(IJK,M) )
      S16_rhs = (TRD_S_C(IJK,M)**2)*ZMAX(  LAMBDA_s_C(IJK,M) )
	 
!
!
      DO L = 1, MMAX
          D_PL = D_P(IJK,L) 
          M_PL = (Pi/6.d0)*D_PL**3 * RO_S(L)
          NU_PL_p = ROP_S(IJK,L)/M_PL
          NU_PL_E = ROP_S(IJKE,L)/M_PL
          NU_PL_W = ROP_S(IJKW,L)/M_PL
          NU_PL_N = ROP_S(IJKN,L)/M_PL
          NU_PL_S = ROP_S(IJKS,L)/M_PL   
          NU_PL_T = ROP_S(IJKT,L)/M_PL
          NU_PL_B = ROP_S(IJKB,L)/M_PL

          T_PL_p = Theta_m(IJK,L)
          T_PL_E = Theta_m(IJKE,L)
          T_PL_W = Theta_m(IJKW,L)
          T_PL_N = Theta_m(IJKN,L)
          T_PL_S = Theta_m(IJKS,L)
          T_PL_T = Theta_m(IJKT,L)
          T_PL_B = Theta_m(IJKB,L)

          UsL_e = U_S(IJK,L)
          UsL_w = U_S(IMJK,L)
          VsL_n = V_S(IJK,L)
          VsL_s = V_S(IJMK,L)
          WsL_t = W_S(IJK,L)
          WsL_b = W_S(IJKM,L)
          UsL_p = AVG_X_E(U_S(IMJK,L),U_S(IJK,L),I)
          VsL_p = AVG_Y_N(V_S(IJMK,L),V_S(IJK,L) )
          WsL_p = AVG_Z_T(W_S(IJKM,L),W_S(IJK,L) )
!
!
!         Energy dissipation by collisions: Sum(Nip)
!              SUM( EDT_s_ip )
          S20_sum_lhs = S20_sum_lhs + EDT_s_ip(IJK,M,L)
!
!
!         Energy dissipation by collisions: SUM(Nip)
!              SUM( EDvel_sL_ip* div(vp) ) !Modified by sof to include trace of V_s_L
          
	  S11a_sum_lhs = S11a_sum_lhs + ZMAX(-EDvel_sL_ip(IJK,M,L)* &
               TRD_S_C(IJK,L) ) * VOL(IJK)
          
	  S11a_sum_rhs = S11a_sum_rhs + ZMAX( EDvel_sL_ip(IJK,M,L)* &
               TRD_S_C(IJK,L) ) * VOL(IJK)
!
!
!         Energy dissipation by collisions: Sum(Nip)
!              SUM( EDvel_sM_ip* div(vi) ) !Modified by sof to include trace of V_s_M
          
	  S12a_sum_lhs = S12a_sum_lhs + ZMAX(-EDvel_sM_ip(IJK,M,L)*&
               TRD_S_C(IJK,M) ) * VOL(IJK)
          
	  S12a_sum_rhs = S12a_sum_rhs + ZMAX( EDvel_sM_ip(IJK,M,L)*&
               TRD_S_C(IJK,M) ) * VOL(IJK)
!
!                    
          IF (M .NE. L) THEN
	        LM = FUNLM(L,M)
!
!
!              Production by shear: (S:grad(vi))  
!                   SUM(2*Mu_sL_ip*tr(Dk*Di) )
               S17_sum_lhs = S17_sum_lhs + 2.d0*MU_sL_ip(IJK,M,L)*&
                    ZMAX( - TRD_s2_ip(IJK,M,L) )
               S17_sum_rhs = S17_sum_rhs + 2.d0*MU_sL_ip(IJK,M,L)*&
                    ZMAX( TRD_s2_ip(IJK,M,L) )
!
! These two terms can be treated explicitly here by uncommenting the following
! two lines. They are currently treated by PEA algorithm in solve_granular_energy
!
!               S13_sum_lhs = S13_sum_lhs + ED_ss_ip(IJK,LM)*Theta_m(IJK,M)
!               S13_sum_rhs = S13_sum_rhs + ED_ss_ip(IJK,LM)*Theta_m(IJK,L)
!
!
!              Production by shear: (S:grad(vi))  
!                   SUM( (Xi_sL_ip-(2/3)*Mu_sL_ip)*tr(Dk)tr(Di) )
               S18_sum_lhs = S18_sum_lhs + ZMAX(-1.d0* &
                    (Xi_sL_ip(IJK,M,L)-(2.d0/3.d0)*Mu_sL_ip(IJK,M,L))*&
                    TRD_S_C(IJK,M)*TRD_S_C(IJK,L) )
               S18_sum_rhs = S18_sum_rhs + ZMAX(&
                    (Xi_sL_ip(IJK,M,L)-(2.d0/3.d0)*Mu_sL_ip(IJK,M,L))*&
                    TRD_S_C(IJK,M)*TRD_S_C(IJK,L) )
!
!
!              Part of Heat Flux: div (q)
!                   Kth_sL_ip*[grad(Tp)]
!              Note for L=M S21 terms cancel with similar term arising from
!                   grad(Ti)
               Kth_sL_e = AVG_X_S(Kth_sL_ip(IJK,M,L), Kth_sL_ip(IJKE,M,L),I)
               Kth_sL_w = AVG_X_S(Kth_sL_ip(IJKW,M,L),Kth_sL_ip(IJK,M,L), IM)
               Kth_sL_n = AVG_Y_S(Kth_sL_ip(IJK,M,L), Kth_sL_ip(IJKN,M,L),J)
               Kth_sL_s = AVG_Y_S(Kth_sL_ip(IJKS,M,L),Kth_sL_ip(IJK,M,L), JM)
               Kth_sL_t = AVG_Z_S(Kth_sL_ip(IJK,M,L), Kth_sL_ip(IJKT,M,L),K)
               Kth_sL_b = AVG_Z_S(Kth_sL_ip(IJKB,M,L),Kth_sL_ip(IJK,M,L), KM)

               S21a_sum_rhs = S21a_sum_rhs + ( (Kth_sL_e*(T_PL_E-T_PL_p) )*&
                    ODX_E(I)*AYZ(IJK) - (Kth_sL_w*(T_PL_p-T_PL_W) )*ODX_E(IM)*&
                    AYZ(IMJK) )
               S21b_sum_rhs = S21b_sum_rhs + ( (Kth_sL_n*(T_PL_N-T_PL_p) )*&
                    ODY_N(J)*AXZ(IJK) - (Kth_sL_s*(T_PL_p-T_PL_S) )*ODY_N(JM)*&
                    AXZ(IJMK) )
               S21c_sum_rhs = S21c_sum_rhs + ( (Kth_sL_t*(T_PL_T-T_PL_p) )*&
                    ODZ_T(K)*OX(I)*AXY(IJK) - (Kth_sL_b*(T_PL_p-T_PL_B) )*&
                    ODZ_T(KM)*OX(I)*AXY(IJKM) )
!
!
!              Part of Heat Flux: div (q)
!                   Knu_s_ip*[ni*grad(np)-np*grad(ni)]
!              Note S14 terms should evaluate to zero for particles from the same phase
               Knu_sL_e = AVG_X_S(Knu_sL_ip(IJK,M,L), Knu_sL_ip(IJKE,M,L),I)
               Knu_sM_e = AVG_X_S(Knu_sM_ip(IJK,M,L), Knu_sM_ip(IJKE,M,L),I)
               Knu_sL_w = AVG_X_S(Knu_sL_ip(IJKW,M,L),Knu_sL_ip(IJK,M,L), IM)
               Knu_sM_w = AVG_X_S(Knu_sM_ip(IJKW,M,L),Knu_sM_ip(IJK,M,L), IM)
               Knu_sL_n = AVG_Y_S(Knu_sL_ip(IJK,M,L), Knu_sL_ip(IJKN,M,L),J)
               Knu_sM_n = AVG_Y_S(Knu_sM_ip(IJK,M,L), Knu_sM_ip(IJKN,M,L),J)
               Knu_sL_s = AVG_Y_S(Knu_sL_ip(IJKS,M,L),Knu_sL_ip(IJK,M,L), JM)
               Knu_sM_s = AVG_Y_S(Knu_sM_ip(IJKS,M,L),Knu_sM_ip(IJK,M,L), JM)
               Knu_sL_t = AVG_Z_S(Knu_sL_ip(IJK,M,L), Knu_sL_ip(IJKT,M,L),K)
               Knu_sM_t = AVG_Z_S(Knu_sM_ip(IJK,M,L), Knu_sM_ip(IJKT,M,L),K)
               Knu_sL_b = AVG_Z_S(Knu_sL_ip(IJKB,M,L),Knu_sL_ip(IJK,M,L), KM)
               Knu_sM_b = AVG_Z_S(Knu_sM_ip(IJKB,M,L),Knu_sM_ip(IJK,M,L), KM)

               S14a_sum_rhs = S14a_sum_rhs + ( (Knu_sL_e*(NU_PL_E-NU_PL_p) - &
                    Knu_sM_e*(NU_PM_E-NU_PM_p) )*ODX_E(I)*AYZ(IJK) - (Knu_sL_w*&
                    (NU_PL_p-NU_PL_W) - Knu_sM_w*(NU_PM_p-NU_PM_W) )*ODX_E(IM)*&
                    AYZ(IMJK) )
               S14b_sum_rhs = S14b_sum_rhs + ( (Knu_sL_n*(NU_PL_N-NU_PL_p) - &
                    Knu_sM_n*(NU_PM_N-NU_PM_p) )*ODY_N(J)*AXZ(IJK) - (Knu_sL_s*&
                    (NU_PL_p-NU_PL_S) - Knu_sM_s*(NU_PM_p-NU_PM_S) )*ODY_N(JM)*&
                    AXZ(IJMK) )
               S14c_sum_rhs = S14c_sum_rhs + ( (Knu_sL_t*(NU_PL_T-NU_PL_p) - &
                    Knu_sM_t*(NU_PM_T-NU_PM_p) )*ODZ_T(K)*OX(I)*AXY(IJK) - &
                    (Knu_sL_b*(NU_PL_p-NU_PL_B) - Knu_sM_b*(NU_PM_p-NU_PM_B) )*&
                    ODZ_T(KM)*OX(I)*AXY(IJKM) )
!
!
!              Part of Heat Flux: div (q)
!                   Kvel_s_ip*[vi-vp]
!              Note S9 terms should evaluate to zero for particles from the same phase
               Kvel_s_e = AVG_X_H(Kvel_s_ip(IJK,M,L), Kvel_s_ip(IJKE,M,L),I)
               Kvel_s_w = AVG_X_H(Kvel_s_ip(IJKW,M,L),Kvel_s_ip(IJK,M,L), IM)
               Kvel_s_n = AVG_Y_H(Kvel_s_ip(IJK,M,L), Kvel_s_ip(IJKN,M,L),J)
               Kvel_s_s = AVG_Y_H(Kvel_s_ip(IJKS,M,L),Kvel_s_ip(IJK,M,L), JM)
               Kvel_s_t = AVG_Z_H(Kvel_s_ip(IJK,M,L), Kvel_s_ip(IJKT,M,L),K)
               Kvel_s_b = AVG_Z_H(Kvel_s_ip(IJKB,M,L),Kvel_s_ip(IJK,M,L), KM)
 
               S9_sum_rhs = S9_sum_rhs + ( Kvel_s_e*(UsM_e-UsL_e)*AYZ(IJK) - &
                    Kvel_s_w*(UsM_w-UsL_w)*AYZ(IMJK) + Kvel_s_n*(VsM_n-VsL_n)*AXZ(IJK)-&
                    Kvel_s_s*(VsM_s-VsL_s)*AXZ(IJMK) + Kvel_s_t*(WsM_t-WsL_t)*AXY(IJK)-&
                    Kvel_s_b*(WsM_b-WsL_b)*AXY(IJKM) )
!
!
          ENDIF    ! (IF M.NE.L)
!
      ENDDO
!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!
! The commented two lines in LHS and RHS caused convergence issues. These
! need to be revisited and debugged before release (sof, Nov 07 2006)
!
      SOURCELHS = ( (S11a_sum_lhs+S11b_sum_lhs+S11c_sum_lhs+S12a_sum_lhs+&
          S12b_sum_lhs+S12c_sum_lhs) + (S10_lhs+S16_lhs+S17_sum_lhs+&
          S18_sum_lhs-S20_sum_lhs + S13_sum_lhs)*VOL(IJK) + &
!	  ZMAX(S21a_sum_rhs+S21b_sum_rhs+S21c_sum_rhs)+ &
	  ZMAX(S14a_sum_rhs+S14b_sum_rhs+S14c_sum_rhs)+ ZMAX(S9_sum_rhs) ) / &
          Theta_m(IJK,M)
!
!
      SOURCERHS = ( S10_rhs+ S15_rhs + S16_rhs + S17_sum_rhs+S18_sum_rhs  + S13_sum_rhs) * VOL(IJK) + &
          S11a_sum_rhs+S11b_sum_rhs+S11c_sum_rhs+S12a_sum_rhs+S12b_sum_rhs+S12c_sum_rhs + &
          ZMAX(- (S14a_sum_rhs+S14b_sum_rhs+S14c_sum_rhs) ) + ZMAX(-S9_sum_rhs)! + &
!	  ZMAX(- (S21a_sum_rhs+S21b_sum_rhs+S21c_sum_rhs) )
	  
!	  
!
      RETURN  
      END SUBROUTINE SOURCE_IA_NONEP_GRANULAR_ENERGY
!-----------------------------------------------  

