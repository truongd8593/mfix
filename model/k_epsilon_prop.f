!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: K_Epsilon_PROP(IER)                                    C
!  Purpose: Calculate diffusion coefficeint and sources for K &        C
!           Epsilon equations                                          C
!                                                                      C
!  Author:                                                    Date:    C
!  Reviewer:                                                           C
!  Modified: S. Benyahia                         Date:May-13-04        C
!                                                                      C
!                                                                      C
!  Literature/Document References: Wilcox, D.C., Turbulence Modeling   C
!  for CFD. DCW Industries, Inc. La Canada, Ca. 1994.                  C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE K_Epsilon_PROP( IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1
      USE parallel
      USE physprop
      USE drag
      USE run
      USE output
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE trace
      USE indices
      USE constant
      Use vshear
      USE turb
      USE toleranc 
      USE compar
      USE TAU_G
      USE sendrecv
      USE time_cpu 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
        DIMENSION USX1(IMAX, JMAX), USX2(IMAX), USTIME(IMAX)
        DIMENSION USX3(IMAX), USX4(IMAX), D_g(3,3)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                       Error index
      INTEGER          IER

      INTEGER          L
!                 
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2./3. 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
                      IM, JM, KM
      INTEGER          IMJPK, IMJMK, IMJKP, IMJKM, IPJKM, IPJMK, IJMKP, &
                      IJMKM, IJPKM, II, I1, J1, K1
!		    
!                      Solids phase
      INTEGER          M  
!
!
!                      Strain rate tensor components for mth solids phase
      DOUBLE PRECISION D_g
!
!                      U_g at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION U_g_N
!
!                      U_g at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION U_g_S
!
!                      U_g at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION U_g_T
!
!                      U_g at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION U_g_B
!
!                      U_g at the center of the THETA cell-(i, j, k)
!                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION U_g_C
!
!                      V_g at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION V_g_E
!
!                      V_g at the west face of the THETA cell-(i-1/2, j, k)
      DOUBLE PRECISION V_g_W
!
!                      V_g at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION V_g_T
!
!                      V_g at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION V_g_B
!
!                      W_g at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION W_g_E
!
!                      W_g at the west face of the THETA cell-(1-1/2, j, k)
      DOUBLE PRECISION W_g_W
!
!                      W_g at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION W_g_N
!
!                      W_g at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION W_g_S
!
!                      W_g at the center of the THETA cell-(i, j, k).
!                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION W_g_C
!
!                      Second invariant of the deviator of D_g
      DOUBLE PRECISION I2_devD_g
! 
!                      Cell center value of U_g 
      DOUBLE PRECISION UGC 
! 
!                      Cell center value of V_g 
      DOUBLE PRECISION VGC 
! 
!                      Cell center value of W_g 
      DOUBLE PRECISION WGC 
!
! Use Information from calc_mu_s.f to compute the strain rate tensor
!
!
!                      Strain rate tensor components for mth solids phase
      DOUBLE PRECISION D_s(3,3)
!
!                      U_s at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION U_s_N
!
!                      U_s at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION U_s_S
!
!                      U_s at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION U_s_T
!
!                      U_s at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION U_s_B
!
!                      U_s at the center of the THETA cell-(i, j, k)
!                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION U_s_C
!
!                      V_s at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION V_s_E
!
!                      V_s at the west face of the THETA cell-(i-1/2, j, k)
      DOUBLE PRECISION V_s_W
!
!                      V_s at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION V_s_T
!
!                      V_s at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION V_s_B
!
!                      W_s at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION W_s_E
!
!                      W_s at the west face of the THETA cell-(1-1/2, j, k)
      DOUBLE PRECISION W_s_W
!
!                      W_s at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION W_s_N
!
!                      W_s at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION W_s_S
!
!                      W_s at the center of the THETA cell-(i, j, k).
!                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION W_s_C
!                      Cell center value of U_sm 
      DOUBLE PRECISION USCM 
! 
!                      Cell center value of U_g 
      DOUBLE PRECISION UGC 
! 
!                      Cell center value of V_sm 
      DOUBLE PRECISION VSCM 
! 
!                      Cell center value of V_g 
      DOUBLE PRECISION VGC 
! 
!                      Cell center value of W_sm 
      DOUBLE PRECISION WSCM 
! 
!                      Cell center value of W_g 
      DOUBLE PRECISION WGC 
! 
!                      Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION VREL 
!
!
!                      Local DO-LOOP counters
      INTEGER          I2
!
!                      Second invariant of the deviator of D_s
      DOUBLE PRECISION I2_devD_s, SRT, Trace_G, Trace_S, Mu_gas_t
! 
!
! 
!
      DOUBLE PRECISION D_12 
! 
      DOUBLE PRECISION Tau_12 , T_SOF
      DOUBLE PRECISION Tau_1 
      DOUBLE PRECISION C_beta, Cos_Theta
      DOUBLE PRECISION Zeta_r, Etha_12, SIGMA_12, B, CV,V_12V_dr,C_Eps_3
      DOUBLE PRECISION Diss , C_Eps_Pope, Xsi_Pope
      DOUBLE PRECISION K_1, K_2 
      DOUBLE PRECISION Tau_12_st, Tau_12_st_2, Re_s, CD_s
      DOUBLE PRECISION C_eps 
      DOUBLE PRECISION VREL_new
!  Production of Turb. Due to shear, Turb Visc, and Constants.
!  See Wilcox PP. 89
      DOUBLE PRECISION Tauij_gDUi_gODxj, C_MU, Sigma_k, Sigma_E, Kappa
      DOUBLE PRECISION Pos_Tauij_gDUi_gODxj, Neg_Tauij_gDUi_gODxj
      DOUBLE PRECISION Pos_Tauij_sDUi_sODxj, Neg_Tauij_sDUi_sODxj
      DOUBLE PRECISION Ceps_1, Ceps_2,PI_kq_1,PI_kq_2 
      DOUBLE PRECISION Pos_PI_kq_2, Neg_PI_kq_2
      DOUBLE PRECISION MU_12_T, MU_2_T, SIGMA_c, Tau_2_c, X_21
      DOUBLE PRECISION Tauij_sDUi_sODxj
      DOUBLE PRECISION Diss_2, Zeta_c
      DOUBLE PRECISION K_2_T, PI_q_2, Pos_PI_q_2, Neg_PI_q_2
! Modif. for Sof Local Var.      
      DOUBLE PRECISION USX1, USX3, USX4
      DOUBLE PRECISION USX2, USTIME, XKSTEP, RES_TIME     
!
      DOUBLE PRECISION S_1_dev(3,3), UG(3,3),S_12(3,3), trace_sg
      DOUBLE PRECISION S_12_dev(3,3), M_12(3,3), Us(3,3)
      
      DOUBLE PRECISION U_11_U_21,U_12_U_21,U_13_U_21
      DOUBLE PRECISION U_11_U_22,U_12_U_22,U_13_U_22
      DOUBLE PRECISION U_11_U_23,U_12_U_23,U_13_U_23
      
      DOUBLE PRECISION Tau_gi_sj, Pos_Tau_gi_sj, Neg_Tau_gi_sj
      DOUBLE PRECISION PI_q_12, Pos_PI_q_12, Neg_PI_q_12, A1, A2, A3,Check_Log
      
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
!
      IF( .NOT. K_Epsilon) RETURN
!
! M should be forced to be equal to one to get some info. from solids phase.
        M = ONE
!
! Add constants. Some of them may not be used now but will be necessary when the Simonin
! model is implemented in the standard MFIX model.
! Most of these constants have the same names and values as the ones defined in Wilcox
! book (turbulence modeling for CFD)
        C_MU = 9.0D-02
        Kappa = 0.42D+0
        Sigma_k = 1.0D0
        Sigma_E = 1.3D0
        Ceps_1 = 1.44D0 !should be 1.6 for axisymmetric cases with no Pope Correction.
        Ceps_2 = 1.92D0
        SIGMA_12 = 0.67D0
        CV = 0.5D0 !Added mass coef.
        C_Eps_3 = 1.2D0 ! for Simonin model
        C_Eps_Pope = 0.79D0
! C_Beta as defined by Eq. 20-4-77 from Fluent6 manual with Cos(Theta)=1
        C_Beta=  0.45D0	

!
! loezos
      IF (SHEAR) THEN           
      SRT=(2d0*V_sh/XLENGTH)
!$omp parallel do private(IJK)
       Do IJK= ijkstart3, ijkend3         
          IF (FLUID_AT(IJK)) THEN  
           V_s(ijk,m)=V_s(IJK,m)+VSH(IJK)
          END IF
       END DO
      END IF
! loezo   

!	
!  ---  Remember to include all the local variables here for parallel
!  ---- processing
!$omp  parallel do private(ijk, L)
      DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN
!	  
!	
!------------------
! Part copied from calc_mu_g.f to calculate the rate of strain Sii tensor
! of the gas phase.
! Sof@fluent.com
!------------------
!
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IM = IM1(I) 
            JM = JM1(J) 
            KM = KM1(K) 

            IMJK = IM_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKM = KM_OF(IJK)
            IJKP = KP_OF(IJK) 
            IMJPK = IM_OF(IJPK) 
            IMJMK = IM_OF(IJMK) 
            IMJKP = IM_OF(IJKP) 
            IMJKM = IM_OF(IJKM) 
            IPJKM = IP_OF(IJKM) 
            IPJMK = IP_OF(IJMK) 
            IJMKP = JM_OF(IJKP) 
            IJMKM = JM_OF(IJKM) 
            IJPKM = JP_OF(IJKM) 
!
!         Find fluid velocity values at faces of the cell
            U_G_N = AVG_Y(AVG_X_E(U_G(IMJK),U_G(IJK),I),AVG_X_E(U_G(IMJPK),U_G(&
               IJPK),I),J)                       !i, j+1/2, k 
            U_G_S = AVG_Y(AVG_X_E(U_G(IMJMK),U_G(IJMK),I),AVG_X_E(U_G(IMJK),U_G&
               (IJK),I),JM)                      !i, j-1/2, k 
            U_G_T = AVG_Z(AVG_X_E(U_G(IMJK),U_G(IJK),I),AVG_X_E(U_G(IMJKP),U_G(&
               IJKP),I),K)                       !i, j, k+1/2 
            U_G_B = AVG_Z(AVG_X_E(U_G(IMJKM),U_G(IJKM),I),AVG_X_E(U_G(IMJK),U_G&
               (IJK),I),KM)                      !i, j, k-1/2            
	    V_G_E = AVG_X(AVG_Y_N(V_G(IJMK),V_G(IJK)),AVG_Y_N(V_G(IPJMK),V_G(&
               IPJK)),I)                         !i+1/2, j, k 
            V_G_W = AVG_X(AVG_Y_N(V_G(IMJMK),V_G(IMJK)),AVG_Y_N(V_G(IJMK),V_G(&
               IJK)),IM)                         !i-1/2, j, k 
            V_G_T = AVG_Z(AVG_Y_N(V_G(IJMK),V_G(IJK)),AVG_Y_N(V_G(IJMKP),V_G(&
               IJKP)),K)                         !i, j, k+1/2 
            V_G_B = AVG_Z(AVG_Y_N(V_G(IJMKM),V_G(IJKM)),AVG_Y_N(V_G(IJMK),V_G(&
               IJK)),KM)                         !i, j, k-1/2 
            W_G_N = AVG_Y(AVG_Z_T(W_G(IJKM),W_G(IJK)),AVG_Z_T(W_G(IJPKM),W_G(&
               IJPK)),J)                         !i, j+1/2, k 
            W_G_S = AVG_Y(AVG_Z_T(W_G(IJMKM),W_G(IJMK)),AVG_Z_T(W_G(IJKM),W_G(&
               IJK)),JM)                         !i, j-1/2, k 
            W_G_E = AVG_X(AVG_Z_T(W_G(IJKM),W_G(IJK)),AVG_Z_T(W_G(IPJKM),W_G(&
               IPJK)),I)                         !i+1/2, j, k 
            W_G_W = AVG_X(AVG_Z_T(W_G(IMJKM),W_G(IMJK)),AVG_Z_T(W_G(IJKM),W_G(&
               IJK)),IM)                         !i-1/2, j, k 
!
            IF (CYLINDRICAL) THEN 
!                                                !i, j, k
               U_G_C = AVG_X_E(U_G(IMJK),U_G(IJK),I) 
!                                                !i, j, k
               W_G_C = AVG_Z_T(W_G(IJKM),W_G(IJK)) 
            ELSE 
               U_G_C = ZERO 
               W_G_C = ZERO 
            ENDIF 
!
!         Calculate velocity components at i, j, k to be used in wall functions
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I) 
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK)) 
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))
	     
!
!	Velocity derivives computed for the Pope Approximation.
!
          UG(1,1) = (U_G(IJK)-U_G(IMJK))*ODX(I)
          UG(1,2) = (U_G_N - U_G_S)*ODY(J)
          UG(1,3) = (U_G_T-U_G_B)*(OX(I)*ODZ(K))-W_G_C*OX(I)
          UG(2,1) = (V_G_E-V_G_W)*ODX(I)
          UG(2,2) = (V_G(IJK)-V_G(IJMK))*ODY(J)
          UG(2,3) = (V_G_T-V_G_B)*(OX(I)*ODZ(K))
          UG(3,1) = (W_G_E - W_G_W)*ODX(I)
          UG(3,2) = (W_G_N-W_G_S)*ODY(J)
          UG(3,3) = (W_G(IJK)-W_G(IJKM))*(OX(I)*ODZ(K)) + U_G_C*OX(I)
!
!
!         Find components of fluid phase strain rate
!         tensor, D_g, at center of the cell - (i, j, k)
            D_G(1,1) = (U_G(IJK)-U_G(IMJK))*ODX(I) 
            D_G(1,2) = HALF*((U_G_N - U_G_S)*ODY(J)+(V_G_E-V_G_W)*ODX(I)) 
            D_G(1,3) = HALF*((W_G_E - W_G_W)*ODX(I)+(U_G_T-U_G_B)*(OX(I)*ODZ(K)&
               )-W_G_C*OX(I)) 
            D_G(2,1) = D_G(1,2) 
            D_G(2,2) = (V_G(IJK)-V_G(IJMK))*ODY(J) 
            D_G(2,3)=HALF*((V_G_T-V_G_B)*(OX(I)*ODZ(K))+(W_G_N-W_G_S)*ODY(J)) 
            D_G(3,1) = D_G(1,3) 
            D_G(3,2) = D_G(2,3) 
            D_G(3,3) = (W_G(IJK)-W_G(IJKM))*(OX(I)*ODZ(K)) + U_G_C*OX(I) 
!!!
            Trace_g = D_G(1,1) + D_G(2,2) + D_G(3,3)

!
! Sof@fluent.com
!
! Pope's correction in 2-D to the Epsilon Equation for a round-Jet
! from Wilcox book, Page 103.
          Xsi_Pope = ZERO
          DO I1 = 1,3
            DO J1 = 1,3
              DO K1 = 1,3
                Xsi_Pope = Xsi_Pope + (UG(I1,J1) - UG(J1,I1))&
                         *(UG(J1,K1) - UG(K1,J1))*(UG(K1,I1) + UG(I1,K1))
              END DO
            END DO
          END DO
!

          Xsi_Pope = Xsi_Pope/6. * K_Turb_G(IJK)**2/E_Turb_G(IJK) 
!
        X_21 = Ep_s(IJK,M)*RO_s(M)/(EP_g(IJK)*RO_g(IJK))
!
! This IF statment is to ensure that we are using the turbulent viscosity
! and NOT the effective viscosity.
!	
        IF (MU_GT(IJK) .GE. MU_g(IJK)) THEN
          Mu_gas_t = MU_GT(IJK) - MU_g(IJK)
        ELSE
          Mu_gas_t = ZERO
        ENDIF  
!
!        Calculate Tau(i,j)*dUi/dXj (production term in the K Equation
   
         Tauij_gDUi_gODxj = 2D0*Mu_gas_t*(                             &
                            D_G(1,1) * D_G(1,1) +                      &
                            D_G(1,2) * (U_G_N - U_G_S)*ODY(J) +        &
                            D_G(1,3) * ((U_G_T-U_G_B)*                 &
                            (OX(I)*ODZ(K))-W_G_C*OX(I)) +              &
                            D_G(2,1) * (V_G_E-V_G_W)*ODX(I) +          &
                            D_G(2,2) * D_G(2,2) +                      &
                            D_G(2,3) * (V_G_T-V_G_B)*(OX(I)*ODZ(K)) +  &
                            D_G(3,1) * (W_G_E - W_G_W)*ODX(I) +        &
                            D_G(3,2) * (W_G_N-W_G_S)*ODY(J) +          &
                            D_G(3,3) * D_G(3,3)) -                     &
                            F2O3 * RO_G(IJK) * K_Turb_G(IJK)*Trace_g   &
                           - F2O3 * Mu_gas_t * Trace_g**2
! To avoid very small negative numbers
!

        IF (Tauij_gDUi_gODxj .GE. ZERO) THEN
          Pos_Tauij_gDUi_gODxj = Tauij_gDUi_gODxj
          Neg_Tauij_gDUi_gODxj = ZERO
        ELSE
          Pos_Tauij_gDUi_gODxj = ZERO
          Neg_Tauij_gDUi_gODxj = Tauij_gDUi_gODxj
        ENDIF
!
! These are the components of the fluid Reynolds stress to be used in
! the computation of the cross-correlation Reynolds stresses for Simonin model	
!
        S_1_dev(1,1) = D_G(1,1) - Trace_g /3.d0
  
        S_1_dev(1,2) = D_G(1,2)
	
        S_1_dev(2,1) = S_1_dev(1,2)
        
	S_1_dev(1,3) = D_G(1,3)
	
	S_1_dev(3,1) = S_1_dev(1,3)
	
        S_1_dev(2,2) = D_G(2,2) - Trace_g /3.d0   
	
	S_1_dev(2,3) = D_G(2,3)
	
	S_1_dev(3,2) = S_1_dev(2,3)
	
        S_1_dev(3,3) = D_G(3,3) - Trace_g /3.d0
!!!
!------------------
! Part copied from calc_mu_s.f to calculate the rate of strain Sii tensor
! of the Mth solids phase.
! Sof@fluent.com
!------------------
!  
          U_s_N = AVG_Y(                                   &   !i, j+1/2, k
                  AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I) ,&
                  AVG_X_E(U_s(IMJPK, M), U_s(IJPK, M), I) , J&
                 )
 
          U_s_S = AVG_Y(                                   &   !i, j-1/2, k
                   AVG_X_E(U_s(IMJMK, M), U_s(IJMK, M), I),&
                   AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I), JM&
                 )
          U_s_T = AVG_Z(                                 &     !i, j, k+1/2
                   AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I),&
                   AVG_X_E(U_s(IMJKP, M), U_s(IJKP, M), I), K&
                 )
          U_s_B = AVG_Z(                                 &     !i, j, k-1/2
                   AVG_X_E(U_s(IMJKM, M), U_s(IJKM, M), I),&
                   AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I), KM&
                 )
! start loezos
        IF (SHEAR)  THEN
!
          V_s_E = AVG_X(                                 &     !i+1/2, j, k
                   AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),&
                   AVG_Y_N((V_s(IPJMK, M)-VSH(IPJMK)+VSH(IJMK)&
                   +SRT*1d0/oDX_E(I)),&
                   (V_s(IPJK, M)-VSH(IPJK)+VSH(IJK)&
                   +SRT*1d0/oDX_E(I))), I&
                 )
          V_s_W = AVG_X(                                 &     !i-1/2, j, k
                   AVG_Y_N((V_s(IMJMK, M)-VSH(IMJMK)+VSH(IJMK)&
                   -SRT*1d0/oDX_E(IM1(I))),&
                   (V_s(IMJK, M)-VSH(IMJK)+VSH(IJK)&
                   -SRT*1d0/oDX_E(IM1(I)))),&
                   AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), IM)

        ELSE   

          V_s_E = AVG_X(                                 &     !i+1/2, j, k
                   AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),&
                   AVG_Y_N(V_s(IPJMK, M), V_s(IPJK, M)), I&
                 )
          V_s_W = AVG_X(                                 &     !i-1/2, j, k
                   AVG_Y_N(V_s(IMJMK, M), V_s(IMJK, M)),&
                   AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), IM&
                 )
        ENDIF ! end loezos

          V_s_T = AVG_Z(                                 &     !i, j, k+1/2
                   AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),&
                   AVG_Y_N(V_s(IJMKP, M), V_s(IJKP, M)), K&
                 )
          V_s_B = AVG_Z(                                 &    !i, j, k-1/2
                   AVG_Y_N(V_s(IJMKM, M), V_s(IJKM, M)),&
                   AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), KM&
                 )
!  
   
          W_s_N = AVG_Y(                                 &    !i, j+1/2, k
                   AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)),&
                   AVG_Z_T(W_s(IJPKM, M), W_s(IJPK, M)), J&
                 )
          W_s_S = AVG_Y(                                 &   !i, j-1/2, k
                   AVG_Z_T(W_s(IJMKM, M), W_s(IJMK, M)),&
                   AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)), JM&
                 )
          
          W_s_E = AVG_X(                                 &   !i+1/2, j, k
                   AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)),&
                   AVG_Z_T(W_s(IPJKM, M), W_s(IPJK, M)), I&
                 )
          W_s_W = AVG_X(                                 &   !i-1/2, j, k
                   AVG_Z_T(W_s(IMJKM, M), W_s(IMJK, M)),&
                   AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)), IM&
                 )
!
          IF(CYLINDRICAL) THEN
            U_s_C = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)    !i, j, k
    
            W_s_C = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))    !i, j, k
    
    
          ELSE
            U_s_C = ZERO
            W_s_C = ZERO
          ENDIF
!
!
          Us(1,2) = (U_s_N - U_s_S)*ODY(J)
          Us(1,3) = (U_s_T-U_s_B)*(OX(I)*ODZ(K))-W_s_C*OX(I)
          Us(2,1) = (V_s_E-V_s_W)*ODX(I)
          Us(2,3) = (V_s_T-V_s_B)*(OX(I)*ODZ(K))
          Us(3,1) = (W_s_E - W_s_W)*ODX(I)
          Us(3,2) = (W_s_N-W_s_S)*ODY(J)
!
!
!         Find components of Mth solids phase continuum strain rate
!         tensor, D_s, at center of THETA cell-(i, j, k)
          D_s(1,1) = ( U_s(IJK,M) - U_s(IMJK,M) ) * oDX(I)         !du1odx1
          
  
          D_s(1,2) = HALF * ( (U_s_N - U_s_S) * oDY(J) +&          !du1dx2
                               (V_s_E - V_s_W) * oDX(I) )          !du2odx1
  
          D_s(1,3) = HALF * ( (W_s_E - W_s_W) * oDX(I) +&          !du3odx1
                               (U_s_T - U_s_B) * (oX(I)*oDZ(K)) -& !du1odx3
                                 W_s_C * oX(I) )
          
          D_s(2,1) = D_s(1,2)
          D_s(2,2) = ( V_s(IJK,M) - V_s(IJMK,M) ) * oDY(J)         !du2odx2
          
          D_s(2,3) = HALF * ( (V_s_T - V_s_B) * (oX(I)*oDZ(K)) +&  !du2odx3
                               (W_s_N - W_s_S) * oDY(J) )          !du3odx2
         
          D_s(3,1) = D_s(1,3)
          D_s(3,2) = D_s(2,3)
          
          D_s(3,3) = ( W_s(IJK,M) - W_s(IJKM,M) ) * (oX(I)*oDZ(K)) +& !du3odx3
                      U_s_C * oX(I)
		       
! trace of solids strain rate tensor
          trace_s = D_s(1,1)+ D_s(2,2)+  D_s(3,3) 
! The Simonin model Ahmadi model etc. will be implemented in the near future.		       
!
! Interaction terms in the K-Epsilon equations, are set to zero for now.
          PI_kq_1 = ZERO
          Pos_PI_kq_2 = ZERO
          Neg_PI_kq_2 = ZERO

        IF(K_Turb_G(IJK) > Small_number .AND. E_Turb_G(IJK) > Small_number) THEN
!
! Start Adding source terms to the K equation
!
             K_Turb_G_c (IJK) = (EP_g(IJK) * Pos_Tauij_gDUi_gODxj +   &
                                  Pos_PI_kq_2 )
!
             K_Turb_G_p (IJK) =(-EP_g(IJK) * Neg_Tauij_gDUi_gODxj+    &
                                 EP_g(IJK)*RO_G(IJK)*E_Turb_G(IJK)+   &
                                 Neg_PI_kq_2)/K_Turb_G(IJK)
!
! Implementing wall functions targeted to fluid cells next to walls...
! Setting Source and sink terms in the K equation since the production 
! in the K eq. due to shear should include the LOG law of the wall.
!
          IF(WALL_AT(JP_OF(IJK)).OR.WALL_AT(JM_OF(IJK))) THEN
!
            Check_Log = LOG(9.81*C_mu**0.25*                          &
                            RO_G(IJK)*SQRT(K_Turb_G(IJK))*DY(J)/      &
                            2.0D+0/Mu_g(IJK))
            IF(Check_Log .LE. ZERO) THEN
	      K_Turb_G_c (IJK) = zero
	      K_Turb_G_p (IJK) = zero
	    ELSE
              K_Turb_G_c (IJK) = SQRT(C_mu)*2.D+0/DY(J)*               &
              MAX(ABS(UGC),ABS(WGC))*                        &
	      EP_g(IJK)*RO_G(IJK)*K_Turb_G(IJK)                        &
                            /Check_Log !+ Pos_PI_kq_2
 
              K_Turb_G_p (IJK) = ((EP_g(IJK)*RO_G(IJK)                 &
                *C_mu**0.75*K_Turb_G(IJK)**1.5/DY(J)*2.0D+0/Kappa)     &
                 )/K_Turb_G(IJK)
!(PI_kq_1-Neg_PI_kq_2)
	    ENDIF !for check_log less than zero
	  
          ELSE IF(WALL_AT(KP_OF(IJK)).OR.WALL_AT(KM_OF(IJK))) THEN
!!
            Check_Log = LOG(9.81*C_mu**0.25*                            &
                        RO_G(IJK)*SQRT(K_Turb_G(IJK))/                  &
                        (ODZ(K)*OX(I)*2.0D+0)/Mu_g(IJK))
            IF(Check_Log .LE. ZERO) THEN
	      K_Turb_G_c (IJK) = zero
	      K_Turb_G_p (IJK) = zero
	    ELSE 
              K_Turb_G_c (IJK) = SQRT(C_mu)*(ODZ(K)*OX(I)*2.0D+0)*      &
               MAX(ABS(UGC),ABS(VGC))                         &
	       *EP_g(IJK)*RO_G(IJK)*K_Turb_G(IJK)                       &
                        /Check_Log !+ Pos_PI_kq_2

              K_Turb_G_p (IJK) = ((EP_g(IJK)*RO_G(IJK)                  &
                 *C_mu**0.75*K_Turb_G(IJK)**1.5/Kappa*                  &
                 (ODZ(K)*OX(I)*2.0D+0))                                 &
                 )/K_Turb_G(IJK) 
!(PI_kq_1-Neg_PI_kq_2)
	    ENDIF !for check_log less than zero
	  
          ENDIF !For walls			    
!
!				
! For Cylindrical cases, wall_at (IP) is a wall cell, but wall_at (IM) is
! the axis of symmetry where wall functions obviously don't apply.
!          
          IF(CYLINDRICAL) THEN
	    IF (WALL_AT(IP_OF(IJK)))  THEN
!!
              Check_Log = LOG(9.81*C_mu**0.25*                          &
                     RO_G(IJK)*SQRT(K_Turb_G(IJK))*DX(I)/               &
                       2.0D+0/Mu_g(IJK))
		               
	      IF(Check_Log .LE. ZERO) THEN
	        K_Turb_G_c (IJK) = zero
	        K_Turb_G_p (IJK) = zero
	      ELSE		       
                K_Turb_G_c (IJK) = SQRT(C_mu)*2.D+0/DX(I)*              &
                            MAX(ABS(VGC),ABS(WGC))            &
	                   *EP_g(IJK)*RO_G(IJK)*K_Turb_G(IJK)           &
                          /Check_Log !+ Pos_PI_kq_2
 
                K_Turb_G_p (IJK) = ((EP_g(IJK)*RO_G(IJK)                &
                  *C_mu**0.75*K_Turb_G(IJK)**1.5/DX(I)*2.0D+0/Kappa)    &
		  )/K_Turb_G(IJK)
		
              ENDIF !for check_log less than zero	    
!(PI_kq_1-Neg_PI_kq_2)          
	    ENDIF  ! for wall cells in I direction
          ELSE IF (WALL_AT(IP_OF(IJK)).OR.WALL_AT(IM_OF(IJK))) THEN
!!
            Check_Log = LOG(9.81*C_mu**0.25*                            &
                     RO_G(IJK)*SQRT(K_Turb_G(IJK))*DX(I)/               &
                       2.0D+0/Mu_g(IJK))
		               
	    IF(Check_Log .LE. ZERO) THEN
	      K_Turb_G_c (IJK) = zero
	      K_Turb_G_p (IJK) = zero
	    ELSE	
              K_Turb_G_c (IJK) = SQRT(C_mu)*2.D+0/DX(I)*                &
               MAX(ABS(VGC),ABS(WGC))*                        &
	       EP_g(IJK)*RO_G(IJK)*K_Turb_G(IJK)                        &
                              /Check_Log  !+ Pos_PI_kq_2

              K_Turb_G_p (IJK) = ((EP_g(IJK)*RO_G(IJK)                  &
                  *C_mu**0.75*K_Turb_G(IJK)**1.5/DX(I)*2.0D+0/Kappa)    &
		)/K_Turb_G(IJK)
!(PI_kq_1-Neg_PI_kq_2)
            ENDIF !for check_log less than zero
          ENDIF ! for cylindrical  

!	 
!            Diffusion coefficient for turbulent kinetic energy (K)
!            
             Dif_K_Turb_G(IJK) = EP_g(IJK)* (MU_G(IJK) + Mu_gas_t /Sigma_k)
!
! Add here Dissipation of Turbulence source terms
!	     
             E_Turb_G_c (IJK) = (Ceps_1 *EP_g(IJK)*Pos_Tauij_gDUi_gODxj &
                                *E_Turb_G(IJK)/K_Turb_G(IJK) +          &
                                C_Eps_3*Pos_PI_kq_2*                    &
                                E_Turb_G(IJK)/K_Turb_G(IJK))           !&
!				
! Pope Correction in Xsi_Pope, Add it to E_Turb_G_c to use this option.
!
	                       ! + C_Eps_Pope*RO_G(IJK)*EP_g(IJK)*&
			       !  ZMAX(Xsi_Pope)
!
             E_Turb_G_p (IJK) = -Ceps_1 *EP_g(IJK)*Neg_Tauij_gDUi_gODxj &
                                 /K_Turb_G(IJK) +                       &
                                Ceps_2 * EP_g(IJK) *RO_G(IJK)           &
                                *E_Turb_G(IJK)/K_Turb_G(IJK)            &
                                + C_Eps_3*(Neg_PI_kq_2)                 &
                                 /K_Turb_G(IJK)
!
!
!            Diffusion coefficient for Dissipation of turbulent energy (Epsilon)
!            
             Dif_E_Turb_G(IJK) =EP_g(IJK)* (MU_G(IJK) + Mu_gas_t /Sigma_E)

        ELSE
	  K_Turb_G_c (IJK) = zero
	  K_Turb_G_p (IJK) = zero
	  E_Turb_G_c (IJK) = zero
	  E_Turb_G_p (IJK) = zero
	  Dif_K_Turb_G(IJK) = zero
	  Dif_E_Turb_G(IJK) = zero
	  
        ENDIF !for K_turb_G and E_Turb_G having very small numbers
	
         ENDIF 
      END DO
!
!
      RETURN  
      END SUBROUTINE K_Epsilon_PROP 
