!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_MU_s(M, IER)                                      C
!  Purpose: Calculate granular stress terms: THETA, P_s, LAMBDA_s, MU_sC
!                                                                      C
!  Author: W. Rogers                                  Date: 04-mar-92  C
!  Reviewer: M. Syamlal                               Date: 16-MAR-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modifications for cylindrical geometry                     C
!  Author: M. Syamlal                                 Date: 15-MAY-92  C
!  Revision Number: 2                                                  C
!  Purpose: Add volume-weighted averaging statement functions for      C
!           variable grid capability                                   C
!  Author:  W. Rogers                                 Date: 21-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!  Revision Number: 3                                                  C
!  Purpose: Add plastic-flow stress terms                              C
!  Author: M. Syamlal                                 Date: 10-FEB-93  C
!  Revision Number: 4                                                  C
!  Purpose: Add Boyle-Massoudi stress terms                            C
!  Author: M. Syamlal                                 Date: 2-NOV-95   C
!  Revision Number: 5                                                  C
!  Purpose: MFIX 2.0 mods  (old name CALC_THETA)                       C
!  Author: M. Syamlal                                 Date: 24-APR_96  C
!  Author: Kapil Agrawal, Princeton University        Date: 6-FEB-98   C
!  Revision Number: 6                                                  C
!  Purpose: Add calculation of viscosities and conductivities for use  C
!           with granular temperature PDE. New common block contained  C
!           in 'trace.inc' contains trD_s_C(DIMENSION_3, DIMENSION_M)  C
!           and trD_s2(DIMENSION_3, DIMENSION_M)                       C
!  Author: Anuj Srivastava, Princeton University      Date: 20-APR-98  C
!  Revision Number:7                                                   C
!  Purpose: Add calculation of frictional stress terms                 C
!  Author: Sofiane Benyahia, Fluent, Inc.				C
!  Revision Number:8							C
!  Purpose: Add Turbulent Viscosity from Simonin Model			C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: U_s, V_s, W_s, IMAX2, JMAX2, KMAX2, DX, DY,   C
!                        DZ, IMJPK, IMJK, IPJMK, IPJK, IJMK, IJKP,     C
!                        IMJKP, IPJKM, IJKM, IJMKP, IJPK, IJPKM, IJMK, C
!                        M,  RO_s, C_e, D_p, Pi, G_0, X                C
!                                                                      C
!  Variables modified: I, J, K, IJK, MU_s, LAMBDA_s, P_s               C
!                                                                      C
!  Local variables: K_1m, K_2m, K_3m, K_4m, D_s, U_s_N, U_s_S, V_s_E,  C
!                   V_s_W, U_s_T, U_s_B, W_s_E, W_s_W, V_s_T, V_s_B,   C
!                   W_s_N, W_s_S, trD_s_C, W_s_C                       C
!                   trD_s2, EP_s2xTHETA, EP_sxSQRTHETA, I1, I2, U_s_C, C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_MU_s(M, IER)
!
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE drag
      USE run
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE trace
      USE indices
      USE constant
      Use vshear
      USE scalars
      USE compar
      USE sendrecv
      IMPLICIT NONE
!                      Maximum value of solids viscosity in poise
      DOUBLE PRECISION MAX_MU_s
      PARAMETER (MAX_MU_s = 1000.)
!
!  Function subroutines
!
      DOUBLE PRECISION G_0
!
!     Local Variables
!
!                      Error index
      INTEGER          IER
!
!                      Constant in equation for mth solids phase pressure
      DOUBLE PRECISION K_1m
!
!                      Constant in equation for mth solids phase bulk viscosity
      DOUBLE PRECISION K_2m
!
!                      Constant in equation for mth solids phase viscosity
      DOUBLE PRECISION K_3m
!
!                      Constant in equation for mth solids phase dissipation
      DOUBLE PRECISION K_4m
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
!
!                      Value of EP_s * SQRT( THETA )for Mth solids phase
!                      continuum
      DOUBLE PRECISION EP_sxSQRTHETA
!
!                      Value of EP_s * EP_s * THETA for Mth solids phase
!                      continuum
      DOUBLE PRECISION EP_s2xTHETA
!
!                      Local DO-LOOP counters
      INTEGER          I1, I2
!
!                      Second invariant of the deviator of D_s
      DOUBLE PRECISION I2_devD_s
!
!                      Factor in plastic-flow stress terms
      DOUBLE PRECISION qxP_s
!
!                      Coefficients of quadratic equation
      DOUBLE PRECISION aq, bq, cq
!
!                      Constant in Boyle-Massoudi stress term
      DOUBLE PRECISION K_5m
!
!                      d(EP_sm)/dX
      DOUBLE PRECISION DEP_soDX
!
!                      d(EP_sm)/dY
      DOUBLE PRECISION DEP_soDY
!
!                      d(EP_sm)/XdZ
      DOUBLE PRECISION DEP_soXDZ
!
!                      Solids volume fraction gradient tensor
      DOUBLE PRECISION M_s(3,3)
!
!                      Trace of M_s
      DOUBLE PRECISION trM_s
!
!                      Trace of (D_s)(M_s)
      DOUBLE PRECISION trDM_s
!
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP,&
                       IJKW, IJKE, IJKS, IJKN, IJKB, IJKT,&
                       IM, JM, KM
      INTEGER          IMJPK, IMJMK, IMJKP, IMJKM, IPJKM, IPJMK, IJMKP,&
                       IJMKM, IJPKM
!
!                      Solids phase
      INTEGER          M
!
! start anuj 04/20
!                      Used to compute frictional terms
      DOUBLE PRECISION Chi, Pc, Mu_f, Lambda_f, Pf, Mu_zeta,Phin,PfoPc
      DOUBLE PRECISION ZETA
! end anuj 04/20
!                      Use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION Mu, Mu_b, Mu_star, MUs, Kth, Kth_star
! start sof@fluent.com
!
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
      DOUBLE PRECISION D_12 
! 
      DOUBLE PRECISION Tau_12 
      DOUBLE PRECISION Tau_1 
      DOUBLE PRECISION C_beta, zeta_c
      DOUBLE PRECISION Zeta_r, Etha_12, SIGMA_12, B, CV
      DOUBLE PRECISION Diss 
      DOUBLE PRECISION K_1, K_2 
      DOUBLE PRECISION Tau_12_st
      DOUBLE PRECISION C_mu 
      DOUBLE PRECISION MU_12_T, MU_2_T, SIGMA_c, Tau_2_c 
       
 
!     SWITCH enables us to turn on/off the modification to the
!     particulate phase viscosity. If we want to simulate gas-particle
!     flow then SWITCH=1 to incorporate the effect of drag on the
!     particle viscosity. If we want to simulate granular flow
!     without the effects of an interstitial gas, SWITCH=0.
!     (Same for conductivity)
 
!  Function subroutines
!
!
!                      dg0/dep
      DOUBLE PRECISION DG_0DNU,SRT
!                  
      
!
!     Include statement functions
!
!
      INCLUDE 's_pr1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 's_pr2.inc'
     
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
!
!

!!$omp  parallel do &
!!$omp& private(IMJPK, I, J, K, IJK,  IMJK, IPJK, IJMK, IJPK, IJKM, &
!!$omp&  IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT, IM, JM, KM, &
!!$omp&  U_s_N, U_s_S, U_s_T, U_s_B, V_s_E, V_s_W, V_s_T, V_s_B, W_s_N, &
!!$omp&  W_s_S, W_s_E, W_s_W, U_s_C, W_s_C, D_s, I2_devD_s, trD_s_C, &
!!$omp&  qxP_s, trD_s2, K_1m, K_2m, K_3m, K_4m, K_5m, aq, bq, cq, &
!!$omp&  DEP_soDX, DEP_soDY, DEP_soXDZ, M_s, trM_s, trDM_s, I1, I2, &
!!$omp&  KTH_STAR,KTH,CHI,PFOPC,PC,ZETA,MU_ZETA,PF,LAMBDA_F,MU_F,MUS, &
!!$omp&  MU_STAR,MU_B,MU,M,IJPKM,IJMKM,IJMKP,IPJMK,IPJKM,IMJKM,IMJKP,IMJMK, &
!!$omp&  EP_sxSQRTHETA, EP_s2xTHETA )  
      DO 200 IJK = ijkstart3, ijkend3       
!
        IF ( FLUID_AT(IJK) ) THEN
 	
!
!------------------------------------------------------------------------
!          CALL SET_INDEX1(IJK, I, J, K, IMJK, IPJK, IJMK, IJPK,
!     &                       IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,
!     &                       IJKB, IJKT, IM, JM, KM)
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)
          IM = Im1(I)
          JM = Jm1(J)
          KM = Km1(K)
          IJKW  = WEST_OF(IJK)
          IJKE  = EAST_OF(IJK)
          IJKS  = SOUTH_OF(IJK)
          IJKN  = NORTH_OF(IJK)
          IJKB  = BOTTOM_OF(IJK)
          IJKT  = TOP_OF(IJK)
          IMJK  = IM_OF(IJK)
          IPJK  = IP_OF(IJK)
          IJMK  = JM_OF(IJK)
          IJPK  = JP_OF(IJK)
          IJKM  = KM_OF(IJK)
          IJKP  = KP_OF(IJK)
          IMJPK = IM_OF(IJPK)
          IMJMK = IM_OF(IJMK)
          IMJKP = IM_OF(IJKP)
          IMJKM = IM_OF(IJKM)
          IPJKM = IP_OF(IJKM)
          IPJMK = IP_OF(IJMK)
          IJMKP = JM_OF(IJKP)
          IJMKM = JM_OF(IJKM)
          IJPKM = JP_OF(IJKM)

	

          U_s_N = AVG_Y(                                   &   !i, j+1/2, k
                   AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I),&
                   AVG_X_E(U_s(IMJPK, M), U_s(IJPK, M), I), J&
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
	END IF
! end loezos
          V_s_T = AVG_Z(                                 &     !i, j, k+1/2
                   AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),&
                   AVG_Y_N(V_s(IJMKP, M), V_s(IJKP, M)), K&
                 )
          V_s_B = AVG_Z(                                 &    !i, j, k-1/2
                   AVG_Y_N(V_s(IJMKM, M), V_s(IJKM, M)),&
                   AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), KM&
                 )
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
!         Find components of Mth solids phase continuum strain rate
!         tensor, D_s, at center of THETA cell-(i, j, k)
          D_s(1,1) = ( U_s(IJK,M) - U_s(IMJK,M) ) * oDX(I)
          D_s(1,2) = HALF * ( (U_s_N - U_s_S) * oDY(J) +&
                               (V_s_E - V_s_W) * oDX(I) )
          D_s(1,3) = HALF * ( (W_s_E - W_s_W) * oDX(I) +&
                               (U_s_T - U_s_B) * (oX(I)*oDZ(K)) -&
                                W_s_C * oX(I) )
          D_s(2,1) = D_s(1,2)
          D_s(2,2) = ( V_s(IJK,M) - V_s(IJMK,M) ) * oDY(J)
          D_s(2,3) = HALF * ( (V_s_T - V_s_B) * (oX(I)*oDZ(K)) +&
                               (W_s_N - W_s_S) * oDY(J) )
          D_s(3,1) = D_s(1,3)
          D_s(3,2) = D_s(2,3)
          D_s(3,3) = ( W_s(IJK,M) - W_s(IJKM,M) ) * (oX(I)*oDZ(K)) +&
                      U_s_C * oX(I)


!
!         Calculate the trace of D_s
          trD_s_C(IJK,M) = D_s(1,1) + D_s(2,2) + D_s(3,3)

!
!         Calculate trace of the square of D_s
          trD_s2(IJK,M) = 0.D0  !Initialize the totalizer
          DO 20 I1 = 1,3
            DO 10 I2 = 1,3
              trD_s2(IJK,M) = trD_s2(IJK,M) + D_s(I1,I2)*D_s(I1,I2)
   10       CONTINUE
   20     CONTINUE
 
!Start Sof
!
	C_MU = 9.0D-02
	C_Beta=  0.45D0
	CV = 0.5D0
				
	 IF(Ep_s(IJK,M) .GT. SMALL_NUMBER	&
	 .AND. Scalar(IJK, 1) .GT. SMALL_NUMBER  		&
	 .AND. Scalar(IJK, 2) .GT. SMALL_NUMBER			&
	 .AND. Scalar(IJK, 3) .GT. SMALL_NUMBER) THEN
!
! Start definition of Relative Velocity
!
!         Calculate velocity components at i, j, k
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I) 
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK)) 
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK)) 
            USCM = AVG_X_E(U_S(IMJK,1),U_S(IJK,1),I) 
            VSCM = AVG_Y_N(V_S(IJMK,1),V_S(IJK,1)) 
            WSCM = AVG_Z_T(W_S(IJKM,1),W_S(IJK,1)) 
!
!         magnitude of gas-solids relative velocity
!
            VREL = SQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)
!
!
! Scalar1 is the gas phase turbulent kinetic energy (K).
	    K_1  = Scalar(IJK, 1)
!
! Scalar2 is the gas phase turbulence dissipation (Epsilon).	    
	    Diss = Scalar(IJK, 2)
!
! Scalar3 is the solids phase turbulent kinetic energy
	    K_2 = Scalar(IJK, 3)
	    
	    K_12(IJK) = Scalar(IJK,4)
!	
! Particle Relaxation time: Tau_12_st
! Just to make sure that the drag coef. isn't zero.

	IF (F_GS(IJK,1) .LT. SMALL_NUMBER) THEN
	  F_GS(IJK,1) = SMALL_NUMBER
	END IF
	
            Tau_12_st = Ep_s(IJK,M)*RO_g(IJK)/F_GS(IJK,1)* 	&
	    		(RO_s(M)/RO_g(IJK)+CV)
!
! This is Zeta_r**2 as defined by Simonin and fluent6 man.
! No need to take the square root of it and raise it to the power two again.
            Zeta_r    = 3.0d0 * VREL**2 / (2.0d0*K_1)
!
!
! Time scale of turbulent eddies: Tau_1 	    
            Tau_1     = 3.d0/2.d0*C_MU*K_1/Diss
!	    
! Lagrangian Integral time scale: Tau_12	    
            Tau_12    = Tau_1/sqrt(ONE+C_Beta*Zeta_r)
!
! Defining the inter-particle collision time
	    Tau_2_c = D_P(M)/6.d0/Ep_s(IJK,M)/SQRT(32.d0/3.d0/PI*K_2)
	    
! Coefficient B used below, CV is added mass coef.
	    B = (ONE + CV)/(RO_s(M)/RO_g(IJK)+CV)  
!
! 
	  Sigma_c = (ONE+ C_e)*(3.d0-C_e)/5.d0
!	
! Zeta_c: const. to be used in the K_2 Diffusion coefficient.
	  Zeta_c  = (ONE+ C_e)*(49.d0-33.d0*C_e)/100.d0
!
! Defining K_12(IJK)(IJK) Turbulent Viscosity: MU_12_T
	  MU_12_T = 1.d0/3.d0*RO_s(M)*K_12(IJK)*Tau_12
!
! Defining the Solids Turbulent Viscosity: MU_2_T
	  MU_2_T = (MU_12_T+RO_s(M)/3.d0*Tau_12_st*K_2) / 	 &
		 (ONE + Tau_12_st/2.d0*SIGMA_c/Tau_2_c)
!
	ELSE
	  MU_2_T = ZERO 
	END IF 
 
! start anuj 04/20
 
!  GRANULAR_ENERGY
!  .FALSE.
!         EP_g < EP_star   -->    plastic
!         EP_g >= EP_star  -->    viscous (algebraic)
 
!  GRANULAR_ENERGY
!  .TRUE.
!
!        FRICTION
!        .TRUE.
!              EP_s(IJK,M) > EPS_f_min  -->  friction + viscous(pde)
!              EP_s(IJK,M) < EP_f_min   -->  viscous (pde)
 
!        FRICTION
!       .FALSE.
!              EP_g < EP_star  -->  plastic + viscous(pde)
!              EP_g >= EP_star -->  viscous (pde)
 
! end anuj 04/20
 
          Mu_s(IJK,M)     = ZERO
          LAMBDA_s(IJK,M) = ZERO
!
 
! start anuj 4/20
         IF (.NOT.FRICTION) THEN
! end anuj 4/20
 
 
 
          IF(EP_g(IJK) .LT. EP_star) THEN
!             P_star(IJK) = Neg_H(EP_g(IJK))
!
!  Plastic-flow stress tensor
!
!           Calculate the second invariant of the deviator of D_s
            I2_devD_s = ( (D_s(1,1)-D_s(2,2))**2&
                          +(D_s(2,2)-D_s(3,3))**2&
                          +(D_s(3,3)-D_s(1,1))**2 )/6.&
                          + D_s(1,2)**2 + D_s(2,3)**2 + D_s(3,1)**2
!
!-----------------------------------------------------------------------
!            Gray and Stiles (1988)
!            IF(Sin2_Phi .GT. SMALL_NUMBER) THEN
!              qxP_s = SQRT( (4. * Sin2_Phi) * I2_devD_s
!     &                       + trD_s_C(IJK,M) * trD_s_C(IJK,M))
!              MU_s(IJK, M)     = P_star(IJK) * Sin2_Phi
!     &                         / (qxP_s + SMALL_NUMBER)
!              MU_s(IJK, M)     = MIN(MU_s(IJK, M), MAX_MU_s)
!              LAMBDA_s(IJK, M) = P_star(IJK) * F_Phi
!     &                         / (qxP_s + SMALL_NUMBER)
!              LAMBDA_s(IJK, M) = MIN(LAMBDA_s(IJK, M), MAX_MU_s)
!            ELSE
!              MU_s(IJK, M)     = ZERO
!              LAMBDA_s(IJK, M) = ZERO
!            ENDIF
!-----------------------------------------------------------------------
!           Schaeffer (1987)
!
            qxP_s            = SQRT( (4. * Sin2_Phi) * I2_devD_s)
            MU_s(IJK, M)     = P_star(IJK) * Sin2_Phi&
                               / (qxP_s + SMALL_NUMBER)
            MU_s(IJK, M)     = MIN(MU_s(IJK, M), MAX_MU_s)
 
            LAMBDA_s(IJK, M) = ZERO
            ALPHA_s(IJK, M)  = ZERO
            P_s(IJK, M)  = ZERO
	    THETA_m(IJK, M) = ZERO
          ENDIF
         ENDIF
 
 
 
 
 
!
!  Viscous-flow stress tensor
!
          IF(.NOT.GRANULAR_ENERGY) THEN  !algebraic granular energy equation
            IF(EP_g(IJK) .GE. EP_star) THEN
!
!
!             Calculate K_1m, K_2m, K_3m, K_4m
              K_1m = 2.D0 * (ONE + C_e) * RO_s(M) * G_0(IJK, M,M)
              K_3m = HALF * D_p(M) * RO_s(M) * (&
                  ( (SQRT_PI / (3.D0*(3.D0 - C_e))) *&
                  (ONE + 0.4D0*(ONE + C_e)*(3.D0*C_e - ONE)*&
                  EP_s(IJK,M)*G_0(IJK, M,M)) ) +&
                  8.D0*EP_s(IJK,M)*G_0(IJK, M,M)*(ONE + C_e) /&
                  (5.D0*SQRT_PI) )
              K_2m = 4.D0 * D_p(M) * RO_s(M) * (ONE + C_e) *&
                  EP_s(IJK,M) * G_0(IJK, M,M) / (3.D0 * SQRT_PI) -&
                  2.D0/3.D0 * K_3m
              K_4m = 12.D0 * (ONE - C_e*C_e) *&
                  RO_s(M) * G_0(IJK, M,M) / (D_p(M) * SQRT_PI)
              aq   = K_4m*EP_s(IJK,M)
              bq   = K_1m*EP_s(IJK,M)*trD_s_C(IJK,M)
              cq   = -(K_2m*trD_s_C(IJK,M)*trD_s_C(IJK,M)&
                     + 2.D0*K_3m*trD_s2(IJK,M))
!
!             Boyle-Massoudi Stress term
!
              IF(V_ex .NE. ZERO) THEN
                K_5m = 0.4 * (ONE + C_e) * G_0(IJK, M,M) * RO_s(M) *&
                  ( (V_ex * D_p(M)) / (ONE - EP_s(IJK,M) * V_ex) )**2
                DEP_soDX  = ( EP_s(IJKE, M) - EP_s(IJK, M) ) * oDX_E(I)&
                 * ( ONE / ( (oDX_E(IM)/oDX_E(I)) + ONE ) ) +&
                 ( EP_s(IJK, M) - EP_s(IJKW, M) ) * oDX_E(IM)&
                 * ( ONE / ( (oDX_E(I)/oDX_E(IM)) + ONE ) )
                DEP_soDY  = ( EP_s(IJKN, M) - EP_s(IJK, M) ) * oDY_N(J)&
                 * ( ONE / ( (oDY_N(JM)/oDY_N(J)) + ONE ) ) +&
                 ( EP_s(IJK, M) - EP_s(IJKS, M) ) * oDY_N(JM)&
                 * ( ONE / ( (oDY_N(J)/oDY_N(JM)) + ONE ) )
                DEP_soXDZ  = (( EP_s(IJKT, M) - EP_s(IJK, M) )&
                 * oX(I)*oDZ_T(K)&
                 * ( ONE / ( (oDZ_T(KM)/oDZ_T(K)) + ONE ) ) +&
                 ( EP_s(IJK, M) - EP_s(IJKB, M) ) * oX(I)*oDZ_T(KM)&
                 * ( ONE / ( (oDZ_T(K)/oDZ_T(KM)) + ONE ) ) ) /&
                 X(I)
                M_s(1,1) = DEP_soDX * DEP_soDX
                M_s(1,2) = DEP_soDX * DEP_soDY
                M_s(1,3) = DEP_soDX * DEP_soXDZ
                M_s(2,1) = DEP_soDX * DEP_soDY
                M_s(2,2) = DEP_soDY * DEP_soDY
                M_s(2,3) = DEP_soDY * DEP_soXDZ
                M_s(3,1) = DEP_soDX * DEP_soXDZ
                M_s(3,2) = DEP_soDY * DEP_soXDZ
                M_s(3,3) = DEP_soXDZ * DEP_soXDZ
                trM_s    = M_s(1,1) + M_s(2,2) + M_s(3,3)
                trDM_s = ZERO
                DO 40 I1 = 1,3
                  DO 30 I2 = 1,3
                    trDM_s = trDM_s + D_s(I1,I2)*M_s(I1,I2)
   30             CONTINUE
   40           CONTINUE
                bq   = bq + EP_s(IJK,M) * K_5m * (trM_s + 2. * trDM_s)
              ELSE
                K_5m = ZERO
              ENDIF
!
!             Calculate EP_sxSQRTHETA and EP_s2xTHETA
              EP_sxSQRTHETA = (-bq + SQRT(bq**2 - 4. * aq * cq ))&
                             / ( 2. * K_4m )
              EP_s2xTHETA = EP_sxSQRTHETA * EP_sxSQRTHETA
 
              IF(EP_s(IJK,M) > SMALL_NUMBER)THEN
!start      kapil&anuj 01/19/98
!               Find pseudo-thermal temperature in the Mth solids phase
                THETA_m(IJK,M) = EP_s2xTHETA/(EP_s(IJK,M)*EP_s(IJK,M))
!end      kapil&anuj 01/19/98
              ELSE
                THETA_m(IJK,M) = ZERO
              ENDIF
!
!             Find pressure in the Mth solids phase
              P_s(IJK,M) = K_1m * EP_s2xTHETA
!
!             bulk viscosity in Mth solids phase
              LAMBDA_s(IJK,M) = K_2m * EP_sxSQRTHETA +2d0/3d0*Mu_2_T*EP_s(IJK,M)
!
!             shear viscosity in Mth solids phase
              MU_s(IJK,M) = K_3m * EP_sxSQRTHETA + Mu_2_T*EP_s(IJK,M)
!
!             Boyle-Massoudi stress coefficient
              ALPHA_s(IJK, M) = -K_5m * EP_s2xTHETA
	    ENDIF
 
          ELSE	!granular energy transport equation
 
 
!           Find pressure in the Mth solids phase
            P_s(IJK,M) = ROP_s(IJK,M)*(1d0+ 4. * Eta *&
                           EP_s(IJK,M)*G_0(IJK,M,M))*Theta_m(IJK,M)
!
            Mu = (5d0*DSQRT(Pi*Theta_m(IJK,M))*D_p(M)*RO_s(M))/96d0
 
            Mu_b = (256d0*Mu*EP_s(IJK,M)*EP_s(IJK,M)*G_0(IJK,M,M))&
                     /(5d0*Pi)
 
            IF(SWITCH*F_gs(IJK,M) .LT. SMALL_NUMBER)THEN
              Mu_star = Mu
		
	    ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
              Mu_star = ZERO
	
	    ELSE
              Mu_star = Mu/(1.+(2d0*SWITCH*F_gs(IJK,M)*Mu/&
                        (RO_s(M)*RO_s(M)*EP_s(IJK,M)*EP_s(IJK,M)&
                        *G_0(IJK,M,M)*Theta_m(IJK,M))))
            ENDIF
!
!           shear viscosity in Mth solids phase  (add to plastic part)
            Mus =&
                   ((2d0+ALPHA)/3d0)*((Mu_star/(Eta*(2d0-Eta)*&
                   G_0(IJK,M,M)))*(1d0+1.6d0*Eta*EP_s(IJK,M)*&
                   G_0(IJK,M,M))*(1d0+1.6d0*Eta*(3d0*Eta-2d0)*&
                   EP_s(IJK,M)*G_0(IJK,M,M))+(0.6d0*Mu_b*Eta))
 
 
! start anuj 04/20
!           calculate frictional stress
 
            Mu_f = ZERO
            Lambda_f = ZERO
            Pf = ZERO
 
            IF (FRICTION) THEN
               IF (EP_s(IJK,M) .GT. EPS_f_min) THEN
 
                  IF (SAVAGE.EQ.1) THEN !form of Savage
            	     Mu_zeta =&
                           ((2d0+ALPHA)/3d0)*((Mu/(Eta*(2d0-Eta)*&
                           G_0(IJK,M,M)))*(1d0+1.6d0*Eta*EP_s(IJK,M)*&
                           G_0(IJK,M,M))*(1d0+1.6d0*Eta*(3d0*Eta-2d0)*&
                           EP_s(IJK,M)*G_0(IJK,M,M))+(0.6d0*Mu_b*Eta))
 
 
                     ZETA =&
                            ((48d0*Eta*(1d0-Eta)*RO_s(M)*EP_s(IJK,M)*&
                            EP_s(IJK,M)*G_0(IJK,M,M)*&
                            (Theta_m(IJK,M)**1.5d0))/&
                            (SQRT_Pi*D_p(M)*2d0*Mu_zeta))**0.5d0
 
                  ELSEIF (SAVAGE.EQ.0) THEN  !S:S form
                     ZETA = (SMALL_NUMBER +&
                             trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                             trD_s_C(IJK,M))/3.d0))**0.5d0
 
                  ELSE  !combined form
 
                     ZETA = ((Theta_m(IJK,M)/(D_p(M)*D_p(M))) +&
                            (trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                             trD_s_C(IJK,M))/3.d0)))**0.5d0
 
                  ENDIF
 
 
                  IF (EP_s(IJK,M) .GT. EPS_max) THEN
                     Pc = 1d25*((EP_s(IJK,M)- EPS_max)**10d0)
                  ELSE
                     Pc = Fr*((EP_s(IJK,M) - EPS_f_min)**N_Pc)/&
                          ((EPS_max - EP_s(IJK,M) + SMALL_NUMBER)&
                           **D_Pc)
                  ENDIF
 
 
                  IF ((trD_s_Co(IJK,M)/(ZETA*N_Pf*DSQRT(2d0)&
                      *Sin_Phi))&
                       .GT. 1d0) THEN
                     Pf =ZERO
                     PfoPc = ZERO
                  ELSE
 
                     Pf = Pc*(1d0 - (trD_s_Co(IJK,M)/(ZETA&
                          *N_Pf*&
                          DSQRT(2d0)*Sin_Phi)))**(N_Pf-1d0)
 
                     PfoPc = (1d0 - (trD_s_Co(IJK,M)/(ZETA&
                          *N_Pf*&
                          DSQRT(2d0)*Sin_Phi)))**(N_Pf-1d0)
                  ENDIF
 
 
 
                  Chi = DSQRT(2d0)*Pf*Sin_Phi*(N_Pf - (N_Pf-1d0)*&
                                      (PfoPc)**(1d0/(N_Pf-1d0)))
 
                  IF (Chi < ZERO) THEN
                     Pf = Pc*((N_Pf/(N_Pf-1d0))**(N_Pf-1d0))
                     Chi = ZERO
                  ENDIF
 
 
                  Mu_f = Chi/(2d0*ZETA)
                  Lambda_f = - 2d0*Mu_f/3d0
 
 
               ENDIF
            ENDIF
 
            Mu_s_c(IJK,M) = MUs
            Mu_s(IJK,M) = Mu_s(IJK,M) + MUs + Mu_f + MU_2_T
 
!
!           bulk viscosity in Mth solids phase   (add to plastic part)
 
            LAMBDA_s_c(IJK,M)= Eta*Mu_b - (2d0*MUs/3d0)
            LAMBDA_s(IJK,M) = LAMBDA_s(IJK,M) &
                                + Eta*Mu_b - (2d0*MUs/3d0) + Lambda_f &
				+2.D0/3.D0*MU_2_T
 
            P_s_c(IJK,M) = P_s(IJK,M)
            P_s(IJK,M) = P_s(IJK,M) + Pf      !add to P_s
! end anuj 04/20
 
            Kth=75*RO_s(M)*D_p(M)*DSQRT(Pi*Theta_m(IJK,M))/&
                  (48*Eta*(41d0-33*Eta))
 
            IF(SWITCH*F_gs(IJK,M) .LT. SMALL_NUMBER)THEN
              Kth_star=Kth
		
	    ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
              Kth_star = ZERO
	
	    ELSE
              Kth_star=Kth/(1.+(1.2d0*SWITCH*F_gs(IJK,M)*Kth/&
                       (RO_s(M)*RO_s(M)*EP_s(IJK,M)*EP_s(IJK,M)&
                       *G_0(IJK,M,M)*Theta_m(IJK,M))))
            ENDIF
!
!           granular conductivity in Mth solids phase
            Kth_s(IJK,M) = Kth_star*(&
                  ( (1d0/G_0(IJK,M,M)) + (12d0/5.)*Eta*EP_s(IJK,M) )&
                  * ( (1d0) + (12d0/5.)*Eta*Eta*(4d0*Eta-3d0)&
                      *EP_s(IJK,M)*G_0(IJK,M,M) )&
                  + (64d0/(25d0*Pi)) * (41d0-33d0*Eta) *&
                     (Eta*EP_s(IJK,M))**2 * G_0(IJK,M,M)&
              )
 
!
!     granular 'conductivity' in the Mth solids phase associated
!     with gradient in volume fraction
 
!--------------------------------------------------------------------
!  Kphi_s has been set to zero.  To activate the feature uncomment the
!  following lines and also the lines in source_granular_energy.
            Kphi_s(IJK,M) = ZERO
!     &            (Kth_star/(G_0(IJK,M,M)))*(12d0/5.)*Eta*(Eta-1.)*
!     &            (2.*Eta-1.)*(1.+(12d0/5.)*Eta*EP_s(IJK,M)*
!     &            G_0(IJK,M,M))*(EP_s(IJK,M)*
!     &            DG_0DNU(EP_s(IJK,M))
!     &            + 2*G_0(IJK,M,M))*Theta_m(IJK,M)
!--------------------------------------------------------------------
 
!
!           Boyle-Massoudi stress coefficient
            ALPHA_s(IJK, M) = ZERO
          ENDIF
!
        ENDIF
!
  200 CONTINUE
!
! loezos 
   IF (SHEAR) THEN
!$omp parallel do private(IJK)
      Do IJK= ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  	 
	   V_s(IJK,m)=V_s(IJK,m)-VSH(IJK)	
	 END IF
      END DO
   END IF
! loezos
      RETURN
      END

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
