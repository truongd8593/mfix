!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Scalar_PROP(IER)                                       C
!  Purpose: Calculate diffusion coefficeint and sources for user-defined
!           scalars
!                                                                      C
!  Author:                                                    Date:    C
!  Reviewer:
!  Modified: Sof@fluent.com                         Date:05/29/2002    C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SCALAR_PROP( IER) 
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
      USE scalars
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
	DIMENSION USX3(IMAX), USX4(IMAX)
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
                      IJMKM, IJPKM, II
!		    
!                      Solids phase
      INTEGER          M  
!
!
!                      Strain rate tensor components for mth solids phase
      DOUBLE PRECISION D_g(3,3)
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
      INTEGER          I1, I2
!
!                      Second invariant of the deviator of D_s
      DOUBLE PRECISION I2_devD_s, SRT, Trace_G, Trace_S
! 
!
! 
!
      DOUBLE PRECISION D_12 
! 
      DOUBLE PRECISION Tau_12 
      DOUBLE PRECISION Tau_1 
      DOUBLE PRECISION C_beta 
      DOUBLE PRECISION Zeta_r, Etha_12, SIGMA_12, B, CV,V_12V_dr,C_Eps_3
      DOUBLE PRECISION Diss 
      DOUBLE PRECISION K_1, K_2 
      DOUBLE PRECISION Tau_12_st
      DOUBLE PRECISION C_eps 
      DOUBLE PRECISION VREL_new
!			Production of Turb. Due to shear, Turb Visc, Const.
!  See Wilcox PP. 89
      DOUBLE PRECISION Tauij_gDUi_gODxj, C_MU, Sigma_k, Sigma_E
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
      DOUBLE PRECISION U_11_U_11,U_12_U_11,U_13_U_11
      DOUBLE PRECISION U_12_U_12,U_12_U_13,U_13_U_13
      
      DOUBLE PRECISION U_11_U_21,U_12_U_21,U_13_U_21
      DOUBLE PRECISION U_11_U_22,U_12_U_22,U_13_U_22
      DOUBLE PRECISION U_11_U_23,U_12_U_23,U_13_U_23
      
      DOUBLE PRECISION Tau_gi_sj, Pos_Tau_gi_sj, Neg_Tau_gi_sj
      DOUBLE PRECISION PI_q_12, Pos_PI_q_12, Neg_PI_q_12
      
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
!

      IF(NScalar == 0) RETURN
!
!
!!! writing to these files to check the values of variables
!
	!OPEN (UNIT=71,FILE='Scal1.dat',STATUS='UNKNOWN')
	!OPEN (UNIT=72,FILE='Scal2.dat',STATUS='UNKNOWN')
!
! When I don't state that M = 1, I found that the code is giving M = 2!!
! Since M = 2 doesn't mean anything, all I got was garbage.
! So for now, M should be forced to be equal to one.
	  M = ONE
!	  
	  K_1 = SCALAR(IJK,1)
	  Diss = SCALAR(IJK,2)
	  K_2 = SCALAR(IJK,3)
! Overwrite values computed for K_12	  
	  K_12(IJK) = SCALAR(IJK,4)
!	  
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
!
!        Calculate Tau(i,j)*dUi/dXj (production term in the K Equation
! Sof@fluent.com
!
	C_MU = 9.0D-02
	Sigma_k = 1.0D0
	Sigma_E = 1.3D0
	Ceps_1 = 1.44D0
	Ceps_2 = 1.92D0
	SIGMA_12 = 0.67D0
	CV = 0.5D0
	C_Eps_3 = 1.2D0
! C_Beta as defined by Eq. 20-4-77 from Fluent6 manual with Cos(Theta)=1
	C_Beta=  0.45D0	
	X_21 = Ep_s(IJK,M)*RO_s(M)/EP_g(IJK)/RO_g(IJK)
!!
	Trace_g = D_G(1,1) + D_G(2,2) + D_G(3,3)
!		   
	Tauij_gDUi_gODxj = 2D0*MU_GT(IJK)*(				&
		       D_G(1,1) * D_G(1,1) +				&
		       D_G(1,2) * (U_G_N - U_G_S)*ODY(J) +		&
		       D_G(1,3) * ((U_G_T-U_G_B)*			&
		       		  (OX(I)*ODZ(K))-W_G_C*OX(I)) +		&
		       D_G(2,1) * (V_G_E-V_G_W)*ODX(I) +		&
		       D_G(2,2) * D_G(2,2) +				&
		       D_G(2,3) * (V_G_T-V_G_B)*(OX(I)*ODZ(K)) +	&
		       D_G(3,1) * (W_G_E - W_G_W)*ODX(I) + 		&
		       D_G(3,2) * (W_G_N-W_G_S)*ODY(J) + 		&
		       D_G(3,3) * D_G(3,3)) -				&
		       2.d0/3.d0 * RO_G(IJK) * Scalar(IJK, 1)*Trace_g	&
		       -2.d0/3.d0 * MU_GT(IJK) * Trace_g**2
	
	U_11_U_11 = -2D0*MU_GT(IJK)/RO_g(IJK)*D_G(1,1) + 2.d0/3.d0 *	&
		    (Scalar(IJK, 1)+MU_GT(IJK)/RO_G(IJK) * Trace_g)
		    
	U_12_U_11 = -2D0*MU_GT(IJK)/RO_g(IJK)* 				&
	           ((U_G_N - U_G_S)*ODY(J)+(V_G_E-V_G_W)*ODX(I))
		   
	U_13_U_11 = -2D0*MU_GT(IJK)/RO_g(IJK)*				&
		   ((W_G_E - W_G_W)*ODX(I)+(U_G_T-U_G_B)*(OX(I)*ODZ(K)	&
                    )-W_G_C*OX(I))
	
	U_12_U_12 = -2D0*MU_GT(IJK)/RO_g(IJK)*D_G(2,2) + 2.d0/3.d0 *	&
		    (Scalar(IJK, 1)+MU_GT(IJK)/RO_G(IJK) * Trace_g)	     
	
	U_12_U_13 = -2D0*MU_GT(IJK)/RO_g(IJK)* 				&
	           ((V_G_T-V_G_B)*(OX(I)*ODZ(K))+(W_G_N-W_G_S)*ODY(J))
		   
	U_13_U_13 = -2D0*MU_GT(IJK)/RO_g(IJK)*D_G(3,3) + 2.d0/3.d0 *	&
		    (Scalar(IJK, 1)+MU_GT(IJK)/RO_G(IJK) * Trace_g)
		      
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
!!!
!------------------
! Part copied from calc_mu_s.f to calculate the rate of strain Sii tensor
! of the Mth solids phase.
! Sof@fluent.com
!------------------
!
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
          D_s(1,1) = ( U_s(IJK,M) - U_s(IMJK,M) ) * oDX(I)         !du1odx1
          D_s(1,2) = HALF * ( (U_s_N - U_s_S) * oDY(J) +&          !du1dx2
                               (V_s_E - V_s_W) * oDX(I) )	   !du2odx1
          D_s(1,3) = HALF * ( (W_s_E - W_s_W) * oDX(I) +&          !du3odx1
                               (U_s_T - U_s_B) * (oX(I)*oDZ(K)) -& !du1odx3
                                W_s_C * oX(I) )
          D_s(2,1) = D_s(1,2)
          D_s(2,2) = ( V_s(IJK,M) - V_s(IJMK,M) ) * oDY(J)	   !du2odx2
          D_s(2,3) = HALF * ( (V_s_T - V_s_B) * (oX(I)*oDZ(K)) +&  !du2odx3
                               (W_s_N - W_s_S) * oDY(J) )	   !du3odx2
          D_s(3,1) = D_s(1,3)
          D_s(3,2) = D_s(2,3)
          D_s(3,3) = ( W_s(IJK,M) - W_s(IJKM,M) ) * (oX(I)*oDZ(K)) +& !du3odx3
                      U_s_C * oX(I)


!
!         Calculate the trace of D_s
          trace_S = D_s(1,1) + D_s(2,2) + D_s(3,3)
!
! Before calculating the Reynold's stresses in the Solids phase following
! Simonin's approach, I'll first allow the calculation of some characteristic
! Time scales since I'll need to use them.
!	
!!!
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
! The drift velocity is defined only if the solid phase and turbulence
! quantities exist.
!
	 IF(Ep_s(IJK,M) .GT. SMALL_NUMBER	 		&
	 .AND. Scalar(IJK, 1) .GT. SMALL_NUMBER  		&
	 .AND. Scalar(IJK, 2) .GT. SMALL_NUMBER			&
	 .AND. Scalar(IJK, 3) .GT. SMALL_NUMBER) THEN
!
! Scalar1 is the gas phase turbulent kinetic energy (K).
	    K_1  = Scalar(IJK, 1)
!
! Scalar2 is the gas phase turbulence dissipation (Epsilon).	    
	    Diss = Scalar(IJK, 2)
!
! Scalar3 is the solids phase turbulent kinetic energy
	    K_2 = Scalar(IJK, 3)
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
!Definition of ratio Tau_12/Tau_12_st: Etha_12
	    Etha_12 =  Tau_12/Tau_12_st 
!
! Defining the inter-particle collision time
	    Tau_2_c = D_P(M)/6.d0/Ep_s(IJK,M)/SQRT(32.d0/3.d0/PI*K_2)
	    
! Coefficient B used below, CV is added mass coef.
	    B = (ONE + CV)/(RO_s(M)/RO_g(IJK)+CV)  
!
! K_12(IJK) Should be Scalar4, cross correlation between gas and solid turb.
! For now, I'm using definition from Fluent6			
!            K_12(IJK)      = 2.D0*K_1*(B+Etha_12)/(1.d0+Etha_12)
!	    
! K_2 will be theta_m, gran. temp. Not used for now
! Below is a definition of K_2 from Fluent6 man.	    
!	    K_2       = K_1*(B**2+Etha_12)/(1.d0+Etha_12) !not used	    
!
!
! Def of Dissipation in the Solids turbulence equation.
	Diss_2 = 1.d0/3.d0*(ONE- C_e**2)/Tau_2_c*K_2
!				
! C_e Rest. Coef.
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
	MU_2_T = (MU_12_T+RO_s(M)/3.d0*Tau_12_st*K_2) / 		&
		 (ONE + Tau_12_st/2.d0*SIGMA_c/Tau_2_c)
!



	ELSE
	  MU_2_T = ZERO

	END IF  ! For Eps_S, K and Epsilon being very small
!	 
            VREL_new  = SQRT((UGC + U_dx(IJK) - USCM)**2 + &
                        (VGC + U_dy(IJK) - VSCM)**2 + (WGC + U_dz(IJK) - WSCM)**2)
!			 
!The last term in Eq. 20-4-72 Fluent6 manual
	    V_12V_dr= (UGC - USCM)*U_dx(IJK)+(VGC - VSCM)*U_dy(IJK)+	&
	    		(WGC - WSCM)*U_dz(IJK)
!
!		   
	Tauij_sDUi_sODxj = 2D0*MU_2_T*(					&
		       D_s(1,1) * D_s(1,1) +				&
		       D_s(1,2) * (U_s_N - U_s_S)*ODY(J) +		&
		       D_s(1,3) * ((U_s_T-U_s_B)*			&
		       		  (OX(I)*ODZ(K))-W_s_C*OX(I)) +		&
		       D_s(2,1) * (V_s_E-V_s_W)*ODX(I) +		&
		       D_s(2,2) * D_s(2,2) +				&
		       D_s(2,3) * (V_s_T-V_s_B)*(OX(I)*ODZ(K)) +	&
		       D_s(3,1) * (W_s_E - W_s_W)*ODX(I) + 		&
		       D_s(3,2) * (W_s_N-W_s_S)*ODY(J)  + 		&
		       D_s(3,3) * D_s(3,3)) 				&
		       -2.d0/3.d0 * RO_s(M) * Scalar(IJK, 3)* trace_s   &
		       -2.d0/3.d0 * MU_2_T * Trace_s**2
		      
! To avoid very small negative numbers.
	IF(Tauij_sDUi_sODxj .GE. ZERO) THEN
	  Pos_Tauij_sDUi_sODxj = Tauij_sDUi_sODxj
	  Neg_Tauij_sDUi_sODxj = ZERO
	ELSE
	  Pos_Tauij_sDUi_sODxj = ZERO
	  Neg_Tauij_sDUi_sODxj = Tauij_sDUi_sODxj
!
	END IF   
!
!! Implementing the Reynold's stresses in the K_12 Equation.. Good luck!
!
	U_11_U_21 = ONE/3.d0*K_12(IJK) + Etha_12/(ONE+Etha_12)* (	&
		   U_11_U_11-2.d0/3.d0*SCALAR(IJK,1)) - MU_2_T/RO_s(M)  &
		   /(ONE+Etha_12)*( D_G(1,1)+D_s(1,1)			&
		   -ONE/3.d0*trace_g- ONE/3.d0*trace_s )
		   
	U_12_U_21 = Etha_12/(ONE+Etha_12)* U_12_U_11 - MU_2_T/RO_s(M)   &
		    /(ONE+Etha_12)*( (V_g_E - V_g_W) * oDX(I) +	   	&
		    (U_s_N - U_s_S) * oDY(J) )
		    
	U_13_U_21 = Etha_12/(ONE+Etha_12)* U_13_U_11 - MU_2_T/RO_s(M)   &
		    /(ONE+Etha_12)*((W_g_E - W_g_W) * oDX(I) +	   	&
		    (U_s_T - U_s_B) * (oX(I)*oDZ(K)) -W_s_C * oX(I))
		    
	U_11_U_22 = Etha_12/(ONE+Etha_12)* U_12_U_11 - MU_2_T/RO_s(M)   &
		    /(ONE+Etha_12)*((U_G_N - U_G_S)*ODY(J)  +	   	&
		    (V_s_E - V_s_W) * oDX(I))
	
	U_12_U_22 = ONE/3.d0*K_12(IJK) + Etha_12/(ONE+Etha_12)* (	&
		   U_12_U_12-2.d0/3.d0*SCALAR(IJK,1)) - MU_2_T/RO_s(M)  &
		   /(ONE+Etha_12)* ( D_G(2,2)+D_s(2,2)			&
		   -ONE/3.d0*trace_g- ONE/3.d0*trace_s )
		   
	U_13_U_22 = Etha_12/(ONE+Etha_12)* U_12_U_13 - MU_2_T/RO_s(M)   &
		    /(ONE+Etha_12)*( (W_G_N-W_G_S)*ODY(J) +		&
		    (V_s_T - V_s_B) * (oX(I)*oDZ(K)) )
		    
	U_11_U_23 = Etha_12/(ONE+Etha_12)* U_13_U_11 - MU_2_T/RO_s(M)   &
		    /(ONE+Etha_12)*((U_G_T-U_G_B)*(OX(I)*ODZ(K))-	&
		    W_G_C*OX(I) + (W_s_E - W_s_W) * oDX(I) )
		    
	U_12_U_23 = Etha_12/(ONE+Etha_12)* U_12_U_13 - MU_2_T/RO_s(M)   &
		    /(ONE+Etha_12)*( (V_G_T-V_G_B)*(OX(I)*ODZ(K))+	&
		    (W_s_N - W_s_S) * oDY(J) )
		    
	U_13_U_23 = ONE/3.d0*K_12(IJK) + Etha_12/(ONE+Etha_12)* (	&
		   U_13_U_13-2.d0/3.d0*SCALAR(IJK,1)) - MU_2_T/RO_s(M)  &
		   /(ONE+Etha_12)* ( D_G(3,3)+D_s(3,3)			&
		   -ONE/3.d0*trace_g- ONE/3.d0*trace_s )
		   
	Tau_gi_sj = U_11_U_21 * ( ( U_s(IJK,M) - U_s(IMJK,M))*oDX(I)	&
		    + (U_G(IJK)-U_G(IMJK))*ODX(I) ) +			&
		    U_12_U_21 * ( (V_s_E - V_s_W) * oDX(I) +		&
		    (U_G_N - U_G_S)*ODY(J) ) +				&
		    U_13_U_21 * ( (W_s_E - W_s_W) * oDX(I) +		&
		    (U_G_T-U_G_B)*(OX(I)*ODZ(K))-W_G_C*OX(I) ) +	&
		    U_11_U_22 * ( (U_s_N - U_s_S) * oDY(J) +		&
		    (V_G_E-V_G_W)*ODX(I) ) +				&
		    U_12_U_22 * ( (V_s(IJK,M) - V_s(IJMK,M))*oDY(J) +   &
		    (V_G(IJK)-V_G(IJMK))*ODY(J) ) +			&
		    U_13_U_22 * ( (W_s_N - W_s_S) * oDY(J) +		&
		    (V_G_T-V_G_B)*(OX(I)*ODZ(K)) ) +			&
		    U_11_U_23 * ( (U_s_T - U_s_B)*(oX(I)*oDZ(K)) -	&
		    W_s_C * oX(I) + (W_G_E - W_G_W)*ODX(I) ) +		&
		    U_12_U_23 * ( (V_s_T - V_s_B) * (oX(I)*oDZ(K)) +	&
		    (W_G_N-W_G_S)*ODY(J) ) +				&
		    U_13_U_23 * ( (W_s(IJK,M)-W_s(IJKM,M))*		&
                    (oX(I)*oDZ(K)) + U_s_C * oX(I) +			&
		    (W_G(IJK)-W_G(IJKM))*(OX(I)*ODZ(K)) + U_G_C*OX(I) )
		    
	Tau_gi_sj = -Ep_s(IJK,M)*RO_s(M)*Tau_gi_sj
	
! To avoid very small negative numbers.
	IF(Tau_gi_sj .GE. ZERO) THEN
	  Pos_Tau_gi_sj = Tau_gi_sj
	  Neg_Tau_gi_sj = ZERO
	ELSE
	  Pos_Tau_gi_sj = ZERO
	  Neg_Tau_gi_sj = Tau_gi_sj
!
	END IF   
!		     
		   
	

! This is a write statment just to check the values of the data
!

 3      FORMAT(9(E11.5,1X))	
!		       
!
           DO L = 1, NScalar 
!            d (Scalar)/dt = S
!            S is linearized as S = Scalar_c - Scalar_p * Scalar
!            Scalar_c and Scalar_p must be >= 0 to insure that scalar>0
!!
! PI_kq	turbulence exchange betweeen gas and solids as defined by 
! Eq. 20-4-71 in fluent6 manual.
!
! PI_kq	has been separtated into terms 1 & 2 to ensure that term2 (which can
! be negative) is put in the right place.
	  PI_kq_1 = F_GS(IJK,1)*2.0D0*Scalar(IJK, 1)
	  PI_kq_2 = F_GS(IJK,1)*(K_12(IJK) + V_12V_dr)
!
	  IF(PI_kq_2 .GE. ZERO) THEN
	    Pos_PI_kq_2 = PI_kq_2
	    Neg_PI_kq_2 = ZERO
	  ELSE
	    Pos_PI_kq_2 = ZERO
	    Neg_PI_kq_2 = PI_kq_2 
	  END IF
	  
	  !PI_kq_1 = zero
	  !PI_kq_2 = 	  
!
! F_GS(IJK,1) : Drag Coef. computed by MFIX through the drag subroutine.
!	  			       	
!
!	  			       	
	IF (L.EQ.1) THEN
	IF (SCALAR(IJK,1).GT. SMALL_NUMBER .AND.  			&
	    SCALAR(IJK,2) .GT. SMALL_NUMBER) THEN
!
             Scalar_c (IJK, L) = (EP_g(IJK) * Pos_Tauij_gDUi_gODxj +	&
	     	                  Pos_PI_kq_2 )
!
             Scalar_p (IJK, L) =(-EP_g(IJK) * Neg_Tauij_gDUi_gODxj/	&
	     			Scalar(IJK, 1) +			&
				EP_g(IJK)*RO_G(IJK)*Scalar(IJK, 2) /	&
	     			Scalar(IJK, 1)	+	      		&
				(PI_kq_1-Neg_PI_kq_2)/Scalar(IJK, 1))	
!	
	ELSE
	     Scalar_c (IJK, L) = ZERO
	     Scalar_p (IJK, L) = ONE
	 END IF  !for scalars having very small values    
	  !WRITE(71,3) (SCALAR(IJK,1),SCALAR(IJK,2)) ! open files first
!
			  
!
!
!            Diffusion coefficient for User-defined Scalars
!            
             Dif_Scalar(IJK, L) =EP_g(IJK)*MU_GT(IJK) /Sigma_k
!
! Dissipation of Turbulence
!	     
	ELSE IF(L.EQ.2) THEN
	IF (SCALAR(IJK,1).GT. SMALL_NUMBER .AND.  			&
	    SCALAR(IJK,2) .GT. SMALL_NUMBER) THEN	
!
             Scalar_c (IJK, L) =  (Ceps_1 *EP_g(IJK) *Pos_Tauij_gDUi_gODxj &
	     	               *Scalar(IJK, 2)/Scalar(IJK, 1) +		&
			       C_Eps_3*Pos_PI_kq_2*			&
			       Scalar(IJK, 2)/Scalar(IJK, 1))
!
             Scalar_p (IJK, L) = (-Ceps_1 *EP_g(IJK) *Neg_Tauij_gDUi_gODxj &
	     	               /Scalar(IJK, 1) +			&
	     			Ceps_2 * EP_g(IJK) *RO_G(IJK)		&
			        *Scalar(IJK, 2)/Scalar(IJK, 1) 	        &
			         + C_Eps_3*(PI_kq_1-Neg_PI_kq_2)	&
				 /Scalar(IJK, 1))
!				 	
!

	ELSE
	     Scalar_c (IJK, L) = ZERO
	     Scalar_p (IJK, L) = ONE
	 END IF  !for scalars having very small values 
!
! *F_GS(IJK,1) : Drag Coef. for later use, sof
!
!            Diffusion coefficient for User-defined Scalars
!            
             Dif_Scalar(IJK, L) =EP_g(IJK)*MU_GT(IJK) /Sigma_E
!
	ELSE IF(L.EQ.3) THEN
!
! This equation (L=3) is not defined if Solids and gas turbulence are zero.
!	
	IF (Ep_s(IJK,M) .GT. SMALL_NUMBER .AND. 			&
	    SCALAR(IJK,1).GT. SMALL_NUMBER .AND.  			&
	    SCALAR(IJK,2) .GT. SMALL_NUMBER.AND.  			&
	    SCALAR(IJK,3) .GT. SMALL_NUMBER) THEN	
!
! PI_q_2: interaction term accounting for the continuous phase influence
! As defined by Simonin (Eq. 32).
	  PI_q_2 = EP_s(IJK,M)*RO_s(M)/Tau_12_st*(K_12(IJK)-2.d0*K_2)
!
	  IF(PI_q_2 .GE. ZERO) THEN
	    Pos_PI_q_2 = PI_q_2
	    Neg_PI_q_2 = ZERO
	  ELSE
	    Pos_PI_q_2 = ZERO
	    Neg_PI_q_2 = PI_q_2
!
	  END IF	
!
	    K_2 = SCALAR(IJK,3)
!
             Scalar_c (IJK, L) = (EP_s(IJK,M) * Pos_Tauij_sDUi_sODxj +	&
	     	                  Pos_PI_q_2 )
!
             Scalar_p (IJK, L) =(-EP_s(IJK,M) * Neg_Tauij_sDUi_sODxj/	&
	     			Scalar(IJK, 3) +			&
				EP_s(IJK,M)*RO_s(M)*Diss_2 /K_2		&
				-Neg_PI_q_2/K_2)	
!		  
!            Diffusion coefficient for User-defined Scalar 3
!            
	     K_2_T = (MU_2_T/Sigma_k+EP_s(IJK,M)*RO_s(M)*10.d0/27.d0*	&
	     	      Tau_12_st*K_2)/	 				&
	             (EP_s(IJK,M)*RO_s(M)*				&
		     (ONE + 5.d0/9.d0*Tau_12_st*zeta_c/Tau_2_c))
	     
             Dif_Scalar(IJK, L) = K_2_T
!
	 !WRITE(71,3) (Scalar_c (IJK, L),Scalar_p (IJK, L)*K_2,   &
	            !  Dif_Scalar(IJK, L),MU_2_T )
	ELSE
	     Scalar_c (IJK, L) = ZERO
	     Scalar_p (IJK, L) = ONE
	     Dif_Scalar(IJK, L) = ZERO
	 END IF  !for scalars having very small values
	 RES_TIME = (INT((TIME + 0.1*DT)/RES_DT) + 1)*RES_DT
	 !WRITE(71,3) (Mu_s(IJK,M), Mu_2_T*Ep_s(IJK,M), Ep_s(IJK,M))
!DIL_EP_S
!
	ELSE IF(L.EQ.4) THEN
!
! This equation (L=4) is not defined if Solids and gas turbulence are zero.
!	
	IF (Ep_s(IJK,M) .GT. SMALL_NUMBER .AND. 			&
	    SCALAR(IJK,1).GT. SMALL_NUMBER .AND.  			&
	    SCALAR(IJK,2) .GT. SMALL_NUMBER.AND.  			&
	    SCALAR(IJK,3) .GT. SMALL_NUMBER.AND.  			&
	    SCALAR(IJK,4) .GT. SMALL_NUMBER) THEN
!
! PI_q_12: interaction term accounting for the continuous phase influence
! As defined by Simonin (Eq. 27).
	  PI_q_12 = -EP_s(IJK,M)*RO_s(M)/Tau_12_st*( (ONE+X_21)		&
	  	    *SCALAR(IJK,4)-2.d0*SCALAR(IJK,1)-2.d0*X_21*	&
		    SCALAR(IJK,3))
!
	  IF(PI_q_12 .GE. ZERO) THEN
	    Pos_PI_q_12 = PI_q_12
	    Neg_PI_q_12 = ZERO
	  ELSE
	    Pos_PI_q_12 = ZERO
	    Neg_PI_q_12 = PI_q_12
!
	  END IF	
!!
             Scalar_c (IJK, L) = (Pos_Tau_gi_sj+Pos_PI_q_12)
!
             Scalar_p (IJK, L) = (EP_s(IJK,M)*RO_s(M)/Tau_12-		&
	     			(Neg_PI_q_12+Neg_Tau_gi_sj)		&
	     			/SCALAR(IJK,4))
!		  
!            Diffusion coefficient for User-defined Scalar 3
!            
	     K_2_T = Ep_s(IJK,M)*RO_s(M)/3.d0*SCALAR(IJK,4)*Tau_12/Sigma_k
	 ELSE
	     Scalar_c (IJK, L) = ZERO
	     Scalar_p (IJK, L) = ONE
	     Dif_Scalar(IJK, L) = ZERO
	 END IF  !for scalars having very small values	
	 !WRITE(71,3) (PI_q_12,Tau_gi_sj,EP_s(IJK,M)*RO_s(M)/Tau_12*SCALAR(IJK,4))
	END IF ! For Scalars with ID 1, 2 ,3 and 4.	     
           END DO 
!
         ENDIF 
      END DO
!OK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!\\Sendrecv operations - just to make sure all the variables computed are
!  are passed and updated locally - fool-proof approach - Sreekanth - 102199

       call send_recv(Scalar_c,4)
       call send_recv(Scalar_p,4)
       call send_recv(Dif_Scalar,4)
      RETURN  
      END SUBROUTINE SCALAR_PROP 
