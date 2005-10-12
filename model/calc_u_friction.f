!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_U_FRICTION(IJK1, IJK2, FCELL, COM, L, M, Gw, Hw,  C
!                               Cw)                                    C
!                                                                      C
!  Purpose: Calculate hw and cw for kinetic theory and frictional      C
!           boundary conditions                                        C
!                                                                      C
!                                                                      C
!  Author: Anuj Srivastava, Princeton University      Date: 12-APR-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_U_FRICTION(IJK1, IJK2, FCELL, COM, L, M, Gw,&
                                 Hw, Cw)
!
      USE param 
      USE param1 
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE turb
      USE visc_s
      USE geometry
      USE indices
      USE bc
      USE compar
      IMPLICIT NONE
!
!  Function subroutines
!
!
!  Local variables
!
!
!                      IJK indices for wall cell and fluid cell
      INTEGER          IJK, IJK1, IJK2
!
!                      Other indices
      INTEGER          IJK2E, IPJK2, IPJKM2, IPJKP2, IPJMK2, IPJPK2
      INTEGER          IJK2N, IJPK2, IJPKM2, IJPKP2, IMJPK2
      INTEGER          IJK2T, IJKP2, IJMKP2, IMJKP2
 
!
!                      Average scalars
      DOUBLE PRECISION EP_avg, EPg_avg, TH_avg, Mu_g_avg, RO_g_avg,Dp_avg
!
!                      Average velocities
      DOUBLE PRECISION WGC1, WGC2, WGCM, VGC1, VGC2, UGC1, UGC2
!
!                      The location (e,w,n...) of fluid cell
      CHARACTER        FCELL
!
!                      Velocity component or granular temp (U, V, W, T)
      CHARACTER        COM
!
!                      Solids phase index
      INTEGER          M, MM
!
!                      Coefficient of 1st term on LHS
      DOUBLE PRECISION Gw
!
!                      Coefficient of 2nd term on LHS
      DOUBLE PRECISION Hw
!
!                      Cw - all terms appearing on RHS
      DOUBLE PRECISION Cw
 
!                      del.u
      DOUBLE PRECISION DEL_DOT_U
!
!                      S:S
      DOUBLE PRECISION S_DDOT_S
!
!                      S_dd (d can be x, y or z)
      DOUBLE PRECISION S_dd
!
!                      values of U_sm, V_sm, W_sm at appropriate place
!                      on boundary wall
      DOUBLE PRECISION USCM, VSCM,WSCM
      DOUBLE PRECISION USCM1,USCM2,VSCM1,VSCM2,WSCM1,WSCM2
!
!                      values of U_g, V_g, W_g at appropriate place
!                      on boundary wall
      DOUBLE PRECISION UGC, VGC, WGC
      DOUBLE PRECISION USC1,USC2,VSC1,VSC2,WSC1,WSC2
!
!                      Magnitude of gas-solids relative velocity
      DOUBLE PRECISION VREL
 
!
!                      radial distribution function at contact
      DOUBLE PRECISION g0
!
!                      Sum of eps*G_0
      DOUBLE PRECISION g0EP_avg
 
!
!                      Error message
      CHARACTER*80     LINE
 
!                      Index corresponding to boundary condition
      INTEGER          L
      
!              Average Radial distribution function
      DOUBLE PRECISION g_0AVG
!
!                      Average Simonin and Ahmadi variables (sof)
      DOUBLE PRECISION K_12_avg, Tau_12_avg, Tau_1_avg
!
!              slip velocity between wall and particles for Jenkins bc (sof)
      DOUBLE PRECISION VSLIP
!
 
!  Function subroutines
 
!  Statement functions
!
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
 
 
      IF (BC_JJ_PS(L).EQ.3) THEN
 
	GW=1d0
	HW=0d0
	CW=0d0
 
      ELSEIF(BC_JJ_PS(L).EQ.2) THEN
 
        GW=0d0
        HW=1d0
 
	IF (COM.EQ.'U') THEN
          CW=BC_Uw_s(L, M)
	ELSEIF (COM.EQ.'V') THEN
          CW=BC_Vw_s(L, M)
        ELSE
	  CW=BC_Ww_s(L, M)
	ENDIF
 
      ELSE
!
!  Note:  EP_s, MU_g, and RO_g are undefined at IJK1 (wall cell).  Hence
!         IJK2 (fluid cell) is used in averages.
!
      IF     (COM .EQ. 'U')THEN
        IF(FCELL .EQ. 'N')THEN
          IJK2E = EAST_OF(IJK2)
          IPJK2 = IP_OF(IJK2)
          IPJMK2 = JM_OF(IPJK2)
 
          g0 = g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK1), M, M)
          EP_avg = AVG_X(EP_s(IJK2, M), EP_s(IJK2E, M), I_OF(IJK1))
	  EPg_avg = AVG_X(EP_g(IJK2), EP_g(IJK2E), I_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), MM, MM) &
	               * AVG_X(EP_s(IJK2, MM), EP_s(IJK2E, MM), I_OF(IJK2))
          ENDDO
	  
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Y(&
                AVG_X(Theta_m(IJK1, M), Theta_m(IPJMK2, M), I_OF(IJK1)),&
                AVG_X(Theta_m(IJK2, M), Theta_m(IPJK2, M), I_OF(IJK2)),&
                J_OF(IJK1))
          ELSE
	    TH_avg = AVG_X(THETA_M(IJK2,M), THETA_M(IJK2E,M), I_OF(IJK2))
          ENDIF
          Mu_g_avg = AVG_X(Mu_g(IJK2), Mu_g(IJK2E), I_OF(IJK2))
          RO_g_avg = AVG_X(RO_g(IJK2), RO_g(IJK2E), I_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_X(K_12(IJK2), K_12(IJK2E), I_OF(IJK2))	    
            Tau_12_avg = AVG_X(Tau_12(IJK2), Tau_12(IJK2E), I_OF(IJK2))	    
            Tau_1_avg = AVG_X(Tau_1(IJK2), Tau_1(IJK2E), I_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF

          DP_avg   = AVG_X(D_P(IJK2,M), D_P(IJK2E,M), I_OF(IJK2))
 
!         Calculate velocity components at i+1/2, j+1/2, k (relative to IJK1)
          UGC  = AVG_Y(U_g(IJK1), U_g(IJK2),J_OF(IJK1))
          VGC  = AVG_X(V_g(IJK1), V_g(IPJMK2),I_OF(IJK1))
          WGC1 = AVG_X(AVG_Z_T(W_g(KM_OF(IJK2)), W_g(IJK2)),&
                       AVG_Z_T(W_g(KM_OF(IPJK2)), W_g(IPJK2)),&
                       I_OF(IJK2))
          WGC2 = AVG_X(AVG_Z_T(W_g(KM_OF(IJK1)), W_g(IJK1)),&
                       AVG_Z_T(W_g(KM_OF(IPJMK2)), W_g(IPJMK2)),&
                       I_OF(IJK1))
          WGC  = AVG_Y(WGC2, WGC1, J_OF(IJK1))
          USCM = AVG_Y(U_s(IJK1,M), U_s(IJK2,M),J_OF(IJK1))
          VSCM = AVG_X(V_s(IJK1,M), V_s(IPJMK2,M),I_OF(IJK1))
          WSCM1= AVG_X(AVG_Z_T(W_s(KM_OF(IJK2),M),W_s(IJK2,M)),&
                       AVG_Z_T(W_s(KM_OF(IPJK2),M),W_s(IPJK2,M)),&
                       I_OF(IJK2))
          WSCM2= AVG_X(AVG_Z_T(W_s(KM_OF(IJK1),M),W_s(IJK1,M)),&
                       AVG_Z_T(W_s(KM_OF(IPJMK2),M),W_s(IPJMK2,M)),&
                       I_OF(IJK1))
          WSCM = AVG_Y(WSCM2, WSCM1, J_OF(IJK1))
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
          DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
 
        ELSEIF(FCELL .EQ. 'S')THEN
           IJK2E= EAST_OF(IJK2)
           IPJK2= IP_OF(IJK2)
           IPJPK2= JP_OF(IPJK2)
 
          g0 = g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), M, M)
          EP_avg = AVG_X(EP_s(IJK2, M), EP_s(IJK2E, M), I_OF(IJK2))
	  EPg_avg = AVG_X(EP_g(IJK2), EP_g(IJK2E), I_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), MM, MM) &
	               * AVG_X(EP_s(IJK2, MM), EP_s(IJK2E, MM), I_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Y(&
                  AVG_X(Theta_m(IJK2, M),Theta_m(IPJK2, M),I_OF(IJK2)),&
                  AVG_X(Theta_m(IJK1, M),Theta_m(IPJPK2, M),I_OF(IJK1)),&
                  J_OF(IJK2))
          ELSE
	    TH_avg =  AVG_X(THETA_M(IJK2,M), THETA_M(IJK2E,M), I_OF(IJK2))
	  ENDIF
          Mu_g_avg = AVG_X(Mu_g(IJK2), Mu_g(IJK2E), I_OF(IJK2))
 
          RO_g_avg = AVG_X(RO_g(IJK2), RO_g(IJK2E), I_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_X(K_12(IJK2), K_12(IJK2E), I_OF(IJK2))	    
            Tau_12_avg = AVG_X(Tau_12(IJK2), Tau_12(IJK2E), I_OF(IJK2))	    
            Tau_1_avg = AVG_X(Tau_1(IJK2), Tau_1(IJK2E), I_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF

          DP_avg   = AVG_X(D_P(IJK2,M), D_P(IJK2E,M), I_OF(IJK2))
 
!         Calculate velocity components at i+1/2, j+1/2, k relative to IJK2
          UGC  = AVG_Y(U_g(IJK2),U_g(IJK1),J_OF(IJK2))
          VGC  = AVG_X(V_g(IJK2),V_g(IPJK2),I_OF(IJK2))
          WGC1 = AVG_X(AVG_Z_T(W_g(KM_OF(IJK2)), W_g(IJK2)),&
                       AVG_Z_T(W_g(KM_OF(IPJK2)), W_g(IPJK2)),&
                       I_OF(IJK2))
          WGC2 = AVG_X(AVG_Z_T(W_g(KM_OF(IJK1)), W_g(IJK1)),&
                       AVG_Z_T(W_g(KM_OF(IPJPK2)), W_g(IPJPK2)),&
                       I_OF(IJK1))
          WGC  = AVG_Y(WGC1, WGC2, J_OF(IJK2))
          USCM = AVG_Y(U_s(IJK2, M),U_s(IJK1, M),J_OF(IJK2))
          VSCM = AVG_X(V_s(IJK2, M),V_s(IPJK2, M),I_OF(IJK2))
          WSCM1= AVG_X(AVG_Z_T(W_s(KM_OF(IJK2),M),W_s(IJK2,M)),&
                       AVG_Z_T(W_s(KM_OF(IPJK2),M),W_s(IPJK2,M)),&
                       I_OF(IJK2))
          WSCM2= AVG_X(AVG_Z_T(W_s(KM_OF(IJK1),M),W_s(IJK1,M)),&
                       AVG_Z_T(W_s(KM_OF(IPJPK2),M),W_s(IPJPK2,M)),&
                       I_OF(IJK1))
          WSCM = AVG_Y(WSCM1, WSCM2, J_OF(IJK2))
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
           DSQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSEIF(FCELL .EQ. 'T')THEN
           IJK2E= EAST_OF(IJK2)
           IPJK2= IP_OF(IJK2)
           IPJKM2= KM_OF(IPJK2)
          g0 = g_0AVG(IJK2, IJK2E, 'X',I_OF(IJK2), M, M)
          EP_avg =&
                AVG_X(EP_s(IJK2, M), EP_s(IJK2E, M),I_OF(IJK2))
	  EPg_avg = AVG_X(EP_g(IJK2), EP_g(IJK2E), I_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), MM, MM) &
	               * AVG_X(EP_s(IJK2, MM), EP_s(IJK2E, MM), I_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Z(&
                  AVG_X(Theta_m(IJK1,M),Theta_m(IPJKM2, M),I_OF(IJK1)),&
                  AVG_X(Theta_m(IJK2,M),Theta_m(IPJK2, M),I_OF(IJK2)),&
                  K_OF(IJK1))
          ELSE
	    TH_avg = AVG_X(THETA_M(IJK2,M), THETA_M(IJK2E,M), I_OF(IJK2))
	  ENDIF
          Mu_g_avg =&
                AVG_X(Mu_g(IJK2), Mu_g(IJK2E), I_OF(IJK2))
 
          RO_g_avg =&
                AVG_X(RO_g(IJK2), RO_g(IJK2E), I_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_X(K_12(IJK2), K_12(IJK2E), I_OF(IJK2))	    
            Tau_12_avg = AVG_X(Tau_12(IJK2), Tau_12(IJK2E), I_OF(IJK2))	    
            Tau_1_avg = AVG_X(Tau_1(IJK2), Tau_1(IJK2E), I_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF

          DP_avg  =& 
                AVG_X(D_P(IJK2,M), D_P(IJK2E,M), I_OF(IJK2))
 
!         Calculate velocity components at i+1/2,j,k-1/2 relative to IJK2
          UGC  = AVG_Z(U_g(IJK1), U_g(IJK2), K_OF(IJK1))
          VGC1 = AVG_X(AVG_Y_N(V_g(JM_OF(IJK2)),V_g(IJK2)),&
                       AVG_Y_N(V_g(JM_OF(IPJK2)),V_g(IPJK2)),&
                       I_OF(IJK2))
          VGC2 = AVG_X(AVG_Y_N(V_g(JM_OF(IJK1)),V_g(IJK1)),&
                       AVG_Y_N(V_g(JM_OF(IPJKM2)),V_g(IPJKM2)),&
                       I_OF(IJK1))
          VGC  = AVG_Z(VGC2,VGC1,K_OF(IJK1))
          WGC  = AVG_X(W_g(IJK1), W_g(IPJKM2),I_OF(IJK1))
          USCM = AVG_Z(U_s(IJK1,M), U_s(IJK2,M), K_OF(IJK1))
          VSCM1= AVG_X(AVG_Y_N(V_s(JM_OF(IJK2),M),V_s(IJK2,M)),&
                       AVG_Y_N(V_s(JM_OF(IPJK2),M),V_s(IPJK2,M)),&
                       I_OF(IJK2))
          VSCM2= AVG_X(AVG_Y_N(V_s(JM_OF(IJK1),M),V_s(IJK1,M)),&
                       AVG_Y_N(V_s(JM_OF(IPJKM2),M),V_s(IPJKM2,M)),&
                       I_OF(IJK1))
          VSCM  = AVG_Z(VSCM2,VSCM1,K_OF(IJK1))
          WSCM = AVG_X(W_s(IJK1,M), W_s(IPJKM2,M), I_OF(IJK1))
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
           DSQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSEIF(FCELL .EQ. 'B')THEN
           IJK2E= EAST_OF(IJK2)
           IPJK2= IP_OF(IJK2)
           IPJKP2= KP_OF(IPJK2)
          g0 = g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), M, M)
          EP_avg =&
                AVG_X(EP_s(IJK2, M), EP_s(IJK2E, M), I_OF(IJK2))
	  EPg_avg = AVG_X(EP_g(IJK2), EP_g(IJK2E), I_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), MM, MM) &
	               * AVG_X(EP_s(IJK2, MM), EP_s(IJK2E, MM), I_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Z(&
                  AVG_X(Theta_m(IJK2,M), Theta_m(IPJK2,M),I_OF(IJK2)),&
                  AVG_X(Theta_m(IJK1,M), Theta_m(IPJKP2,M),I_OF(IJK1)),&
                  K_OF(IJK2))
          ELSE
	    TH_avg = AVG_X(THETA_M(IJK2,M), THETA_M(IJK2E,M), I_OF(IJK2))
	  ENDIF
          Mu_g_avg =&
                AVG_X(Mu_g(IJK2), Mu_g(IJK2E), I_OF(IJK2))
 
          RO_g_avg =&
                AVG_X(RO_g(IJK2), RO_g(IJK2E), I_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_X(K_12(IJK2), K_12(IJK2E), I_OF(IJK2))	    
            Tau_12_avg = AVG_X(Tau_12(IJK2), Tau_12(IJK2E), I_OF(IJK2))	    
            Tau_1_avg = AVG_X(Tau_1(IJK2), Tau_1(IJK2E), I_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF

          DP_avg = AVG_X(D_P(IJK2,M), D_P(IJK2E,M), I_OF(IJK2))
 
!         Calculate velocity components at i+1/2,j,k-1/2 relative to IJK1
          UGC  = AVG_Z(U_g(IJK2), U_g(IJK1), K_OF(IJK2))
          VGC1 = AVG_X(AVG_Y_N(V_g(JM_OF(IJK2)),V_g(IJK2)),&
                       AVG_Y_N(V_g(JM_OF(IPJK2)),V_g(IPJK2)),&
                       I_OF(IJK2))
          VGC2 = AVG_X(AVG_Y_N(V_g(JM_OF(IJK1)),V_g(IJK1)),&
                       AVG_Y_N(V_g(JM_OF(IPJKP2)),V_g(IPJKP2)),&
                       I_OF(IJK1))
          VGC  = AVG_Z(VGC1, VGC2, K_OF(IJK2))
          WGC  = AVG_X(W_g(IJK2), W_g(IPJK2),I_OF(IJK2))
          USCM = AVG_Z(U_s(IJK2, M), U_s(IJK1, M), K_OF(IJK2))
          VSCM1= AVG_X(AVG_Y_N(V_s(JM_OF(IJK2),M),V_s(IJK2,M)),&
                       AVG_Y_N(V_s(JM_OF(IPJK2),M),V_s(IPJK2,M)),&
                       I_OF(IJK2))
          VSCM2= AVG_X(AVG_Y_N(V_s(JM_OF(IJK1),M),V_s(IJK1,M)),&
                       AVG_Y_N(V_s(JM_OF(IPJKP2),M),V_s(IPJKP2,M)),&
                       I_OF(IJK1))
          VSCM = AVG_Z(VSCM1, VSCM2, K_OF(IJK2))
          WSCM = AVG_X(W_s(IJK2, M), W_s(IPJK2, M),I_OF(IJK2))
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
            DSQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
         CALL WRITE_ERROR('CALC_U_FRICTION', LINE, 1)
        ENDIF
 
      ELSEIF (COM .EQ. 'V')THEN
        IF(FCELL .EQ. 'T')THEN
          IJK2N = NORTH_OF(IJK2)
          IJPK2 = JP_OF(IJK2)
          IJPKM2 = KM_OF(IJPK2)
          g0 = g_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), M, M)
          EP_avg =&
                  AVG_Y(EP_s(IJK2, M), EP_s(IJK2N, M), J_OF(IJK2))
	  EPg_avg = AVG_Y(EP_g(IJK2), EP_g(IJK2N), J_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), MM, MM) &
	               * AVG_Y(EP_s(IJK2, MM), EP_s(IJK2N, MM), J_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Y(&
             AVG_Z(Theta_m(IJK1, M), Theta_m(IJK2, M), K_OF(IJK1)),&
             AVG_Z(Theta_m(IJPKM2, M), Theta_m(IJPK2, M), K_OF(IJPKM2)),&
             J_OF(IJK1))
          ELSE
	    TH_avg = AVG_Y(THETA_M(IJK2,M), THETA_M(IJK2N,M), J_OF(IJK2))
          ENDIF
          Mu_g_avg =&
                     AVG_Y(Mu_g(IJK2), Mu_g(IJK2N), J_OF(IJK2))
 
          RO_g_avg =&
                     AVG_Y(RO_g(IJK2), RO_g(IJK2N), J_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_Y(K_12(IJK2), K_12(IJK2N), J_OF(IJK2))	    
            Tau_12_avg = AVG_Y(Tau_12(IJK2), Tau_12(IJK2N), J_OF(IJK2))	    
            Tau_1_avg = AVG_Y(Tau_1(IJK2), Tau_1(IJK2N), J_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF

          DP_avg   = AVG_Y(D_P(IJK2,M), D_P(IJK2N,M), J_OF(IJK2))
 
!         Calculate velocity components at i,j+1/2,k+1/2 (relative to IJK1)
          UGC1 = AVG_X_E(&
                         AVG_Y(U_g(IM_OF(IJK1)), U_g(IM_OF(IJPKM2)),&
                         J_OF(IM_OF(IJK1))),&
                         AVG_Y(U_g(IJK1), U_g(IJPKM2),J_OF(IJK1)),&
                         I_OF(IJK1))
          UGC2 = AVG_X_E(&
                         AVG_Y(U_g(IM_OF(IJK2)), U_g(IM_OF(IJPK2)),&
                         J_OF(IM_OF(IJK2))),&
                         AVG_Y(U_g(IJK2), U_g(IJPK2),J_OF(IJK2)),&
                         I_OF(IJK2))
          UGC  = AVG_Z(UGC1, UGC2, K_OF(IJK1))
          VGC  = AVG_Z(V_g(IJK1), V_g(IJK2),K_OF(IJK1))
          WGC  = AVG_Y(W_g(IJK1), W_g(IJPKM2), J_OF(IJK1))
          USCM1= AVG_X_E(&
                       AVG_Y(U_s(IM_OF(IJK1),M),U_s(IM_OF(IJPKM2),M),&
                       J_OF(IM_OF(IJK1))),&
                       AVG_Y(U_s(IJK1,M), U_s(IJPKM2,M),J_OF(IJK1)),&
                       I_OF(IJK1))
          USCM2= AVG_X_E(&
                       AVG_Y(U_s(IM_OF(IJK2),M),U_s(IM_OF(IJPK2),M),&
                       J_OF(IM_OF(IJK2))),&
                       AVG_Y(U_s(IJK2,M), U_s(IJPK2,M),J_OF(IJK2)),&
                       I_OF(IJK2))
          USCM = AVG_Z(USCM1, USCM2, K_OF(IJK1))
          VSCM = AVG_Z(V_s(IJK1,M), V_s(IJK2,M),K_OF(IJK1))
          WSCM = AVG_Y(W_s(IJK1,M), W_s(IJPKM2,M), J_OF(IJK1))
 
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
          DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSEIF(FCELL .EQ. 'B')THEN
          IJK2N = NORTH_OF(IJK2)
          IJPK2 = JP_OF(IJK2)
          IJPKP2 = KP_OF(IJPK2)
          g0 = g_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), M, M)
          EP_avg =&
                   AVG_Y(EP_s(IJK2, M), EP_s(IJK2N, M), J_OF(IJK2))
	  EPg_avg = AVG_Y(EP_g(IJK2), EP_g(IJK2N), J_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), MM, MM) &
	               * AVG_Y(EP_s(IJK2, MM), EP_s(IJK2N, MM), J_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Y(&
                AVG_Z(Theta_m(IJK2, M), Theta_m(IJK1, M), K_OF(IJK2)),&
                AVG_Z(Theta_m(IJPK2, M), Theta_m(IJPKP2, M), K_OF(IJPK2)),&
                J_OF(IJK2))
          ELSE
	    TH_avg = AVG_Y(THETA_M(IJK2,M), THETA_M(IJK2N,M), J_OF(IJK2))
          ENDIF
          Mu_g_avg =&
                     AVG_Y(Mu_g(IJK2), Mu_g(IJK2N), J_OF(IJK2))
! chang by rong 
          RO_g_avg =&
                     AVG_Z(RO_g(IJK2), RO_g(IJK2N), J_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_Y(K_12(IJK2), K_12(IJK2N), J_OF(IJK2))	    
            Tau_12_avg = AVG_Y(Tau_12(IJK2), Tau_12(IJK2N), J_OF(IJK2))	    
            Tau_1_avg = AVG_Y(Tau_1(IJK2), Tau_1(IJK2N), J_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF

          DP_avg   = AVG_Y(D_P(IJK2,M), D_P(IJK2N,M), J_OF(IJK2))
! change by rong
 
!         Calculate velocity components at i,j+1/2,k+1/2 (relative to IJK2)
          UGC1 = AVG_X_E(&
                         AVG_Y(U_g(IM_OF(IJK1)), U_g(IM_OF(IJPKP2)),&
                         J_OF(IM_OF(IJK1))),&
                         AVG_Y(U_g(IJK1), U_g(IJPKP2),J_OF(IJK1)),&
                         I_OF(IJK1))
          UGC2 = AVG_X_E(&
                         AVG_Y(U_g(IM_OF(IJK2)), U_g(IM_OF(IJPK2)),&
                         J_OF(IM_OF(IJK2))),&
                         AVG_Y(U_g(IJK2), U_g(IJPK2),J_OF(IJK2)),&
                         I_OF(IJK2))
          UGC  = AVG_Z(UGC2, UGC1, K_OF(IJK2))
          VGC  = AVG_Z(V_g(IJK2), V_g(IJK1),K_OF(IJK2))
          WGC  = AVG_Y(W_g(IJK2), W_g(IJPK2), J_OF(IJK2))
          USCM1= AVG_X_E(&
                       AVG_Y(U_s(IM_OF(IJK1),M),U_s(IM_OF(IJPKP2),M),&
                       J_OF(IM_OF(IJK1))),&
                       AVG_Y(U_s(IJK1,M), U_s(IJPKP2,M),J_OF(IJK1)),&
                       I_OF(IJK1))
          USCM2= AVG_X_E(&
                       AVG_Y(U_s(IM_OF(IJK2),M),U_s(IM_OF(IJPK2),M),&
                       J_OF(IM_OF(IJK2))),&
                       AVG_Y(U_s(IJK2,M), U_s(IJPK2,M),J_OF(IJK2)),&
                       I_OF(IJK2))
          USCM = AVG_Z(USCM2, USCM1, K_OF(IJK2))
          VSCM = AVG_Z(V_s(IJK2,M), V_s(IJK1,M),K_OF(IJK2))
          WSCM = AVG_Y(W_s(IJK2,M), W_s(IJPK2,M), J_OF(IJK2))
 
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
           DSQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSEIF(FCELL .EQ. 'E')THEN
           IJK2N= NORTH_OF(IJK2)
           IJPK2= JP_OF(IJK2)
           IMJPK2= IM_OF(IJPK2)
          g0 = g_0AVG(IJK2, IJK2N, 'Y',J_OF(IJK2), M, M)
          EP_avg =&
                AVG_Y(EP_s(IJK2, M), EP_s(IJK2N, M),J_OF(IJK2))
	  EPg_avg = AVG_Y(EP_g(IJK2), EP_g(IJK2N), J_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), MM, MM) &
	               * AVG_Y(EP_s(IJK2, MM), EP_s(IJK2N, MM), J_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Y(&
                  AVG_X(Theta_m(IJK1,M),Theta_m(IJK2, M),I_OF(IJK1)),&
                  AVG_X(Theta_m(IMJPK2,M),Theta_m(IJPK2, M),I_OF(IMJPK2)),&
                  J_OF(IJK1))
          ELSE
	    TH_avg = AVG_Y(THETA_M(IJK2,M), THETA_M(IJK2N,M), J_OF(IJK2))
	  ENDIF
          Mu_g_avg =&
                AVG_Y(Mu_g(IJK2), Mu_g(IJK2N), J_OF(IJK2))
 
          RO_g_avg =&
                AVG_Y(RO_g(IJK2), RO_g(IJK2N), J_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_Y(K_12(IJK2), K_12(IJK2N), J_OF(IJK2))	    
            Tau_12_avg = AVG_Y(Tau_12(IJK2), Tau_12(IJK2N), J_OF(IJK2))	    
            Tau_1_avg = AVG_Y(Tau_1(IJK2), Tau_1(IJK2N), J_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF

          DP_avg = AVG_Y(D_P(IJK2,M), D_P(IJK2N,M), J_OF(IJK2))
 
!         Calculate velocity components at i+1/2,j+1/2,k relative to IJK1
          UGC  = AVG_Y(U_g(IJK1), U_g(IMJPK2), J_OF(IJK1))
          VGC  = AVG_X(V_g(IJK1), V_g(IJK2), I_OF(IJK1))
          WGC1 = AVG_Y(AVG_Z_T(W_g(KM_OF(IJK1)), W_g(IJK1)),&
                       AVG_Z_T(W_g(KM_OF(IMJPK2)), W_g(IMJPK2)),&
                       J_OF(IJK1))
          WGC2 = AVG_Y(AVG_Z_T(W_g(KM_OF(IJK2)), W_g(IJK2)),&
                       AVG_Z_T(W_g(KM_OF(IJPK2)), W_g(IJPK2)),&
                       J_OF(IJK2))
          WGC  = AVG_X(WGC1, WGC2, I_OF(IJK1))
          USCM = AVG_Y(U_s(IJK1,M), U_s(IMJPK2,M), J_OF(IJK1))
          VSCM = AVG_X(V_s(IJK1, M), V_s(IJK2, M), I_OF(IJK1))
          WSCM1= AVG_Y(AVG_Z_T(W_s(KM_OF(IJK1),M), W_s(IJK1,M)),&
                       AVG_Z_T(W_s(KM_OF(IMJPK2),M), W_s(IMJPK2,M)),&
                       J_OF(IJK1))
          WSCM2 = AVG_Y(AVG_Z_T(W_s(KM_OF(IJK2),M), W_s(IJK2,M)),&
                       AVG_Z_T(W_s(KM_OF(IJPK2),M), W_s(IJPK2,M)),&
                       J_OF(IJK2))
          WSCM  = AVG_X(WSCM1, WSCM2, I_OF(IJK1))
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
           DSQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSEIF(FCELL .EQ. 'W')THEN
           IJK2N= NORTH_OF(IJK2)
           IJPK2= JP_OF(IJK2)
           IPJPK2= IP_OF(IJPK2)
 
          g0 = g_0AVG(IJK2, IJK2N, 'Y',J_OF(IJK2), M, M)
          EP_avg =&
                AVG_Y(EP_s(IJK2, M), EP_s(IJK2N, M),J_OF(IJK2))
	  EPg_avg = AVG_Y(EP_g(IJK2), EP_g(IJK2N), J_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), MM, MM) &
	               * AVG_Y(EP_s(IJK2, MM), EP_s(IJK2N, MM), J_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	     TH_avg = AVG_Y(&
                  AVG_X(Theta_m(IJK2,M),Theta_m(IJK1, M),I_OF(IJK2)),&
                  AVG_X(Theta_m(IJPK2,M),Theta_m(IPJPK2, M),I_OF(IJPK2)),&
                  J_OF(IJK2))
          ELSE
	    TH_avg = AVG_Y(THETA_M(IJK2,M), THETA_M(IJK2N,M), J_OF(IJK2))
	  ENDIF
          Mu_g_avg =&
                AVG_Y(Mu_g(IJK2), Mu_g(IJK2N), J_OF(IJK2))
 
          RO_g_avg =&
                AVG_Y(RO_g(IJK2), RO_g(IJK2N), J_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_Y(K_12(IJK2), K_12(IJK2N), J_OF(IJK2))	    
            Tau_12_avg = AVG_Y(Tau_12(IJK2), Tau_12(IJK2N), J_OF(IJK2))	    
            Tau_1_avg = AVG_Y(Tau_1(IJK2), Tau_1(IJK2N), J_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF

          DP_avg = AVG_Y(D_P(IJK2,M), D_P(IJK2N,M), J_OF(IJK2))
 
!         Calculate velocity components at i+1/2,j+1/2,k relative to IJK2
          UGC  = AVG_Y(U_g(IJK2), U_g(IJPK2), J_OF(IJK2))
          VGC  = AVG_X(V_g(IJK2), V_g(IJK1), I_OF(IJK2))
          WGC1 = AVG_Y(AVG_Z_T(W_g(KM_OF(IJK1)), W_g(IJK1)),&
                       AVG_Z_T(W_g(KM_OF(IPJPK2)), W_g(IPJPK2)),&
                       J_OF(IJK1))
          WGC2 = AVG_Y(AVG_Z_T(W_g(KM_OF(IJK2)), W_g(IJK2)),&
                       AVG_Z_T(W_g(KM_OF(IJPK2)), W_g(IJPK2)),&
                       J_OF(IJK2))
          WGC  = AVG_X(WGC2, WGC1, I_OF(IJK2))
          USCM = AVG_Y(U_s(IJK2,M), U_s(IJPK2,M), J_OF(IJK2))
          VSCM = AVG_X(V_s(IJK2, M), V_s(IJK1, M), I_OF(IJK2))
          WSCM1= AVG_Y(AVG_Z_T(W_s(KM_OF(IJK1),M), W_s(IJK1,M)),&
                       AVG_Z_T(W_s(KM_OF(IPJPK2),M), W_s(IPJPK2,M)),&
                       J_OF(IJK1))
          WSCM2 = AVG_Y(AVG_Z_T(W_s(KM_OF(IJK2),M), W_s(IJK2,M)),&
                       AVG_Z_T(W_s(KM_OF(IJPK2),M), W_s(IJPK2,M)),&
                       J_OF(IJK2))
          WSCM  = AVG_X(WSCM2, WSCM1, I_OF(IJK2))
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
            DSQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
         CALL WRITE_ERROR('CALC_U_FRICTION', LINE, 1)
        ENDIF
 
      ELSEIF (COM .EQ. 'W')THEN
        IF(FCELL .EQ. 'N')THEN
          IJK2T = TOP_OF(IJK2)
          IJKP2 = KP_OF(IJK2)
          IJMKP2 = JM_OF(IJKP2)
          g0 = g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, M)
          EP_avg =&
                   AVG_Z(EP_s(IJK2, M), EP_s(IJK2T, M), K_OF(IJK2))
	  EPg_avg = AVG_Z(EP_g(IJK2), EP_g(IJK2T), K_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), MM, MM) &
	               * AVG_Z(EP_s(IJK2, MM), EP_s(IJK2T, MM), K_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Y(&
                AVG_Z(Theta_m(IJK1, M), Theta_m(IJMKP2, M), K_OF(IJK1)),&
                AVG_Z(Theta_m(IJK2, M), Theta_m(IJKP2, M), K_OF(IJK2)),&
                J_OF(IJK1))
          ELSE
	    TH_avg = AVG_Z(THETA_M(IJK2, M), THETA_M(IJK2T, M), K_OF(IJK2))
	  ENDIF
          Mu_g_avg =&
                     AVG_Z(Mu_g(IJK2), Mu_g(IJK2T), K_OF(IJK2))
 
          RO_g_avg =&
                     AVG_Z(RO_g(IJK2), RO_g(IJK2T), K_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_Z(K_12(IJK2), K_12(IJK2T), K_OF(IJK2))	    
            Tau_12_avg = AVG_Z(Tau_12(IJK2), Tau_12(IJK2T), K_OF(IJK2))	    
            Tau_1_avg = AVG_Z(Tau_1(IJK2), Tau_1(IJK2T), K_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF

           DP_avg   = AVG_Z(D_P(IJK2,M), D_P(IJK2T,M), K_OF(IJK2))

!         Calculate velocity components at i,j+1/2,k+1/2 (relative to IJK1)
          UGC1 = AVG_X_E(&
                         AVG_Z(U_g(IM_OF(IJK1)), U_g(IM_OF(IJMKP2)),&
                         K_OF(IM_OF(IJK1)) ),&
                         AVG_Z(U_g(IJK1), U_g(IJMKP2), K_OF(IJK1)),&
                         I_OF(IJK1))
          UGC2 = AVG_X_E(&
                         AVG_Z(U_g(IM_OF(IJK2)), U_g(IM_OF(IJKP2)),&
                         K_OF(IM_OF(IJK2))),&
                         AVG_Z(U_g(IJK2), U_g(IJKP2), K_OF(IJK2)),&
                         I_OF(IJK2))
          UGC  = AVG_Y(UGC1, UGC2, J_OF(IJK1))
          VGC  = AVG_Z(V_g(IJK1), V_g(IJMKP2),K_OF(IJK1))
          WGC  = AVG_Y(W_g(IJK1), W_g(IJK2), J_OF(IJK1))
          USCM1= AVG_X_E(&
                       AVG_Z(U_s(IM_OF(IJK1),M),U_s(IM_OF(IJMKP2),M),&
                       K_OF(IM_OF(IJK1))),&
                       AVG_Z(U_s(IJK1,M), U_s(IJMKP2,M),K_OF(IJK1)),&
                       I_OF(IJK1))
          USCM2= AVG_X_E(&
                       AVG_Z(U_s(IM_OF(IJK2),M),U_s(IM_OF(IJKP2),M),&
                       K_OF(IM_OF(IJK2))),&
                       AVG_Z(U_s(IJK2,M), U_s(IJKP2,M),K_OF(IJK2)),&
                       I_OF(IJK2))
          USCM = AVG_Y(USCM1, USCM2, J_OF(IJK1))
          VSCM = AVG_Z(V_s(IJK1,M), V_s(IJMKP2,M),K_OF(IJK1))
          WSCM = AVG_Y(W_s(IJK1,M), W_s(IJK2,M), J_OF(IJK1))
 
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
          DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSEIF(FCELL .EQ. 'S')THEN
          IJK2T = TOP_OF(IJK2)
          IJKP2 = KP_OF(IJK2)
          IJPKP2 = JP_OF(IJKP2)
          g0 = g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, M)
          EP_avg =&
                   AVG_Z(EP_s(IJK2, M), EP_s(IJK2T, M), K_OF(IJK2))
	  EPg_avg = AVG_Z(EP_g(IJK2), EP_g(IJK2T), K_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), MM, MM) &
	               * AVG_Z(EP_s(IJK2, MM), EP_s(IJK2T, MM), K_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Y(&
                AVG_Z(Theta_m(IJK2, M), Theta_m(IJKP2, M), K_OF(IJK2)),&
                AVG_Z(Theta_m(IJK1, M), Theta_m(IJPKP2, M), K_OF(IJK1)),&
                J_OF(IJK2))
          ELSE
	    TH_avg = AVG_Z(THETA_M(IJK2, M), THETA_M(IJK2T, M), K_OF(IJK2))
	  ENDIF
          Mu_g_avg =&
                     AVG_Z(Mu_g(IJK2), Mu_g(IJK2T), K_OF(IJK2))
 
          RO_g_avg =&
                     AVG_Z(RO_g(IJK2), RO_g(IJK2T), K_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_Z(K_12(IJK2), K_12(IJK2T), K_OF(IJK2))	    
            Tau_12_avg = AVG_Z(Tau_12(IJK2), Tau_12(IJK2T), K_OF(IJK2))	    
            Tau_1_avg = AVG_Z(Tau_1(IJK2), Tau_1(IJK2T), K_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF
 
          DP_avg   = AVG_Z(D_P(IJK2,M), D_P(IJK2T,M), K_OF(IJK2))

!         Calculate velocity components at i,j+1/2,k+1/2 (relative to IJK2)
          UGC1 = AVG_X_E(&
                         AVG_Z(U_g(IM_OF(IJK1)), U_g(IM_OF(IJPKP2)),&
                         K_OF(IM_OF(IJK1))),&
                         AVG_Z(U_g(IJK1), U_g(IJPKP2), K_OF(IJK1)),&
                         I_OF(IJK1))
          UGC2 = AVG_X_E(&
                         AVG_Z(U_g(IM_OF(IJK2)), U_g(IM_OF(IJKP2)),&
                         K_OF(IM_OF(IJK2))),&
                         AVG_Z(U_g(IJK2), U_g(IJKP2), K_OF(IJK2)),&
                         I_OF(IJK2))
          UGC  = AVG_Y(UGC2, UGC1, J_OF(IJK2))
          VGC  = AVG_Z(V_g(IJK2), V_g(IJKP2),K_OF(IJK2))
          WGC  = AVG_Y(W_g(IJK2), W_g(IJK1), J_OF(IJK2))
          USCM1= AVG_X_E(&
                       AVG_Z(U_s(IM_OF(IJK1),M),U_s(IM_OF(IJPKP2),M),&
                       K_OF(IM_OF(IJK1))),&
                       AVG_Z(U_s(IJK1,M), U_s(IJPKP2,M),K_OF(IJK1)),&
                       I_OF(IJK1))
          USCM2= AVG_X_E(&
                       AVG_Z(U_s(IM_OF(IJK2),M),U_s(IM_OF(IJKP2),M),&
                       K_OF(IM_OF(IJK2))),&
                       AVG_Z(U_s(IJK2,M), U_s(IJKP2,M),K_OF(IJK2)),&
                       I_OF(IJK2))
          USCM = AVG_Y(USCM2, USCM1, J_OF(IJK2))
          VSCM = AVG_Z(V_s(IJK2,M), V_s(IJKP2,M),K_OF(IJK2))
          WSCM = AVG_Y(W_s(IJK2,M), W_s(IJK1,M), J_OF(IJK2))
 
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
          DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSEIF(FCELL .EQ. 'E')THEN
           IJK2T= TOP_OF(IJK2)
           IJKP2= KP_OF(IJK2)
           IMJKP2 = IM_OF(IJKP2)
          g0 = g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, M)
          EP_avg =&
                AVG_Z(EP_s(IJK2, M), EP_s(IJK2T, M), K_OF(IJK2))
	  EPg_avg = AVG_Z(EP_g(IJK2), EP_g(IJK2T), K_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), MM, MM) &
	               * AVG_Z(EP_s(IJK2, MM), EP_s(IJK2T, MM), K_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Z(&
                  AVG_X(Theta_m(IJK1,M),Theta_m(IJK2, M),I_OF(IJK1)),&
                  AVG_X(Theta_m(IMJKP2,M),Theta_m(IJKP2, M),I_OF(IMJKP2)),&
                  K_OF(IJK1))
          ELSE
	    TH_avg = AVG_Z(THETA_M(IJK2, M), THETA_M(IJK2T, M), K_OF(IJK2))
	  ENDIF
          Mu_g_avg =&
                AVG_Z(Mu_g(IJK2), Mu_g(IJK2T), K_OF(IJK2))
 
          RO_g_avg =&
                AVG_Z(RO_g(IJK2), RO_g(IJK2T), K_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_Z(K_12(IJK2), K_12(IJK2T), K_OF(IJK2))	    
            Tau_12_avg = AVG_Z(Tau_12(IJK2), Tau_12(IJK2T), K_OF(IJK2))	    
            Tau_1_avg = AVG_Z(Tau_1(IJK2), Tau_1(IJK2T), K_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF
 
          DP_avg = AVG_Z(D_P(IJK2,M), D_P(IJK2T,M), K_OF(IJK2))
 
!         Calculate velocity components at i+1/2,j,k+1/2 relative to IJK1
          UGC  = AVG_Z(U_g(IJK1), U_g(IMJKP2), K_OF(IJK1))
          VGC1 = AVG_Z(AVG_Y_N(V_g(JM_OF(IJK1)),V_g(IJK1)),&
                       AVG_Y_N(V_g(JM_OF(IMJKP2)),V_g(IMJKP2)),&
                       K_OF(IJK1))
          VGC2 = AVG_Z(AVG_Y_N(V_g(JM_OF(IJK2)),V_g(IJK2)),&
                       AVG_Y_N(V_g(JM_OF(IJKP2)),V_g(IJKP2)),&
                       K_OF(IJK2))
          VGC  = AVG_X(VGC1,VGC2,I_OF(IJK1))
          WGC  = AVG_X(W_g(IJK1), W_g(IJK2),I_OF(IJK1))
          USCM = AVG_Z(U_s(IJK1,M), U_s(IMJKP2,M), K_OF(IJK1))
          VSCM1= AVG_Z(AVG_Y_N(V_s(JM_OF(IJK1),M),V_s(IJK1,M)),&
                       AVG_Y_N(V_s(JM_OF(IMJKP2),M),V_s(IMJKP2,M)),&
                       K_OF(IJK1))
          VSCM2= AVG_Z(AVG_Y_N(V_s(JM_OF(IJK2),M),V_s(IJK2,M)),&
                       AVG_Y_N(V_s(JM_OF(IJKP2),M),V_s(IJKP2,M)),&
                       K_OF(IJK2))
          VSCM  = AVG_X(VSCM1,VSCM2,I_OF(IJK1))
          WSCM = AVG_X(W_s(IJK1,M), W_s(IJK2,M), I_OF(IJK1))
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
           DSQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSEIF(FCELL .EQ. 'W')THEN
           IJK2T= TOP_OF(IJK2)
           IJKP2= KP_OF(IJK2)
           IPJKP2= IP_OF(IJKP2)
          g0 = g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, M)
          EP_avg =&
                AVG_Z(EP_s(IJK2, M), EP_s(IJK2T, M), K_OF(IJK2))
	  EPg_avg = AVG_Z(EP_g(IJK2), EP_g(IJK2T), K_OF(IJK2))
          g0EP_avg = ZERO
	  DO MM = 1, MMAX
	    g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), MM, MM) &
	               * AVG_Z(EP_s(IJK2, MM), EP_s(IJK2T, MM), K_OF(IJK2))
          ENDDO
 
          IF(GRANULAR_ENERGY) THEN
	    TH_avg = AVG_Z(&
                  AVG_X(Theta_m(IJK2,M),Theta_m(IJK1, M),I_OF(IJK2)),&
                  AVG_X(Theta_m(IJKP2,M),Theta_m(IPJKP2, M),I_OF(IJKP2)),&
                  K_OF(IJK2))
          ELSE
	    TH_avg = AVG_Z(THETA_M(IJK2, M), THETA_M(IJK2T, M), K_OF(IJK2))
	  ENDIF
          Mu_g_avg =&
                AVG_Z(Mu_g(IJK2), Mu_g(IJK2T), K_OF(IJK2))
 
          RO_g_avg =&
                AVG_Z(RO_g(IJK2), RO_g(IJK2T), K_OF(IJK2))

	  IF(SIMONIN .OR. AHMADI) THEN
! added for Simonin and Ahmadi model (sof)
            K_12_avg = AVG_Z(K_12(IJK2), K_12(IJK2T), K_OF(IJK2))	    
            Tau_12_avg = AVG_Z(Tau_12(IJK2), Tau_12(IJK2T), K_OF(IJK2))	    
            Tau_1_avg = AVG_Z(Tau_1(IJK2), Tau_1(IJK2T), K_OF(IJK2))
	  ELSE
	    K_12_avg = ZERO    
            Tau_12_avg = ZERO	    
            Tau_1_avg = ZERO
	  ENDIF

          DP_avg= AVG_Z(D_P(IJK2,M), D_P(IJKP2,M), K_OF(IJK2))
 
!         Calculate velocity components at i+1/2,j,k+1/2 relative to IJK2
          UGC  = AVG_Z(U_g(IJK2), U_g(IJKP2), K_OF(IJK2))
          VGC1 = AVG_Z(AVG_Y_N(V_g(JM_OF(IJK1)),V_g(IJK1)),&
                       AVG_Y_N(V_g(JM_OF(IPJKP2)),V_g(IPJKP2)),&
                       K_OF(IJK1))
          VGC2 = AVG_Z(AVG_Y_N(V_g(JM_OF(IJK2)),V_g(IJK2)),&
                       AVG_Y_N(V_g(JM_OF(IJKP2)),V_g(IJKP2)),&
                       K_OF(IJK2))
          VGC  = AVG_X(VGC2,VGC1,I_OF(IJK2))
          WGC  = AVG_X(W_g(IJK2), W_g(IJK1),I_OF(IJK2))
          USCM = AVG_Z(U_s(IJK2,M), U_s(IJKP2,M), K_OF(IJK2))
          VSCM1= AVG_Z(AVG_Y_N(V_s(JM_OF(IJK1),M),V_s(IJK1,M)),&
                       AVG_Y_N(V_s(JM_OF(IPJKP2),M),V_s(IPJKP2,M)),&
                       K_OF(IJK1))
          VSCM2= AVG_Z(AVG_Y_N(V_s(JM_OF(IJK2),M),V_s(IJK2,M)),&
                       AVG_Y_N(V_s(JM_OF(IJKP2),M),V_s(IJKP2,M)),&
                       K_OF(IJK2))
          VSCM  = AVG_X(VSCM2,VSCM1,I_OF(IJK2))
          WSCM = AVG_X(W_s(IJK2,M), W_s(IJK1,M), I_OF(IJK2))
!
!         magnitude of gas-solids relative velocity
!
          VREL =&
           DSQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)
!
! slip velocity for use in Jenkins bc (sof)	  
	  VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
	                + (WSCM-BC_WW_S(L,M))**2 )
 
          CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
                             S_DDOT_S, S_dd)
 
          CALL CALC_Gw_Hw_Cw(g0, EP_avg,EPg_avg, g0EP_avg,TH_avg,Mu_g_avg,RO_g_avg, &
                             DP_avg, VREL, M, DEL_DOT_U, S_DDOT_S, S_dd, USCM,&
			     K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
                             BC_Uw_s(L, M), gw, hw,cw)
 
        ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
         CALL WRITE_ERROR('CALC_U_FRICTION', LINE, 1)
        ENDIF
 
      ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown COM'
         CALL WRITE_ERROR('CALC_U_FRICTION', LINE, 1)
      ENDIF
!
      ENDIF
 
      RETURN
      END
 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Gw_Hw_Cw(g0, EPS, EPG, g0EP_avg TH, Mu_g_avg,     C
!                            RO_g_avg, DP_avg, VREL, M,                C 
!                             DEL_U, S_S, S_dd, VEL, W_VEL,            C
!                             gw, hw,cw)                               C
!                       +K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP  (sof)  C
!                                                                      C
!  Purpose: Calculate Gw, Hw, and Cw                                   C
!                                                                      C
!  Author: A. Srivastava & K. Agrawal, Princeton Univ. Date: 10-APR-98 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: EPS, TH, C_e, RO_s                            C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:Dp_avg,RO_g_avg,Mu_g_avg                            C
!                                                                      C
!  Modified: Sofiane Benyahia, Fluent Inc.             Date: 03-FEB-05 C
!  Purpose: Include conductivity defined by Simonin and Ahmadi         C
!           Also included Jenkins small frictional limit               C
!                                                                      C
!  Literature/Document References: See calcmu_s.f for ref. on Simonin  C
!  and Ahmadi models; for Jenkins BC: Jenkins and Louge, Phys. fluids  C
!  9 (10), 2835. See equation (2) in the paper                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
 
      SUBROUTINE CALC_Gw_Hw_Cw(g0,EPS,EPG,g0EP_avg,TH,Mu_g_avg,RO_g_avg, &
                               DP_avg, VREL, M, DEL_U, S_S, S_dd, VEL, &
			       K_12_avg, Tau_12_avg, Tau_1_avg,VSLIP,&
			       W_VEL, gw, hw,cw)

      USE param 
      USE param1 
      USE constant
      USE physprop
      USE fldvar
      USE bc
      USE run
      USE mpi_utility
      IMPLICIT NONE
 
      INTEGER          M
!
!              Average solids and gas volume fraction
      DOUBLE PRECISION EPS, EPG
 
!              Average theta_m
      DOUBLE PRECISION TH
 
!              Coefficient of 2nd term
      DOUBLE PRECISION Hw
 
!              Coefficient of 1st term; part of this coefficient
      DOUBLE PRECISION Gw, Mu_s
 
!              Term appearing in coefficient of 2nd term
      DOUBLE PRECISION F_2
 
!              Right-hand side term
      DOUBLE PRECISION Cw
 
!              Relevant solids velocity at wall
      DOUBLE PRECISION VEL
 
!              Relevant wall velocity
      DOUBLE PRECISION W_VEL
 
!              del.u
      DOUBLE PRECISION DEL_U
 
!              S:S
      DOUBLE PRECISION S_S
 
!              S_dd (d can be x,y or z)
      DOUBLE PRECISION S_dd
 
!              Frictional  Pressure Pf
      DOUBLE PRECISION Pf
 
!              Critical pressure Pc
      DOUBLE PRECISION Pc
 
!              Chi (appears in frictional boundary condition)
      DOUBLE PRECISION Chi
 
!              Viscosity
      DOUBLE PRECISION Mu
 
!              Bulk viscosity
      DOUBLE PRECISION Mu_b
 
!              Magnitude of slip velocity
      DOUBLE PRECISION VREL
!
!              slip velocity between wall and particles for Jenkins bc (sof)
      DOUBLE PRECISION VSLIP
 
!              Average gas density
      DOUBLE PRECISION RO_g_avg
 
!              Average cross-correlation K_12 and time-scales (sof)
      DOUBLE PRECISION K_12_avg, Tau_12_avg, Tau_1_avg

! add by rong
!              Average solid diameter
      DOUBLE PRECISION dp_avg

! add by rong 
!              Average gas viscosity
      DOUBLE PRECISION Mu_g_avg
 
!              Reynolds number based on slip velocity
      DOUBLE PRECISION Re_g
 
!              Friction Factor in drag coefficient
      DOUBLE PRECISION C_d
 
!              Drag Coefficient
      DOUBLE PRECISION Beta, DgA
 
!              Square root of S:S or the form suggested by Savage
      DOUBLE PRECISION ZETA
 
!              Viscosity corrected for interstitial fluid effects
      DOUBLE PRECISION Mu_star
 
!              Radial distribution function
      DOUBLE PRECISION g0, G_0
!
!                      Sum of eps*G_0 (sof June 16 2005)
      DOUBLE PRECISION g0EP_avg
 
!              Constants in Simonin model
      DOUBLE PRECISION Sigma_c, Tau_2_c, Tau_12_st, Nu_t
      DOUBLE PRECISION Tau_2, zeta_c_2, MU_2_T_Kin, Mu_2_Col
 
!                      Error message
      CHARACTER*80     LINE
 
!              Other local terms
      DOUBLE PRECISION phin
!
! This is done here similar to bc_theta to avoid small negative values of
! Theta coming most probably from linear solver
      IF(TH .LE. ZERO)THEN
        TH = 1D-8

        if (myPE.eq.PE_IO) then   
	   WRITE(*,*)'Warning: Negative granular temp at wall set to 1e-8'
!          CALL WRITE_ERROR('THETA_HW_CW', LINE, 1)
        end if
      ENDIF
! 
      G_0 = g0
      
! modify F_2 if Jenkins BC is used (sof)    
 
      IF(JENKINS) THEN
!
        IF (VSLIP == ZERO) THEN
! if solids velocity field is initialized to zero, use free slip bc
	  F_2 = zero
!
	ELSE
          IF(AHMADI) THEN
! Ahmadi model uses different solids pressure model
!
! the coefficient mu in Jenkins paper is defined as tan_Phi_w, that's how
! I understand it from soil mechanic papers, i.e., G.I. Tardos, powder
! Tech. 92 (1997), 61-74. See his equation (1). Define Phi_w in mfix.dat!
!
            F_2 = tan_Phi_w*RO_s(M)*EPS* &
	          ((ONE + 4.0D0*g0EP_avg) + HALF*(ONE -C_e*C_e))*TH/VSLIP
!
! here F_2 divided by VSLIP to use the same bc as Johnson&Jackson
!
          ELSE  ! Simonin or granular models use same solids pressure
            F_2 = tan_Phi_w*RO_s(M)*EPS*(1d0+ 4.D0 * Eta *g0EP_avg)*TH/VSLIP
	  ENDIF !for Ahmadi
!
	ENDIF ! for vslip < small_number
!
      ELSE ! no change to the original code if Jenkins BC not used
 
        F_2 = (PHIP*DSQRT(3d0*TH)*Pi*RO_s(M)*EPS*G_0)&
              /(6d0*EPS_max)
!
      ENDIF !for Jenkins
 
      Mu = (5d0*DSQRT(Pi*TH)*Dp_avg*RO_s(M))/96d0
 
      Mu_b = (256d0*Mu*EPS*g0EP_avg)/(5d0*Pi)
 
      Re_g = EPG*RO_g_avg*Dp_avg*VREL/Mu_g_avg
      IF (Re_g.lt.1000d0) THEN
         C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
      ELSE
         C_d = 0.44d0
      ENDIF
      DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(Dp_avg*EPG**(2.65d0))
      IF(VREL == ZERO) DgA = LARGE_NUMBER
      Beta = SWITCH*EPS*DgA
!
! particle relaxation time
      Tau_12_st = RO_s(M)/(DgA+small_number)
 
!     SWITCH enables us to turn on/off the modification to the
!     particulate phase viscosity. If we want to simulate gas-particle
!     flow then SWITCH=1 to incorporate the effect of drag on the
!     particle viscosity. If we want to simulate granular flow
!     without the effects of an interstitial gas, SWITCH=0.
 
      IF(SWITCH == ZERO .OR. Ro_g_avg == ZERO)THEN
        Mu_star = Mu
		
      ELSEIF(TH .LT. SMALL_NUMBER)THEN
        MU_star = ZERO
	
      ELSE
	Mu_star = RO_S(M)*EPS* G_0*TH* Mu/ &
	         (RO_S(M)*g0EP_avg*TH + 2.0d0*SWITCH*DgA/RO_S(M)* Mu)
	
      ENDIF
 
      Mu_s = ((2d0+ALPHA)/3d0)*((Mu_star/(Eta*(2d0-Eta)*&
                   G_0))*(ONE+1.6d0*Eta*g0EP_avg&
                   )*(ONE+1.6d0*Eta*(3d0*Eta-2d0)*&
                   g0EP_avg)+(0.6d0*Mu_b*Eta))
 
      IF(SIMONIN) THEN !see calc_mu_s for explanation of these definitions
!
        Sigma_c = (ONE+ C_e)*(3.d0-C_e)/5.d0
        Tau_2_c = DP_avg/(6.d0*EPS*G_0*DSQRT(16.d0*(TH+Small_number)/PI))
	zeta_c_2= 2.D0/5.D0*(ONE+ C_e)*(3.d0*C_e-ONE)
	Nu_t =  Tau_12_avg/Tau_12_st
        Tau_2 = ONE/(2.D0/Tau_12_st+Sigma_c/Tau_2_c)
!
	MU_2_T_Kin = (2.0D0/3.0D0*K_12_avg*Nu_t + TH * &
                     (ONE+ zeta_c_2*EPS*G_0))*Tau_2
!
	Mu_2_Col = 8.D0/5.D0*EPS*G_0*Eta* (MU_2_T_Kin+ &
                   Dp_avg*DSQRT(TH/PI))
!
	Mu_s = EPS*RO_s(M)*(MU_2_T_Kin + Mu_2_Col)
!
      ELSE IF(AHMADI) THEN
!
	Mu_s = ONE/(ONE+ Tau_1_avg/Tau_12_st * (ONE-EPS/EPS_max)**3)&
	       *0.1045D0*(ONE/G_0+3.2D0*EPS+12.1824D0*G_0*EPS*EPS)  &
	       *Dp_avg*RO_s(M)* DSQRT(TH)
      ENDIF
 
!     Calculating frictional terms
 
      IF (EPS<= EPS_f_min) THEN
         Pf = ZERO
         Chi = ZERO
         ZETA = 1d0
 
      ELSE
 
         IF (SAVAGE.EQ.1) THEN    !form of Savage
            ZETA = ((48d0*Eta*(1d0-Eta)*RO_s(M)*EPS*EPS*G_0*&
                    (TH**1.5d0))/&
                    (SQRT_Pi*Dp_avg*2d0*Mu_s))**0.5d0
 
         ELSEIF (SAVAGE.EQ.0)  THEN !S:S form
            ZETA = DSQRT(S_S)
 
         ELSE
            ZETA = DSQRT(S_S + (TH/(Dp_avg*Dp_avg)))
         ENDIF
 
 
         IF (EPS> EPS_max) THEN
            Pc=  1d25*((EPS-EPS_max)**10d0)
         ELSE
            Pc = Fr*((EPS - EPS_f_min)**N_Pc)/((EPS_max - EPS +&
                  SMALL_NUMBER)**D_Pc)
         ENDIF
 
 
 
         IF ((DEL_U/(ZETA*N_Pf*DSQRT(2d0)*Sin_Phi)) .GT. 1d0) THEN
             Pf = ZERO
         ELSE
             Pf = Pc*(1d0 - (DEL_U/(ZETA*N_Pf*DSQRT(2d0)*Sin_Phi)))**&
                  (N_Pf-1d0)
         ENDIF
 
         Chi =&
                        DSQRT(2d0)*Pf*Sin_Phi*(N_Pf - (N_Pf-1d0)*&
                                      ((Pf/Pc)**(1d0/(N_Pf-1d0))))
 
 
 
         IF (Chi< ZERO) THEN
            Pf = Pc*((N_Pf/(N_Pf-1d0))**(N_Pf-1d0))
            Chi = ZERO
         ENDIF
 
      ENDIF
 
!     Calculating gw, hw, cw
 
      Gw = (MU_s + Chi/(2d0*ZETA))*DABS(VEL - W_VEL)
 
      Hw = F_2*DABS(VEL - W_VEL) + Pf*tan_Phi_w -&
           Chi*S_dd*tan_Phi_w/ZETA
 
      Cw = hw * W_VEL
 
      RETURN
      END
