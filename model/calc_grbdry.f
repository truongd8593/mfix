!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_GRBDRY(IJK1, IJK2, FCELL, COM, M, L, Hw)          C
!  Purpose: Calculate hw and cw for kinetic theory boundary conditions C
!                                                                      C
!  Author: K. Agrawal & A. Srivastava, Princeton Univ. Date: 19-JAN-98 C
!  Reviewer:                                           Date:           C
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
      SUBROUTINE CALC_GRBDRY(IJK1, IJK2, FCELL, COM, M, L, Hw)
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
      DOUBLE PRECISION EP_avg, TH_avg, Mu_g_avg, RO_g_avg
!
!                      Average Simonin and Ahmadi variables (sof)
      DOUBLE PRECISION K_12_avg, Tau_12_avg, Tau_1_avg
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
      INTEGER          M
!
!                      Wall momentum or granular energy coefficient
      DOUBLE PRECISION Hw
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
!              slip velocity between wall and particles for Jenkins bc (sof)
      DOUBLE PRECISION VSLIP
!
!                      radial distribution function at contact
      DOUBLE PRECISION g0
 
!
!                      Error message
      CHARACTER*80     LINE
!
!                      Index corresponding to boundary condition (sof)
      INTEGER          L
 
!              Average Radial distribution function
      DOUBLE PRECISION g_0AVG
!
!  Function subroutines
!
!
!
      DOUBLE PRECISION F_HW
 
!
!  Statement functions
!
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
!  Note:  EP_s, MU_g, and RO_g are undefined at IJK1 (wall cell).  Hence
!         IJK2 (fluid cell) is used in averages.
!
      IF     (COM .EQ. 'U')THEN
        IF(FCELL .EQ. 'N')THEN
          IJK2E = EAST_OF(IJK2)
          IPJK2 = IP_OF(IJK2)
          IPJMK2 = JM_OF(IPJK2)
 
          g0 = g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), M, M)
          EP_avg = AVG_X(EP_s(IJK2, M), EP_s(IJK2E, M), I_OF(IJK2))
          TH_avg = AVG_Y(&
              AVG_X(Theta_m(IJK1, M), Theta_m(IPJMK2, M), I_OF(IJK1)),&
              AVG_X(Theta_m(IJK2, M), Theta_m(IPJK2, M), I_OF(IJK2)),&
              J_OF(IJK1))
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
	  VSLIP= DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
 
        ELSEIF(FCELL .EQ. 'S')THEN
           IJK2E= EAST_OF(IJK2)
           IPJK2= IP_OF(IJK2)
           IPJPK2= JP_OF(IPJK2)
 
          g0 = g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), M, M)
          EP_avg = AVG_X(EP_s(IJK2, M), EP_s(IJK2E, M), I_OF(IJK2))
 
          TH_avg = AVG_Y(&
                AVG_X(Theta_m(IJK2, M),Theta_m(IPJK2, M),I_OF(IJK2)),&
                AVG_X(Theta_m(IJK1, M),Theta_m(IPJPK2, M),I_OF(IJK1)),&
                J_OF(IJK2))
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
	  VSLIP= DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
 
        ELSEIF(FCELL .EQ. 'T')THEN
           IJK2E= EAST_OF(IJK2)
           IPJK2= IP_OF(IJK2)
           IPJKM2= KM_OF(IPJK2)
          g0 = g_0AVG(IJK2, IJK2E, 'X',I_OF(IJK2), M, M)
          EP_avg =&
                AVG_X(EP_s(IJK2, M), EP_s(IJK2E, M),I_OF(IJK2))
 
          TH_avg = AVG_Z(&
                AVG_X(Theta_m(IJK1,M),Theta_m(IPJKM2, M),I_OF(IJK1)),&
                AVG_X(Theta_m(IJK2,M),Theta_m(IPJK2, M),I_OF(IJK2)),&
                K_OF(IJK1))
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
	  VSLIP= DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
 
        ELSEIF(FCELL .EQ. 'B')THEN
           IJK2E= EAST_OF(IJK2)
           IPJK2= IP_OF(IJK2)
           IPJKP2= KP_OF(IPJK2)
          g0 = g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), M, M)
          EP_avg =&
                AVG_X(EP_s(IJK2, M), EP_s(IJK2E, M), I_OF(IJK2))
 
          TH_avg = AVG_Z(&
                AVG_X(Theta_m(IJK2,M), Theta_m(IPJK2,M),I_OF(IJK2)),&
                AVG_X(Theta_m(IJK1,M), Theta_m(IPJKP2,M),I_OF(IJK1)),&
                K_OF(IJK2))
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
	  VSLIP= DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
        ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
         CALL WRITE_ERROR('CALC_GRBDRY', LINE, 1)
        ENDIF
 
      ELSEIF (COM .EQ. 'V')THEN
        IF(FCELL .EQ. 'T')THEN
          IJK2N = NORTH_OF(IJK2)
          IJPK2 = JP_OF(IJK2)
          IJPKM2 = KM_OF(IJPK2)
          g0 = g_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), M, M)
          EP_avg =&
                  AVG_Y(EP_s(IJK2, M), EP_s(IJK2N, M), J_OF(IJK2))
 
          TH_avg = AVG_Y(&
           AVG_Z(Theta_m(IJK1, M), Theta_m(IJK2, M), K_OF(IJK1)),&
           AVG_Z(Theta_m(IJPKM2, M), Theta_m(IJPK2, M), K_OF(IJPKM2)),&
           J_OF(IJK1))
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
	  VSLIP = DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
 
        ELSEIF(FCELL .EQ. 'B')THEN
          IJK2N = NORTH_OF(IJK2)
          IJPK2 = JP_OF(IJK2)
          IJPKP2 = KP_OF(IJPK2)
          g0 = g_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), M, M)
          EP_avg =&
                   AVG_Y(EP_s(IJK2, M), EP_s(IJK2N, M), J_OF(IJK2))
 
          TH_avg = AVG_Y(&
              AVG_Z(Theta_m(IJK2, M), Theta_m(IJK1, M), K_OF(IJK2)),&
              AVG_Z(Theta_m(IJPK2, M), Theta_m(IJPKP2, M), K_OF(IJPK2)),&
              J_OF(IJK2))
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
	  VSLIP = DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
 
        ELSEIF(FCELL .EQ. 'E')THEN
           IJK2N= NORTH_OF(IJK2)
           IJPK2= JP_OF(IJK2)
           IMJPK2= IM_OF(IJPK2)
          g0 = g_0AVG(IJK2, IJK2N, 'Y',J_OF(IJK2), M, M)
          EP_avg =&
                AVG_Y(EP_s(IJK2, M), EP_s(IJK2N, M),J_OF(IJK2))
 
          TH_avg = AVG_Y(&
                AVG_X(Theta_m(IJK1,M),Theta_m(IJK2, M),I_OF(IJK1)),&
                AVG_X(Theta_m(IMJPK2,M),Theta_m(IJPK2, M),I_OF(IMJPK2)),&
                J_OF(IJK1))
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
	  VSLIP = DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
 
        ELSEIF(FCELL .EQ. 'W')THEN
           IJK2N= NORTH_OF(IJK2)
           IJPK2= JP_OF(IJK2)
           IPJPK2= IP_OF(IJPK2)
 
          g0 = g_0AVG(IJK2, IJK2N, 'Y',J_OF(IJK2), M, M)
          EP_avg =&
                AVG_Y(EP_s(IJK2, M), EP_s(IJK2N, M),J_OF(IJK2))
 
          TH_avg = AVG_Y(&
                AVG_X(Theta_m(IJK2,M),Theta_m(IJK1, M),I_OF(IJK2)),&
                AVG_X(Theta_m(IJPK2,M),Theta_m(IPJPK2, M),I_OF(IJPK2)),&
                J_OF(IJK2))
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
	  VSLIP = DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
        ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
         CALL WRITE_ERROR('CALC_GRBDRY', LINE, 1)
        ENDIF
 
      ELSEIF (COM .EQ. 'W')THEN
        IF(FCELL .EQ. 'N')THEN
          IJK2T = TOP_OF(IJK2)
          IJKP2 = KP_OF(IJK2)
          IJMKP2 = JM_OF(IJKP2)
          g0 = g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, M)
          EP_avg =&
                   AVG_Z(EP_s(IJK2, M), EP_s(IJK2T, M), K_OF(IJK2))
 
          TH_avg = AVG_Y(&
              AVG_Z(Theta_m(IJK1, M), Theta_m(IJMKP2, M), K_OF(IJK1)),&
              AVG_Z(Theta_m(IJK2, M), Theta_m(IJKP2, M), K_OF(IJK2)),&
              J_OF(IJK1))
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
	  VSLIP = DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
 
        ELSEIF(FCELL .EQ. 'S')THEN
          IJK2T = TOP_OF(IJK2)
          IJKP2 = KP_OF(IJK2)
          IJPKP2 = JP_OF(IJKP2)
          g0 = g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, M)
          EP_avg =&
                   AVG_Z(EP_s(IJK2, M), EP_s(IJK2T, M), K_OF(IJK2))
 
          TH_avg = AVG_Y(&
              AVG_Z(Theta_m(IJK2, M), Theta_m(IJKP2, M), K_OF(IJK2)),&
              AVG_Z(Theta_m(IJK1, M), Theta_m(IJPKP2, M), K_OF(IJK1)),&
              J_OF(IJK2))
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
	  VSLIP = DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
 
 
        ELSEIF(FCELL .EQ. 'E')THEN
           IJK2T= TOP_OF(IJK2)
           IJKP2= KP_OF(IJK2)
           IMJKP2 = IM_OF(IJKP2)
          g0 = g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, M)
          EP_avg =&
                AVG_Z(EP_s(IJK2, M), EP_s(IJK2T, M), K_OF(IJK2))
 
          TH_avg = AVG_Z(&
                AVG_X(Theta_m(IJK1,M),Theta_m(IJK2, M),I_OF(IJK1)),&
                AVG_X(Theta_m(IMJKP2,M),Theta_m(IJKP2, M),I_OF(IMJKP2)),&
                K_OF(IJK1))
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
	  VSLIP = DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
 
        ELSEIF(FCELL .EQ. 'W')THEN
           IJK2T= TOP_OF(IJK2)
           IJKP2= KP_OF(IJK2)
           IPJKP2= IP_OF(IJKP2)
          g0 = g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, M)
          EP_avg =&
                AVG_Z(EP_s(IJK2, M), EP_s(IJK2T, M), K_OF(IJK2))
 
          TH_avg = AVG_Z(&
                AVG_X(Theta_m(IJK2,M),Theta_m(IJK1, M),I_OF(IJK2)),&
                AVG_X(Theta_m(IJKP2,M),Theta_m(IPJKP2, M),I_OF(IJKP2)),&
                K_OF(IJK2))
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
	  VSLIP = DSQRT( USCM**2 + VSCM**2 + WSCM**2 )
 
          Hw = F_Hw(g0, EP_avg, TH_avg, Mu_g_avg, RO_g_avg,   &
	            K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
        ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
         CALL WRITE_ERROR('CALC_GRBDRY', LINE, 1)
        ENDIF
 
      ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown COM'
         CALL WRITE_ERROR('CALC_GRBDRY', LINE, 1)
      ENDIF
!
      RETURN
      END
 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: F_HW(g0, EPS, TH, Mu_g_avg, RO_g_avg, VREL,M)          C
!  Purpose: Function for hw                                            C
!                                                                      C
!  Author: K. Agrawal & A. Srivastava, Princeton Univ. Date: 24-JAN-98 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: EPS, TH, C_e, RO_s, D_p(M)                    C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables: F_2, Mu_s, Mu, Mu_b, Eta, Mu_g_avg, RO_g_avg,      C
!                   VREL, C_d, Beta                                    C
!                                                                      C
!  Modified: Sofiane Benyahia, Fluent Inc.             Date: 02-FEB-05 C
!  Purpose: Include conductivity defined by Simonin and Ahmadi         C
!           Also included Jenkins small frictional limit               C
!                                                                      C
!  Literature/Document References: See calcmu_s.f for ref. on Simonin  C
!  and Ahmadi models; for Jenkins BC: Jenkins and Louge, Phys. fluids  C
!  9 (10), 2835. See equation (2) in the paper                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION F_HW(g0,EPS,TH,Mu_g_avg,RO_g_avg,&
                                     K_12_avg, Tau_12_avg, Tau_1_avg, &
				     VREL, VSLIP, M)
 
      USE param 
      USE param1 
      USE constant
      USE physprop
      USE run
      USE fldvar
      IMPLICIT NONE
 
      INTEGER          M
!
!              Average solids volume fraction
      DOUBLE PRECISION EPS
 
!              Average theta_m
      DOUBLE PRECISION Th
 
!              Coefficient of 2nd term
      DOUBLE PRECISION F_2
 
!              Coefficient of 1st term
      DOUBLE PRECISION Mu_s
 
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
 
!              Average gas viscosity
      DOUBLE PRECISION Mu_g_avg
 
!              Reynolds number based on slip velocity
      DOUBLE PRECISION Re_g
 
!              Friction Factor in drag coefficient
      DOUBLE PRECISION C_d
 
!              Drag Coefficient
      DOUBLE PRECISION Beta
 
!              Viscosity corrected for interstitial fluid effects
      DOUBLE PRECISION Mu_star
 
!              Radial distribution function 
      DOUBLE PRECISION g0, G_0
 
!              Constants in Simonin model
      DOUBLE PRECISION Sigma_c, Tau_2_c, Tau_12_st, Nu_t
      DOUBLE PRECISION Tau_2, zeta_c_2, MU_2_T_Kin, Mu_2_Col
 
!     In F_2 and Mu a DSQRT(T) has been left out as it appears in both
!     terms and thus cancels out upon dividing the former by the latter
!
!     The above statement was not implemented because Simonin viscosity
!     doesn't have a sqrt(th) directly available to use this simplification.
!     Sof --> 01/31/05
!
      G_0 = g0
      
! modify F_2 if Jenkins BC is used (sof)    
 
      IF(JENKINS) THEN
!
        IF (VSLIP == ZERO) THEN
! if solids velocity field is initialized to zero, use free slip bc
	  F_2 = zero
!
	ELSE IF(AHMADI) THEN
! Ahmadi model uses different solids pressure model
!
! the coefficient mu in Jenkins paper is defined as tan_Phi_w, that's how
! I understand it from soil mechanic papers, i.e., G.I. Tardos, powder
! Tech. 92 (1997), 61-74. See his equation (1). Define Phi_w in mfix.dat!
!
          F_2 = tan_Phi_w*RO_s(M)*EPS* &
	        ((ONE + 4.0*EPS*G_0) + HALF*(ONE -C_e*C_e))*TH/VSLIP
!
! here F_2 divided by VSLIP to use the same bc as Johnson&Jackson
!
        ELSE  ! Simonin or granular models use same solids pressure
          F_2 = tan_Phi_w*RO_s(M)*EPS*(1d0+ 4. * Eta *EPs*G_0)*TH/VSLIP
	ENDIF !VSLIP == ZERO
!
      ELSE ! no change to the original code if Jenkins BC not used
 
        F_2 = (PHIP*DSQRT(3d0*TH)*Pi*RO_s(M)*EPS*G_0)&
              /(6d0*EPS_max)
!
      ENDIF !for Jenkins
 
      Mu = (5d0*DSQRT(Pi*TH)*D_p(M)*RO_s(M))/96d0
 
      Mu_b = (256d0*Mu*EPS*EPS*G_0)/(5d0*Pi)
 
      Re_g = (1d0-EPS)*RO_g_avg*D_p(M)*VREL/Mu_g_avg
      IF (Re_g.lt.1000d0) THEN
         C_d = (24./(Re_g+SMALL_NUMBER))*(1d0 + 0.15 * Re_g**0.687)
      ELSE
         C_d = 0.44d0
      ENDIF
      Beta = 0.75d0*C_d*Ro_g_avg*(1-EPS)*EPS*VREL&
                *((1-EPS)**(-2.65d0))/D_p(M)
! particle relaxation time
      Tau_12_st = EPS*RO_s(M)/(Beta+small_number)
      
      Beta = SWITCH*Beta
 
!     SWITCH enables us to turn on/off the modification to the
!     particulate phase viscosity. If we want to simulate gas-particle
!     flow then SWITCH=1 to incorporate the effect of drag on the
!     particle viscosity. If we want to simulate granular flow
!     without the effects of an interstitial gas, SWITCH=0.
 
      IF(Beta .LT. SMALL_NUMBER)THEN
        Mu_star = Mu
		
      ELSEIF(TH .LT. SMALL_NUMBER)THEN
        MU_star = ZERO
	
      ELSE
        Mu_star = Mu/(1+(2d0*Beta*Mu/(RO_s(M)*RO_s(M)*EPS*EPS*G_0*TH)))
	
      ENDIF
 
      Mu_s = ((2d0+ALPHA)/3d0)*((Mu_star/(Eta*(2d0-Eta)*G_0))*&
             (1d0+1.6d0*Eta*EPS*G_0)*&
             (1d0+1.6d0*Eta*(3d0*Eta-2d0)*EPS*G_0) +&
             (0.6d0*Mu_b*Eta))
 
      IF(SIMONIN) THEN !see calc_mu_s for explanation of these definitions
!
        Sigma_c = (ONE+ C_e)*(3.d0-C_e)/5.d0
        Tau_2_c = D_P(M)/(6.d0*EPS*G_0*DSQRT(16.d0*(TH+Small_number)/PI))
	zeta_c_2= 2./5.*(ONE+ C_e)*(3.d0*C_e-ONE)
	Nu_t =  Tau_12_avg/Tau_12_st
        Tau_2 = ONE/(2./Tau_12_st+Sigma_c/Tau_2_c)
!
	MU_2_T_Kin = (2.0/3.0*K_12_avg*Nu_t + TH * &
                     (ONE+ zeta_c_2*EPS*G_0))*Tau_2
!
	Mu_2_Col = 8./5.*EPS*G_0*Eta* (MU_2_T_Kin+ &
                   D_p(M)*DSQRT(TH/PI))
!
	Mu_s = EPS*RO_s(M)*(MU_2_T_Kin + Mu_2_Col)
!
      ELSE IF(AHMADI) THEN
!
	Mu_s = ONE/(ONE+ Tau_1_avg/Tau_12_st * (ONE-EPS/EPS_max)**3)&
	       *0.1045*(ONE/G_0+3.2*EPS+12.1824*G_0*EPS*EPS)  &
	       *D_p(M)*RO_s(M)* DSQRT(TH)
      ENDIF
        
 
      F_HW =  F_2/Mu_s
 
      RETURN
      END
      
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
