!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BC_THETA( M, A_m, B_m, IER)                            C
!                                                                      C
!  Purpose: Implementation of Johnson & Jackson boundary conditions    C
!  for the pseudo-thermal temperature.                                 C
!                                                                      C
!                                                                      C
!  Author: Kapil Agrawal, Princeton University        Date: 14-MAR-98  C
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
      SUBROUTINE BC_THETA(M, A_m, B_m, IER)
!
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE toleranc 
      USE run
      USE physprop
      USE fldvar
      USE visc_s
      USE geometry
      USE output
      USE indices
      USE bc
      USE compar         !//d
      use mpi_utility    !//d
      IMPLICIT NONE
!
!  Function subroutines
!
!
!  Local variables
!
!
!                      Error index
      INTEGER          IER
!
!                      Boundary condition
      INTEGER          L
!
!                      Indices
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK,&
                       IM, JM, KM
!
!                      Solids phase
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
 
!
!                      Wall momentum or granular energy coefficient
      DOUBLE PRECISION Gw, Hw, Cw
!
!  Function subroutines
!
!
!  Statement functions
!
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
 
!
!  Setup Johnson and Jackson Pseudo-thermal temp B.C.
!
      DO 500 L = 1, DIMENSION_BC
        IF( BC_DEFINED(L) ) THEN
          IF(BC_TYPE(L) .EQ. 'NO_SLIP_WALL' .OR.&
             BC_TYPE(L) .EQ. 'FREE_SLIP_WALL' .OR.&
             BC_TYPE(L) .EQ. 'PAR_SLIP_WALL' ) THEN
            I1 = BC_I_w(L)
            I2 = BC_I_e(L)
            J1 = BC_J_s(L)
            J2 = BC_J_n(L)
            K1 = BC_K_b(L)
            K2 = BC_K_t(L)
            IF(BC_JJ_PS(L).GT.0) THEN
              DO 120 K = K1, K2
              DO 110 J = J1, J2
              DO 100 I = I1, I2
!// 360 0105 Check if current i,j,k resides on this PE	    
	       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE	      
                IJK   = FUNIJK(I, J, K)
                IM    = Im1(I)
                JM    = Jm1(J)
                KM    = Km1(K)
                A_m(IJK, e, M) =  ZERO
                A_m(IJK, w, M) =  ZERO
                A_m(IJK, n, M) =  ZERO
                A_m(IJK, s, M) =  ZERO
                A_m(IJK, t, M) =  ZERO
                A_m(IJK, b, M) =  ZERO
                A_m(IJK, 0, M) = -ONE
                b_m(IJK, M)    =  ZERO
                IF(FLUID_AT(EAST_OF(IJK)))THEN
                   IF(EP_s(EAST_OF(IJK),M).LE.DIL_EP_s)THEN
                     A_m(IJK, e, M) = -ONE
                   ELSE
                      IF (BC_JJ_PS(L).EQ.3) THEN
                       Gw = 1d0
                       Hw = 0d0
                       Cw = 0d0
 
                      ELSE
                       CALL CALC_THETA_BC(IJK,EAST_OF(IJK),'E',M,L,gw,&
                                        hw,cw)
                      ENDIF
 
                      A_m(IJK, e, M) = -(HALF*hw - oDX_E(I)*gw)
                      A_m(IJK, 0, M) = -(HALF*hw + oDX_E(I)*gw)
 
                      IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M)    = -cw
                      ELSE
                        b_m(IJK,M) = ZERO
                      ENDIF
                   ENDIF
 
                ELSEIF(FLUID_AT(WEST_OF(IJK)))THEN
                   IF(EP_s(WEST_OF(IJK),M).LE.DIL_EP_s)THEN
                     A_m(IJK, w, M) = -ONE
                   ELSE
                     IF (BC_JJ_PS(L).EQ.3) THEN
                      Gw = 1d0
                      Hw = 0d0
                      Cw = 0d0
                     ELSE
                      CALL CALC_THETA_BC(IJK,WEST_OF(IJK),'W',M,L,gw,&
                                        hw,cw)
                     ENDIF
 
                     A_m(IJK, w, M) = -(HALF*hw - oDX_E(IM)*gw)
                     A_m(IJK, 0, M) = -(HALF*hw + oDX_E(IM)*gw)
 
                     IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M)    = -cw
                     ELSE
                        b_m(IJK,M) = ZERO
                     ENDIF
                   ENDIF
 
                ELSEIF(FLUID_AT(NORTH_OF(IJK)))THEN
                   IF(EP_s(NORTH_OF(IJK),M).LE.DIL_EP_s)THEN
                     A_m(IJK, n, M) = -ONE
                   ELSE
                     IF (BC_JJ_PS(L).EQ.3) THEN
                      Gw = 1d0
                      Hw = 0d0
                      Cw = 0d0
                     ELSE
 
                      CALL CALC_THETA_BC(IJK,NORTH_OF(IJK),'N',M,L,gw,&
                                        hw,cw)
                     ENDIF
 
                     A_m(IJK, n, M) = -(HALF*hw - oDY_N(J)*gw)
                     A_m(IJK, 0, M) = -(HALF*hw + oDY_N(J)*gw)
 
                     IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M)    = -cw
                     ELSE
                        b_m(IJK,M) = ZERO
                     ENDIF
                   ENDIF
 
                ELSEIF(FLUID_AT(SOUTH_OF(IJK)))THEN
                   IF(EP_s(SOUTH_OF(IJK),M).LE.DIL_EP_s)THEN
                     A_m(IJK, s, M) = -ONE
                   ELSE
                     IF (BC_JJ_PS(L).EQ.3) THEN
                      Gw = 1d0
                      Hw = 0d0
                      Cw = 0d0
                     ELSE
                      CALL CALC_THETA_BC(IJK,SOUTH_OF(IJK),'S',M,L,gw,&
                                        hw,cw)
 
                     ENDIF
 
                     A_m(IJK, s, M) = -(HALF*hw - oDY_N(JM)*gw)
                     A_m(IJK, 0, M) = -(HALF*hw + oDY_N(JM)*gw)
 
                     IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M)    = -cw
                     ELSE
                        b_m(IJK,M) = ZERO
                     ENDIF
                   ENDIF
 
                ELSEIF(FLUID_AT(TOP_OF(IJK)))THEN
                   IF(EP_s(TOP_OF(IJK),M).LE.DIL_EP_s)THEN
                     A_m(IJK, t, M) = -ONE
                   ELSE
                     IF (BC_JJ_PS(L).EQ.3) THEN
                      Gw = 1d0
                      Hw = 0d0
                      Cw = 0d0
                     ELSE
                      CALL CALC_THETA_BC(IJK,TOP_OF(IJK),'T',M,L,gw,&
                                        hw,cw)
                     ENDIF
 
                     A_m(IJK, t, M) = -(HALF*hw - oX(I)*oDZ_T(K)*gw)
                     A_m(IJK, 0, M) = -(HALF*hw + oX(I)*oDZ_T(K)*gw)
 
                     IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M)    = -cw
                     ELSE
                        b_m(IJK,M) = ZERO
                     ENDIF
                   ENDIF
 
                ELSEIF(FLUID_AT(BOTTOM_OF(IJK)))THEN
                   IF(EP_s(BOTTOM_OF(IJK),M).LE.DIL_EP_s)THEN
                     A_m(IJK, b , M) = -ONE
                   ELSE
                     IF (BC_JJ_PS(L).EQ.3) THEN
                      Gw = 1d0
                      Hw = 0d0
                      Cw = 0d0
                     ELSE
                      CALL CALC_THETA_BC(IJK,BOTTOM_OF(IJK),'B',M,L,&
      	                                gw,hw,cw)
                     ENDIF
 
                     A_m(IJK, b, M) = -(HALF*hw - oX(I)*oDZ_T(KM)*gw)
                     A_m(IJK, 0, M) = -(HALF*hw + oX(I)*oDZ_T(KM)*gw)
 
                     IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M)    = -cw
                     ELSE
                        b_m(IJK,M) = ZERO
                     ENDIF
                   ENDIF
                ENDIF
 100          CONTINUE
 110          CONTINUE
 120          CONTINUE
            ENDIF
          ENDIF
        ENDIF
500   CONTINUE
!
      RETURN
      END
 
 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_THETA_BC(IJK1,IJK2,FCELL,M,L,Gw,Hw,Cw)            C
!                                                                      C
!  Purpose: Implementation of Johnson & Jackson boundary conditions    C
!  for the pseudo-thermal temperature.                                 C
!                                                                      C
!                                                                      C
!  Author: Kapil Agrawal, Princeton University        Date: 14-MAR-98  C
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
      SUBROUTINE CALC_THETA_BC(IJK1,IJK2,FCELL,M,L,Gw,Hw,Cw)
!
!
!  Include param.inc file to specify parameter values
!
      USE param 
      USE param1 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE geometry
      USE indices
      USE bc
      USE compar         !//d
      use mpi_utility    !//d
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
      INTEGER          IJK2E, IJK2SE, IJMK2E, IJK2NE, IJPK2E, IJK2BE
      INTEGER          IJKM2E, IJK2N, IJK2NB, IJKM2N, IJK2TE, IJKP2E
      INTEGER          IJK2NT, IJKP2N, IJK2NW, IMJK2N, IPJK2N, IJK2T
      INTEGER          IJK2TS, IJMK2T, IJK2TN, IJPK2T, IJK2TW, IMJK2T
      INTEGER          IPJK2T
!
!                      Average scalars
      DOUBLE PRECISION EP_avg, TH_avg, Mu_g_avg, RO_g_avg
!
!                      The location (e,w,n...) of fluid cell
      CHARACTER        FCELL
!
!                      Solids phase index
      INTEGER          M
!
!                      Wall momentum or granular energy coefficient
      DOUBLE PRECISION Gw, Hw, Cw
!
!                      values of U_sm, V_sm, W_sm at appropriate place
!                      on boundary wall
      DOUBLE PRECISION USCM, VSCM, WSCM
!
!                      values of U_g, V_g, W_g at appropriate place
!                      on boundary wall
      DOUBLE PRECISION UGC, VGC, WGC
!
!                      Magnitude of gas-solids relative velocity
      DOUBLE PRECISION VREL
 
!              Square of slip velocity between wall and particles
      DOUBLE PRECISION VSLIPSQ
 
!
!                      Error message
      CHARACTER*80     LINE(1)
!
!                      Index corresponding to boundary condition
      INTEGER          L
 
 
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
 
      IF(FCELL .EQ. 'N')THEN
 
        EP_avg = EP_s(IJK2,M)
        TH_avg = AVG_Y(Theta_m(IJK1, M),Theta_m(IJK2, M),J_OF(IJK1))
        Mu_g_avg = Mu_g(IJK2)
        RO_g_avg = RO_g(IJK2)
 
!     Calculate velocity components at i, j+1/2, k (relative to IJK1)
        UGC  = AVG_Y(AVG_X_E(U_g(IM_OF(IJK1)),U_g(IJK1),I_OF(IJK1)),&
                     AVG_X_E(U_g(IM_OF(IJK2)),U_g(IJK2),I_OF(IJK2)),&
                     J_OF(IJK1))
        VGC  = V_g(IJK1)
        WGC  = AVG_Y(AVG_Z_T(W_g(KM_OF(IJK1)), W_g(IJK1)),&
                     AVG_Z_T(W_g(KM_OF(IJK2)), W_g(IJK2)),&
                     J_OF(IJK1))
        USCM = AVG_Y(AVG_X_E(U_s(IM_OF(IJK1),M),U_s(IJK1,M),I_OF(IJK1)),&
                     AVG_X_E(U_s(IM_OF(IJK2),M),U_s(IJK2,M),I_OF(IJK2)),&
                     J_OF(IJK1))
        VSCM = V_s(IJK1,M)
        WSCM = AVG_Y(AVG_Z_T(W_s(KM_OF(IJK1),M), W_s(IJK1,M)),&
                     AVG_Z_T(W_s(KM_OF(IJK2),M), W_s(IJK2,M)),&
                     J_OF(IJK1))
!
!     magnitude of gas-solids relative velocity
!
        VREL =&
          DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )
 
        VSLIPSQ=(WSCM-BC_Ww_s(L, M))**2 + (USCM-BC_Uw_s(L, M))**2
 
 
      ELSEIF(FCELL .EQ. 'S')THEN
 
        EP_avg = EP_s(IJK2,M)
        TH_avg = AVG_Y(Theta_m(IJK2, M),Theta_m(IJK1, M),J_OF(IJK2))
        Mu_g_avg = Mu_g(IJK2)
        RO_g_avg = RO_g(IJK2)
 
!     Calculate velocity components at i, j+1/2, k (relative to IJK2)
        UGC  = AVG_Y(AVG_X_E(U_g(IM_OF(IJK2)),U_g(IJK2),I_OF(IJK2)),&
                     AVG_X_E(U_g(IM_OF(IJK1)),U_g(IJK1),I_OF(IJK1)),&
                     J_OF(IJK2))
        VGC  = V_g(IJK2)
        WGC  = AVG_Y(AVG_Z_T(W_g(KM_OF(IJK2)), W_g(IJK2)),&
                     AVG_Z_T(W_g(KM_OF(IJK1)), W_g(IJK1)),&
                     J_OF(IJK2))
        USCM = AVG_Y(AVG_X_E(U_s(IM_OF(IJK2),M),U_s(IJK2,M),I_OF(IJK2)),&
                     AVG_X_E(U_s(IM_OF(IJK1),M),U_s(IJK1,M),I_OF(IJK1)),&
                     J_OF(IJK2))
        VSCM = V_s(IJK2,M)
        WSCM = AVG_Y(AVG_Z_T(W_s(KM_OF(IJK2),M), W_s(IJK2,M)),&
                     AVG_Z_T(W_s(KM_OF(IJK1),M), W_s(IJK1,M)),&
                     J_OF(IJK2))
!
!     magnitude of gas-solids relative velocity
!
        VREL =&
          DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )
 
        VSLIPSQ=(WSCM-BC_Ww_s(L, M))**2 + (USCM-BC_Uw_s(L, M))**2
 
      ELSEIF(FCELL== 'E')THEN
 
        EP_avg = EP_s(IJK2,M)
        TH_avg = AVG_X(Theta_m(IJK1, M),Theta_m(IJK2, M),I_OF(IJK1))
        Mu_g_avg = Mu_g(IJK2)
        RO_g_avg = RO_g(IJK2)
 
!     Calculate velocity components at i+1/2, j, k (relative to IJK1)
        UGC  = U_g(IJK1)
        VGC  = AVG_X(AVG_Y_N(V_g(JM_OF(IJK1)), V_g(IJK1)),&
                     AVG_Y_N(V_g(JM_OF(IJK2)), V_g(IJK2)),&
                     I_OF(IJK1))
        WGC  = AVG_X(AVG_Z_T(W_g(KM_OF(IJK1)), W_g(IJK1)),&
                     AVG_Z_T(W_g(KM_OF(IJK2)), W_g(IJK2)),&
                     I_OF(IJK1))
        USCM = U_s(IJK1,M)
        VSCM = AVG_X(AVG_Y_N(V_s(JM_OF(IJK1),M), V_s(IJK1,M)),&
                     AVG_Y_N(V_s(JM_OF(IJK2),M), V_s(IJK2,M)),&
                     I_OF(IJK1))
        WSCM = AVG_X(AVG_Z_T(W_s(KM_OF(IJK1),M), W_s(IJK1,M)),&
                     AVG_Z_T(W_s(KM_OF(IJK2),M), W_s(IJK2,M)),&
                     I_OF(IJK1))
!
!     magnitude of gas-solids relative velocity
!
        VREL =&
          DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )
 
        VSLIPSQ=(WSCM-BC_Ww_s(L, M))**2 + (VSCM-BC_Vw_s(L, M))**2
 
 
      ELSEIF(FCELL== 'W')THEN
 
        EP_avg = EP_s(IJK2,M)
        TH_avg = AVG_X(Theta_m(IJK2, M),Theta_m(IJK1, M),I_OF(IJK2))
        Mu_g_avg = Mu_g(IJK2)
        RO_g_avg = RO_g(IJK2)
 
!     Calculate velocity components at i+1/2, j, k (relative to IJK2)
        UGC  = U_g(IJK2)
        VGC  = AVG_X(AVG_Y_N(V_g(JM_OF(IJK2)), V_g(IJK2)),&
                     AVG_Y_N(V_g(JM_OF(IJK1)), V_g(IJK1)),&
                     I_OF(IJK2))
        WGC  = AVG_X(AVG_Z_T(W_g(KM_OF(IJK2)), W_g(IJK2)),&
                     AVG_Z_T(W_g(KM_OF(IJK1)), W_g(IJK1)),&
                     I_OF(IJK2))
        USCM = U_s(IJK2,M)
        VSCM = AVG_X(AVG_Y_N(V_s(JM_OF(IJK2),M), V_s(IJK2,M)),&
                     AVG_Y_N(V_s(JM_OF(IJK1),M), V_s(IJK1,M)),&
                     I_OF(IJK2))
        WSCM = AVG_X(AVG_Z_T(W_s(KM_OF(IJK2),M), W_s(IJK2,M)),&
                     AVG_Z_T(W_s(KM_OF(IJK1),M), W_s(IJK1,M)),&
                     I_OF(IJK2))
!
!     magnitude of gas-solids relative velocity
!
        VREL =&
          DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )
 
        VSLIPSQ=(WSCM-BC_Ww_s(L, M))**2 + (VSCM-BC_Vw_s(L, M))**2
 
 
      ELSEIF(FCELL== 'T')THEN
 
        EP_avg = EP_s(IJK2,M)
        TH_avg = AVG_Z(Theta_m(IJK1, M),Theta_m(IJK2, M),K_OF(IJK1))
        Mu_g_avg = Mu_g(IJK2)
        RO_g_avg = RO_g(IJK2)
 
!     Calculate velocity components at i, j, k+1/2 (relative to IJK1)
        UGC  = AVG_Z(AVG_X_E(U_g(IM_OF(IJK1)),U_g(IJK1),I_OF(IJK1)),&
                     AVG_X_E(U_g(IM_OF(IJK2)),U_g(IJK2),I_OF(IJK2)),&
                     K_OF(IJK1))
        VGC  = AVG_Z(AVG_Y_N(V_g(JM_OF(IJK1)), V_g(IJK1)),&
                     AVG_Y_N(V_g(JM_OF(IJK2)), V_g(IJK2)),&
                     K_OF(IJK1))
        WGC  = W_g(IJK1)
        USCM = AVG_Z(AVG_X_E(U_s(IM_OF(IJK1),M),U_s(IJK1,M),I_OF(IJK1)),&
                     AVG_X_E(U_s(IM_OF(IJK2),M),U_s(IJK2,M),I_OF(IJK2)),&
                     K_OF(IJK1))
        VSCM = AVG_Z(AVG_Y_N(V_s(JM_OF(IJK1),M), V_s(IJK1,M)),&
                     AVG_Y_N(V_s(JM_OF(IJK2),M), V_s(IJK2,M)),&
                     K_OF(IJK1))
        WSCM = W_s(IJK1,M)
!
!     magnitude of gas-solids relative velocity
!
        VREL =&
          DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )
 
 
        VSLIPSQ=(VSCM-BC_Vw_s(L, M))**2 + (USCM-BC_Uw_s(L, M))**2
 
 
      ELSEIF(FCELL== 'B')THEN
 
        EP_avg = EP_s(IJK2,M)
        TH_avg = AVG_Z(Theta_m(IJK2, M),Theta_m(IJK1, M),K_OF(IJK2))
        Mu_g_avg = Mu_g(IJK2)
        RO_g_avg = RO_g(IJK2)
 
!     Calculate velocity components at i, j, k+1/2 (relative to IJK2)
        UGC  = AVG_Z(AVG_X_E(U_g(IM_OF(IJK2)),U_g(IJK2),I_OF(IJK2)),&
                     AVG_X_E(U_g(IM_OF(IJK1)),U_g(IJK1),I_OF(IJK1)),&
                     K_OF(IJK2))
        VGC  = AVG_Z(AVG_Y_N(V_g(JM_OF(IJK2)), V_g(IJK2)),&
                     AVG_Y_N(V_g(JM_OF(IJK1)), V_g(IJK1)),&
                     K_OF(IJK2))
        WGC  = W_g(IJK2)
        USCM = AVG_Z(AVG_X_E(U_s(IM_OF(IJK2),M),U_s(IJK2,M),I_OF(IJK2)),&
                     AVG_X_E(U_s(IM_OF(IJK1),M),U_s(IJK1,M),I_OF(IJK1)),&
                     K_OF(IJK2))
        VSCM = AVG_Z(AVG_Y_N(V_s(JM_OF(IJK2),M), V_s(IJK2,M)),&
                     AVG_Y_N(V_s(JM_OF(IJK1),M), V_s(IJK1,M)),&
                     K_OF(IJK2))
        WSCM = W_s(IJK2,M)
!
!     magnitude of gas-solids relative velocity
!
        VREL =&
          DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )
 
        VSLIPSQ=(VSCM-BC_Vw_s(L, M))**2 + (USCM-BC_Uw_s(L, M))**2
 
 
      ELSE
        WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
        CALL WRITE_ERROR('CALC_THETA_BC', LINE, 1)
	call exitMPI(myPE)            !//d
      ENDIF
 
      CALL THETA_Hw_Cw(EP_avg,TH_avg,Mu_g_avg,RO_g_avg,VREL,VSLIPSQ,M,&
                       Gw,Hw,Cw,L)
!
      RETURN
      END
 
 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SUBROUTINE THETA_HW_CW(EPS, TH, Mu_g_avg, RO_g_avg,    C
!                                      VREL,VSLIPSQ,M,GW,HW,CW,L)      C
!  Purpose: Subroutine for hw and cw                                   C
!                                                                      C
!  Author: Kapil Agrawal, Princeton University         Date: 15-MAR-98 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: EPS, TH, C_e, RO_s, D_p(M), Eta, e_w          C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables: F_2, Mu_s, Mu, Mu_b, Mu_g_avg, RO_g_avg,           C
!                   VREL, C_d, Beta                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE THETA_HW_CW(EPS,TH,Mu_g_avg,RO_g_avg,VREL,VSLIPSQ,M,&
                             GW,HW,CW,L)
 

      USE param 
      USE param1 
      USE physprop
      USE constant
      USE fldvar
      USE toleranc 
      USE bc
      USE compar      !//d
      IMPLICIT NONE
 
!
!                      Solids phase index
      INTEGER          M
 
!                      The location (e,w,n...) of fluid cell
      CHARACTER        FCELL
!
!              Average solids volume fraction
      DOUBLE PRECISION EPS
 
!              Average theta_m
      DOUBLE PRECISION Th
 
!              Coefficient of 1st term
      DOUBLE PRECISION K_1
!      DOUBLE PRECISION K_10
 
!              Viscosity
      DOUBLE PRECISION Lambda
!      DOUBLE PRECISION Lambda0
 
!              corrected for interstitial fluid effects
      DOUBLE PRECISION Lambda_star
 
!              Magnitude of slip velocity between two phases
      DOUBLE PRECISION VREL
 
!              Square of slip velocity between wall and particles
      DOUBLE PRECISION VSLIPSQ
 
!              Average gas density
      DOUBLE PRECISION RO_g_avg
 
!              Average gas viscosity
      DOUBLE PRECISION Mu_g_avg
 
!              Reynolds number based on slip velocity
      DOUBLE PRECISION Re_g
 
!              Friction Factor in drag coefficient
      DOUBLE PRECISION C_d
 
!              Drag Coefficient
      DOUBLE PRECISION Beta
      DOUBLE PRECISION Beta0
 
!              Coefficients in boundary conditions
      DOUBLE PRECISION GW, HW, CW
 
 
!              Radial distribution function (Carnahan & Starling)
      DOUBLE PRECISION G_0, G_0EP
 
!
!                      Error message
      CHARACTER*80     LINE
!
!                      Index corresponding to boundary condition
      INTEGER          L
 
      IF(TH .LE. ZERO)THEN
        TH = 1e-8
        if (myPE.eq.PE_IO) then   !//??????  pnicol : on any PE ????
	   WRITE(*,*)'Warning: Negative granular temp at wall set to 1e-8'
!          CALL WRITE_ERROR('THETA_HW_CW', LINE, 1)
        end if
      ENDIF
 
!     In F_2 and Mu a DSQRT(T) has been left out as it appears in both
!     terms and thus cancels out upon dividing the former by the latter
 
      G_0 = G_0EP(EPS)
 
      Lambda = 75*RO_s(M)*D_p(M)*DSQRT(Pi*TH)/(48*Eta*(41d0-33d0*Eta))
 
      Re_g = (1d0-EPS)*RO_g_avg*D_p(M)*VREL/Mu_g_avg
      IF (Re_g.lt.1000d0) THEN
         C_d = (24./(Re_g+SMALL_NUMBER))*(1d0 + 0.15 * Re_g**0.687)
      ELSE
         C_d = 0.44d0
      ENDIF
      Beta = SWITCH*0.75d0*C_d*Ro_g_avg*(1-EPS)*EPS*VREL&
                *((1-EPS)**(-2.65d0))/D_p(M)
 
!     SWITCH enables us to turn on/off the modification to the
!     particulate phase viscosity. If we want to simulate gas-particle
!     flow then SWITCH=1 to incorporate the effect of drag on the
!     particle viscosity. If we want to simulate granular flow
!     without the effects of an interstitial gas, SWITCH=0.
 
 
      IF(Beta .LT. SMALL_NUMBER)THEN
        Lambda_star = Lambda
		
      ELSEIF(TH .LT. SMALL_NUMBER)THEN
        Lambda_star = ZERO
	
      ELSE
        Lambda_star = Lambda/(1+(6d0*Beta*Lambda/(5d0*RO_s(M)*RO_s(M)*&
                    EPS*EPS*G_0*TH)))
 
      ENDIF
 
      K_1 = Lambda_star*(((1d0/G_0) + (12d0/5.)*Eta*EPS)&
            *(1d0 + (12d0/5.)*Eta*Eta*(4d0*Eta-3d0)*EPS*G_0)&
           + (64d0/(25d0*Pi))*(41d0-33d0*Eta)*((Eta*EPS)**2)*G_0)
 
      GW = K_1
 
      HW = (Pi*DSQRT(3d0)/(4.*EPS_max))*(1d0-e_w*e_w)*RO_s(M)*EPS*G_0*&
         DSQRT(TH)
 
      CW = (Pi*DSQRT(3d0)/(6.*EPS_max))*PHIP*RO_s(M)*EPS*G_0*DSQRT(TH)&
         *VSLIPSQ
 
      IF (BC_JJ_PS(L).EQ.2) CW=0d0
 
 
      RETURN
      END
