!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_U_FRICTION
!                                                                      C
!  Purpose: Calculate hw and cw for kinetic theory and frictional      C
!           boundary conditions                                        C
!                                                                      C
!                                                                      C
!  Author: Anuj Srivastava, Princeton University      Date: 12-APR-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Notes:
!   This routine is no longer needed and will likely be removed.
!   The call to the pertinent subroutines were moved into calc_grbdry.
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_U_FRICTION(IJK1, IJK2, FCELL, COM, L, M, Gw,&
                                 Hw, Cw)

!-----------------------------------------------
! Modules
!-----------------------------------------------                                 
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
      USE toleranc
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! IJK indices for wall cell (ijk1) and fluid cell (ijk2)
      INTEGER, INTENT(IN) :: IJK1, IJK2
! The location (e,w,n...) of fluid cell
      CHARACTER, INTENT(IN) :: FCELL
! Velocity component (U, V, W)
      CHARACTER :: COM
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Index corresponding to boundary condition
      INTEGER, INTENT(IN) ::  L      
! Wall momentum coefficients:
! 1st term on LHS
      DOUBLE PRECISION, INTENT(INOUT) :: Gw
! 2nd term on LHS
      DOUBLE PRECISION, INTENT(INOUT) :: Hw
! all terms appearing on RHS
      DOUBLE PRECISION, INTENT(INOUT) :: Cw

!-----------------------------------------------
! Local Variables      
!-----------------------------------------------
! IJK indices for fluid cell
      INTEGER :: IJK
! Other indices
      INTEGER :: IJK2E, IPJK2, IPJKM2, IPJKP2, IPJMK2, IPJPK2
      INTEGER :: IJK2N, IJPK2, IJPKM2, IJPKP2, IMJPK2
      INTEGER :: IJK2T, IJKP2, IJMKP2, IMJKP2
! Solids phase index
      INTEGER :: MM
! Average scalars
      DOUBLE PRECISION :: EP_avg, EPg_avg, TH_avg, Mu_g_avg, RO_g_avg,Dp_avg
      DOUBLE PRECISION :: AVGX1, AVGX2, smallTheta
! void fraction at packing
      DOUBLE PRECISION :: ep_star_avg
! Average Simonin and Ahmadi variables (sof)
      DOUBLE PRECISION :: K_12_avg, Tau_12_avg, Tau_1_avg
! Average velocities
      DOUBLE PRECISION :: WGC1, WGC2, WGCM, VGC1, VGC2, UGC1, UGC2
! values of U_sm, V_sm, W_sm at appropriate place on boundary wall
      DOUBLE PRECISION :: USCM, VSCM,WSCM
      DOUBLE PRECISION :: USCM1,USCM2,VSCM1,VSCM2,WSCM1,WSCM2
! values of U_g, V_g, W_g at appropriate place on boundary wall
      DOUBLE PRECISION :: UGC, VGC, WGC
      DOUBLE PRECISION :: USC1,USC2,VSC1,VSC2,WSC1,WSC2
! velocity variables used to standarize dummy argument for different dirs
      DOUBLE PRECISION :: WVELS
      DOUBLE PRECISION :: VELS
! del.u
      DOUBLE PRECISION :: DEL_DOT_U
! S:S
      DOUBLE PRECISION :: S_DDOT_S
! S_dd (d can be x, y or z)
      DOUBLE PRECISION :: S_dd
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION :: VREL
! slip velocity between wall and particles for Jenkins bc (sof)
      DOUBLE PRECISION :: VSLIP
! radial distribution function at contact
      DOUBLE PRECISION :: g0
! Sum of eps*G_0
      DOUBLE PRECISION :: g0EP_avg
! Average Radial distribution function
      DOUBLE PRECISION :: g_0AVG
! Error message
      CHARACTER*80     LINE
!----------------------------------------------- 
! Include statements functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------     
      smallTheta = (to_SI)**4 * ZERO_EP_S
 
 


! Note: EP_s, MU_g, and RO_g are undefined at IJK1 (wall cell).
!       Hence IJK2 (fluid cell) is used in averages.

!-----------------------------------------------      
! Calculations for U momentum equation
      IF (COM .EQ. 'U')THEN

        IJK2E = EAST_OF(IJK2)
        IPJK2 = IP_OF(IJK2)

        EP_avg = AVG_X(EP_s(IJK2, M), EP_s(IJK2E, M), I_OF(IJK2))
        EPg_avg = AVG_X(EP_g(IJK2), EP_g(IJK2E), I_OF(IJK2))
        ep_star_avg = AVG_X(EP_star_array(IJK2), EP_star_array(IJK2E), I_OF(IJK2))
        Mu_g_avg = AVG_X(Mu_g(IJK2), Mu_g(IJK2E), I_OF(IJK2))
        RO_g_avg = AVG_X(RO_g(IJK2), RO_g(IJK2E), I_OF(IJK2))
        DP_avg   = AVG_X(D_P(IJK2,M), D_P(IJK2E,M), I_OF(IJK2))

        IF(.NOT.GRANULAR_ENERGY) THEN
           TH_avg = AVG_X(THETA_M(IJK2,M), THETA_M(IJK2E,M), I_OF(IJK2))
           IF(TH_avg < ZERO) TH_avg = smallTheta
        ENDIF

        IF(SIMONIN .OR. AHMADI) THEN
           K_12_avg = AVG_X(K_12(IJK2), K_12(IJK2E), I_OF(IJK2))
           Tau_12_avg = AVG_X(Tau_12(IJK2), Tau_12(IJK2E), I_OF(IJK2))
           Tau_1_avg = AVG_X(Tau_1(IJK2), Tau_1(IJK2E), I_OF(IJK2))
        ELSE
           K_12_avg = ZERO    
           Tau_12_avg = ZERO
           Tau_1_avg = ZERO
        ENDIF

        g0 = g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), M, M)          
        g0EP_avg = ZERO
        DO MM = 1, SMAX
           g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), M, MM) &
                     * AVG_X(EP_s(IJK2, MM), EP_s(IJK2E, MM), I_OF(IJK2))
        ENDDO

        WVELS = BC_Uw_s(L,M)


        IF(FCELL .EQ. 'N')THEN
          IPJMK2 = JM_OF(IPJK2)
 
! code modified for some corner cells
          IF(GRANULAR_ENERGY) THEN
            AVGX1 = AVG_X(Theta_m(IJK1,M), Theta_m(IPJMK2,M), I_OF(IJK1))
            AVGX2 = AVG_X(Theta_m(IJK2,M), Theta_m(IPJK2,M), I_OF(IJK2))
            IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
            IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
            IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
              TH_avg = smallTheta
            ELSE
              TH_avg = AVG_Y(AVGX1, AVGX2, J_OF(IJK1))
            ENDIF
          ENDIF

! Calculate velocity components at i+1/2, j+1/2, k (relative to IJK1)
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
          VELS = USCM
 
        ELSEIF(FCELL .EQ. 'S')THEN
          IPJPK2= JP_OF(IPJK2)
 
          IF(GRANULAR_ENERGY) THEN
                 AVGX1 = AVG_X(Theta_m(IJK2,M),Theta_m(IPJK2,M),I_OF(IJK2))
             AVGX2 = AVG_X(Theta_m(IJK1,M),Theta_m(IPJPK2,M),I_OF(IJK1))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Y(AVGX1, AVGX2, J_OF(IJK2))
             ENDIF
          ENDIF

! Calculate velocity components at i+1/2, j+1/2, k relative to IJK2
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
          VELS = USCM
 
        ELSEIF(FCELL .EQ. 'T')THEN
          IPJKM2= KM_OF(IPJK2)

          IF(GRANULAR_ENERGY) THEN
                 AVGX1 = AVG_X(Theta_m(IJK1,M),Theta_m(IPJKM2,M),I_OF(IJK1))
             AVGX2 = AVG_X(Theta_m(IJK2,M),Theta_m(IPJK2,M),I_OF(IJK2))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Z(AVGX1, AVGX2, K_OF(IJK1))
             ENDIF
          ENDIF

! Calculate velocity components at i+1/2,j,k-1/2 relative to IJK2
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
          VELS = USCM
 
        ELSEIF(FCELL .EQ. 'B')THEN
          IPJKP2= KP_OF(IPJK2)
 
          IF(GRANULAR_ENERGY) THEN
             AVGX1 = AVG_X(Theta_m(IJK1,M), Theta_m(IPJKP2,M),I_OF(IJK1))
             AVGX2 = AVG_X(Theta_m(IJK2,M), Theta_m(IPJK2,M),I_OF(IJK2))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Z(AVGX1, AVGX2, K_OF(IJK2))
             ENDIF
          ENDIF
 
! Calculate velocity components at i+1/2,j,k-1/2 relative to IJK1
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
          VELS = USCM

        ELSE
           WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
           CALL WRITE_ERROR('CALC_U_FRICTION', LINE, 1)
           CALL exitMPI(myPE)
        ENDIF
 
!-----------------------------------------------      
! Calculations for V momentum equation

      ELSEIF (COM .EQ. 'V')THEN

        IJK2N = NORTH_OF(IJK2)
        IJPK2 = JP_OF(IJK2)


        EP_avg = AVG_Y(EP_s(IJK2, M), EP_s(IJK2N, M), J_OF(IJK2))
        EPg_avg = AVG_Y(EP_g(IJK2), EP_g(IJK2N), J_OF(IJK2))
        ep_star_avg = AVG_Y(EP_star_array(IJK2), EP_star_array(IJK2N), J_OF(IJK2))
        Mu_g_avg = AVG_Y(Mu_g(IJK2), Mu_g(IJK2N), J_OF(IJK2))
        RO_g_avg = AVG_Y(RO_g(IJK2), RO_g(IJK2N), J_OF(IJK2))
        DP_avg   = AVG_Y(D_P(IJK2,M), D_P(IJK2N,M), J_OF(IJK2))

        IF(.NOT.GRANULAR_ENERGY) THEN
           TH_avg = AVG_Y(THETA_M(IJK2,M), THETA_M(IJK2N,M), J_OF(IJK2))
           IF(TH_avg < ZERO) TH_avg = smallTheta
        ENDIF

        IF(SIMONIN .OR. AHMADI) THEN
           K_12_avg = AVG_Y(K_12(IJK2), K_12(IJK2N), J_OF(IJK2))
           Tau_12_avg = AVG_Y(Tau_12(IJK2), Tau_12(IJK2N), J_OF(IJK2))
           Tau_1_avg = AVG_Y(Tau_1(IJK2), Tau_1(IJK2N), J_OF(IJK2))
        ELSE
           K_12_avg = ZERO    
           Tau_12_avg = ZERO
           Tau_1_avg = ZERO
        ENDIF

        g0 = g_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), M, M)          
        g0EP_avg = ZERO
        DO MM = 1, SMAX
           g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), M, MM) &
              * AVG_Y(EP_s(IJK2, MM), EP_s(IJK2N, MM), J_OF(IJK2))
        ENDDO

        WVELS = BC_Vw_s(L,M)


        IF(FCELL .EQ. 'T')THEN
          IJPKM2 = KM_OF(IJPK2)
 
          IF(GRANULAR_ENERGY) THEN
             AVGX1 = AVG_Z(Theta_m(IJK1,M), Theta_m(IJK2,M), K_OF(IJK1))
             AVGX2 = AVG_Z(Theta_m(IJPKM2,M), Theta_m(IJPK2,M), K_OF(IJPKM2))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Y(AVGX1, AVGX2, J_OF(IJK1))
             ENDIF
          ENDIF
 
! Calculate velocity components at i,j+1/2,k+1/2 (relative to IJK1)
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
          VELS = VSCM 

 
        ELSEIF(FCELL .EQ. 'B')THEN
          IJPKP2 = KP_OF(IJPK2)
 
          IF(GRANULAR_ENERGY) THEN
                 AVGX1 = AVG_Z(Theta_m(IJK2,M), Theta_m(IJK1,M), K_OF(IJK2))
             AVGX2 = AVG_Z(Theta_m(IJPK2,M), Theta_m(IJPKP2,M), K_OF(IJPK2))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Y(AVGX1, AVGX2, J_OF(IJK2))
             ENDIF
          ENDIF

! Calculate velocity components at i,j+1/2,k+1/2 (relative to IJK2)
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
          VELS =VSCM 
          
        ELSEIF(FCELL .EQ. 'E')THEN
          IMJPK2= IM_OF(IJPK2)
 
          IF(GRANULAR_ENERGY) THEN
             AVGX1 = AVG_X(Theta_m(IJK1,M),Theta_m(IJK2,M),I_OF(IJK1))
             AVGX2 = AVG_X(Theta_m(IMJPK2,M),Theta_m(IJPK2,M),I_OF(IMJPK2))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Y(AVGX1, AVGX2, J_OF(IJK1))
             ENDIF
          ENDIF
 
! Calculate velocity components at i+1/2,j+1/2,k relative to IJK1
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
          VELS = VSCM


        ELSEIF(FCELL .EQ. 'W')THEN
          IPJPK2= IP_OF(IJPK2)
 
          IF(GRANULAR_ENERGY) THEN
             AVGX1 = AVG_X(Theta_m(IJK2,M),Theta_m(IJK1,M),I_OF(IJK2))
             AVGX2 = AVG_X(Theta_m(IJPK2,M),Theta_m(IPJPK2,M),I_OF(IJPK2))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Y(AVGX1, AVGX2, J_OF(IJK2))
             ENDIF
          ENDIF
 
! Calculate velocity components at i+1/2,j+1/2,k relative to IJK2
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
          VELS = VSCM

        ELSE
           WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
           CALL WRITE_ERROR('CALC_U_FRICTION', LINE, 1)
           CALL exitMPI(myPE)          
        ENDIF
 
!-----------------------------------------------      
! Calculations for W momentum equation

      ELSEIF (COM .EQ. 'W')THEN

        IJK2T = TOP_OF(IJK2)
        IJKP2 = KP_OF(IJK2)

        EP_avg = AVG_Z(EP_s(IJK2, M), EP_s(IJK2T, M), K_OF(IJK2))
        EPg_avg = AVG_Z(EP_g(IJK2), EP_g(IJK2T), K_OF(IJK2))
        ep_star_avg = AVG_Z(EP_star_array(IJK2), EP_star_array(IJK2T), K_OF(IJK2))
        Mu_g_avg = AVG_Z(Mu_g(IJK2), Mu_g(IJK2T), K_OF(IJK2))
        RO_g_avg =  AVG_Z(RO_g(IJK2), RO_g(IJK2T), K_OF(IJK2))
        DP_avg   = AVG_Z(D_P(IJK2,M), D_P(IJK2T,M), K_OF(IJK2))         

        IF(.NOT.GRANULAR_ENERGY) THEN
           TH_avg = AVG_Z(THETA_M(IJK2,M), THETA_M(IJK2T,M), K_OF(IJK2))
           IF(TH_avg < ZERO) TH_avg = smallTheta
        ENDIF

        IF(SIMONIN .OR. AHMADI) THEN
           K_12_avg = AVG_Z(K_12(IJK2), K_12(IJK2T), K_OF(IJK2))
           Tau_12_avg = AVG_Z(Tau_12(IJK2), Tau_12(IJK2T), K_OF(IJK2))
           Tau_1_avg = AVG_Z(Tau_1(IJK2), Tau_1(IJK2T), K_OF(IJK2))
        ELSE
           K_12_avg = ZERO    
           Tau_12_avg = ZERO
           Tau_1_avg = ZERO
        ENDIF

        g0 = g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, M)                 
        g0EP_avg = ZERO
        DO MM = 1, SMAX
           g0EP_avg = g0EP_avg + g_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, MM) &
              * AVG_Z(EP_s(IJK2, MM), EP_s(IJK2T, MM), K_OF(IJK2))
        ENDDO
      
        WVELS = BC_Ww_s(L,M)


        IF(FCELL .EQ. 'N')THEN
          IJMKP2 = JM_OF(IJKP2)
 
          IF(GRANULAR_ENERGY) THEN
             AVGX1 = AVG_Z(Theta_m(IJK1,M), Theta_m(IJMKP2,M), K_OF(IJK1))
             AVGX2 = AVG_Z(Theta_m(IJK2,M), Theta_m(IJKP2,M), K_OF(IJK2))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Y(AVGX1, AVGX2, J_OF(IJK1))
             ENDIF
             ELSE
          ENDIF

! Calculate velocity components at i,j+1/2,k+1/2 (relative to IJK1)
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
          VELS = WSCM 
 
        ELSEIF(FCELL .EQ. 'S')THEN
          IJPKP2 = JP_OF(IJKP2)
 
          IF(GRANULAR_ENERGY) THEN
             AVGX1 = AVG_Z(Theta_m(IJK2,M), Theta_m(IJKP2,M), K_OF(IJK2))
             AVGX2 = AVG_Z(Theta_m(IJK1,M), Theta_m(IJPKP2,M), K_OF(IJK1))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Y(AVGX1, AVGX2, J_OF(IJK2))
             ENDIF
          ENDIF

! Calculate velocity components at i,j+1/2,k+1/2 (relative to IJK2)
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
          VELS = WSCM          
 
        ELSEIF(FCELL .EQ. 'E')THEN
          IMJKP2 = IM_OF(IJKP2)
 
          IF(GRANULAR_ENERGY) THEN
             AVGX1 = AVG_X(Theta_m(IJK1,M),Theta_m(IJK2,M),I_OF(IJK1))
             AVGX2 = AVG_X(Theta_m(IMJKP2,M),Theta_m(IJKP2,M),I_OF(IMJKP2))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Z(AVGX1, AVGX2, K_OF(IJK1))
             ENDIF
          ENDIF
 
! Calculate velocity components at i+1/2,j,k+1/2 relative to IJK1
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
          VELS = WSCM          

        ELSEIF(FCELL .EQ. 'W')THEN
          IPJKP2= IP_OF(IJKP2)
 
          IF(GRANULAR_ENERGY) THEN
             AVGX1 = AVG_X(Theta_m(IJK2,M),Theta_m(IJK1,M),I_OF(IJK2))
             AVGX2 = AVG_X(Theta_m(IJKP2,M),Theta_m(IPJKP2,M),I_OF(IJKP2))
             IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
             IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
             IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
               TH_avg = smallTheta
             ELSE
               TH_avg = AVG_Z(AVGX1, AVGX2, K_OF(IJK2))
             ENDIF
          ENDIF
 
! Calculate velocity components at i+1/2,j,k+1/2 relative to IJK2
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
          VELS = WSCM
 
        ELSE
           WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
           CALL WRITE_ERROR('CALC_U_FRICTION', LINE, 1)
           CALL exitMPI(myPE)         
        ENDIF
 
      ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown COM'
         CALL WRITE_ERROR('CALC_U_FRICTION', LINE, 1)
         CALL exitMPI(myPE)         
      ENDIF   ! end if com=u,v,w

! magnitude of gas-solids relative velocity
      VREL = DSQRT( (UGC-USCM)**2 + (VGC-VSCM)**2 + (WGC-WSCM)**2 )

! slip velocity for use in Jenkins bc (sof)	  
      VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
         + (WSCM-BC_WW_S(L,M))**2 )
 
      IF (FRICTION .AND. (ONE-EP_G(IJK2))>EPS_F_MIN) THEN         
        CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
           S_DDOT_S, S_dd)
 
        CALL CALC_Gw_Hw_Cw(g0, EP_avg, EPg_avg, ep_star_avg, &
           g0EP_avg, TH_avg, Mu_g_avg, RO_g_avg, &
           DP_avg, K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP,&
           DEL_DOT_U, S_DDOT_S, S_dd, VELS, WVELS, M, gw, hw, cw)
      ENDIF



      RETURN
      END SUBROUTINE CALC_U_FRICTION
 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Gw_Hw_Cw
!                                                                      C
!  Purpose: Calculate Gw, Hw, and Cw                                   C
!                                                                      C
!  Author: A. Srivastava & K. Agrawal, Princeton Univ. Date: 10-APR-98 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
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
 
      SUBROUTINE CALC_Gw_Hw_Cw(g0, EPS, EPG, ep_star_avg, &
            g0EP_avg, TH, Mu_g_avg, RO_g_avg, &
            DP_avg, K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP,&
            DEL_U, S_S, S_dd, VEL, W_VEL, M, gw, hw, cw)

!-----------------------------------------------
! Modules
!-----------------------------------------------              
      USE param 
      USE param1 
      USE constant
      USE physprop
      USE fldvar
      USE bc
      USE run
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------      
! Radial distribution function of solids phase M with each
! other solids phase 
      DOUBLE PRECISION, INTENT(IN) :: g0
! Average solids volume fraction of each solids phase
      DOUBLE PRECISION, INTENT(IN) :: EPS
! Average solids and gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPG, ep_star_avg
! Sum of eps*G_0 
      DOUBLE PRECISION, INTENT(IN) :: g0EP_avg 
! Average theta_m
      DOUBLE PRECISION, INTENT(INOUT) :: TH
! Average gas viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mu_g_avg
! Average gas density
      DOUBLE PRECISION, INTENT(IN) :: RO_g_avg
! Average particle diameter of each solids phase
      DOUBLE PRECISION, INTENT(IN) :: DP_avg
! Average cross-correlation K_12 and lagrangian integral time-scale
      DOUBLE PRECISION, INTENT(IN) :: K_12_avg, Tau_12_avg, Tau_1_avg
! Magnitude of slip velocity between two phases
      DOUBLE PRECISION, INTENT(IN) :: VREL
! Slip velocity between wall and particles
      DOUBLE PRECISION, INTENT(IN) :: VSLIP
! Relevant solids velocity at wall
      DOUBLE PRECISION, INTENT(IN) :: VEL
! Relevant wall velocity
      DOUBLE PRECISION, INTENT(IN) :: W_VEL
! del.u
      DOUBLE PRECISION, INTENT(IN) :: DEL_U
! S:S
      DOUBLE PRECISION, INTENT(IN) :: S_S
! S_dd (d can be x,y or z)
      DOUBLE PRECISION, INTENT(IN) :: S_dd
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Wall momentum coefficients:
! 1st term on LHS
      DOUBLE PRECISION, INTENT(INOUT) :: Gw
! 2nd term on LHS
      DOUBLE PRECISION, INTENT(INOUT) :: Hw
! all terms appearing on RHS
      DOUBLE PRECISION, INTENT(INOUT) :: Cw
!-----------------------------------------------
! Local Variables      
!-----------------------------------------------
! Part of wall momentum coefficient in 1st term of LHS
      DOUBLE PRECISION :: Mu_s
! Term appearing in wall momenutm coefficient of 2nd term LHS
      DOUBLE PRECISION :: F_2
! Frictional Pressure Pf
      DOUBLE PRECISION :: Pf
! Critical pressure Pc
      DOUBLE PRECISION :: Pc
! Chi (appears in frictional boundary condition)
      DOUBLE PRECISION :: Chi, N_Pff
! Viscosity
      DOUBLE PRECISION :: Mu
! Bulk viscosity
      DOUBLE PRECISION :: Mu_b
! Viscosity corrected for interstitial fluid effects
      DOUBLE PRECISION :: Mu_star
! Reynolds number based on slip velocity
      DOUBLE PRECISION :: Re_g
! Friction Factor in drag coefficient
      DOUBLE PRECISION :: C_d
! Drag Coefficient
      DOUBLE PRECISION :: Beta, DgA
! Square root of S:S or the form suggested by Savage
      DOUBLE PRECISION :: ZETA
! Constants in Simonin model
      DOUBLE PRECISION :: Sigma_c, Tau_2_c, Tau_12_st, Nu_t
      DOUBLE PRECISION :: Tau_2, zeta_c_2, MU_2_T_Kin, Mu_2_Col
      DOUBLE PRECISION :: Tmp_Ahmadi_Const
! Error message
      CHARACTER*80     LINE
! Other local terms
      DOUBLE PRECISION :: dpc_dphi
! Local variable for rdf      
      DOUBLE PRECISION :: g0_loc
!----------------------------------------------- 
      
! This is done here similar to bc_theta to avoid small negative values of
! Theta coming most probably from linear solver
      IF(TH .LE. ZERO)THEN
        TH = 1D-8
        IF (myPE.eq.PE_IO) THEN
           WRITE(*,*) &
              'Warning: Negative granular temp at wall set to 1e-8'
!          CALL WRITE_ERROR('THETA_HW_CW', LINE, 1)
        ENDIF
      ENDIF

      g0_loc = g0
      
! modify F_2 if Jenkins BC is used (sof)    
      IF(JENKINS) THEN

        IF (VSLIP == ZERO) THEN
! if solids velocity field is initialized to zero, use free slip bc
           F_2 = zero

        ELSE
           IF(AHMADI) THEN
! Ahmadi model uses different solids pressure model
! the coefficient mu in Jenkins paper is defined as tan_Phi_w, that's how
! I understand it from soil mechanic papers, i.e., G.I. Tardos, powder
! Tech. 92 (1997), 61-74. See his equation (1). Define Phi_w in mfix.dat!
! here F_2 divided by VSLIP to use the same bc as Johnson&Jackson
              F_2 = tan_Phi_w*RO_s(M)*EPS* &
                 ((ONE + 4.0D0*g0EP_avg) + HALF*(ONE -C_e*C_e))*TH/VSLIP

           ELSE
! Simonin or granular models use same solids pressure
              F_2 = tan_Phi_w*RO_s(M)*EPS*(1d0+ 4.D0 * Eta *g0EP_avg)*TH/VSLIP
           ENDIF   ! end if(Ahmadi)

        ENDIF ! endif(vslip==0)

      ELSE    ! if(.not.jenkins)

         F_2 = (PHIP*DSQRT(3d0*TH)*Pi*RO_s(M)*EPS*g0_loc)&
              /(6d0*(ONE-ep_star_avg))

      ENDIF   ! end if(Jenkins)/else 
 

      Mu = (5d0*DSQRT(Pi*TH)*Dp_avg*RO_s(M))/96d0 
      Mu_b = (256d0*Mu*EPS*g0EP_avg)/(5d0*Pi) 

! This is from Wen-Yu correlation, you can put here your own single particle drag 
      Re_g = EPG*RO_g_avg*Dp_avg*VREL/Mu_g_avg
      IF (Re_g.lt.1000d0) THEN
         C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
      ELSE
         C_d = 0.44d0
      ENDIF
      DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(Dp_avg*EPG**(2.65d0))
      IF(VREL == ZERO) DgA = LARGE_NUMBER
      Beta = SWITCH*EPS*DgA
 
! SWITCH enables us to turn on/off the modification to the
! particulate phase viscosity. If we want to simulate gas-particle
! flow then SWITCH=1 to incorporate the effect of drag on the
! particle viscosity. If we want to simulate granular flow
! without the effects of an interstitial gas, SWITCH=0.
      IF(SWITCH == ZERO .OR. Ro_g_avg == ZERO)THEN
         Mu_star = Mu
      ELSEIF(TH .LT. SMALL_NUMBER)THEN
         MU_star = ZERO
      ELSE
         Mu_star = RO_S(M)*EPS* g0_loc*TH* Mu/ &
            (RO_S(M)*g0EP_avg*TH + 2.0d0*SWITCH*DgA/RO_S(M)* Mu)
      ENDIF

      Mu_s = ((2d0+ALPHA)/3d0)*((Mu_star/(Eta*(2d0-Eta)*&
                   g0_loc))*(ONE+1.6d0*Eta*g0EP_avg&
                   )*(ONE+1.6d0*Eta*(3d0*Eta-2d0)*&
                   g0EP_avg)+(0.6d0*Mu_b*Eta))

! particle relaxation time
      Tau_12_st = RO_s(M)/(DgA+small_number)

      IF(SIMONIN) THEN !see calc_mu_s for explanation of these definitions

         Sigma_c = (ONE+ C_e)*(3.d0-C_e)/5.d0
         Tau_2_c = DP_avg/(6.d0*EPS*g0_loc*DSQRT(16.d0*(TH+Small_number)/PI))
         zeta_c_2= 2.D0/5.D0*(ONE+ C_e)*(3.d0*C_e-ONE)
         Nu_t =  Tau_12_avg/Tau_12_st
         Tau_2 = ONE/(2.D0/Tau_12_st+Sigma_c/Tau_2_c)
         MU_2_T_Kin = (2.0D0/3.0D0*K_12_avg*Nu_t + TH * &
                     (ONE+ zeta_c_2*EPS*g0_loc))*Tau_2
         Mu_2_Col = 8.D0/5.D0*EPS*g0_loc*Eta* (MU_2_T_Kin+ &
                   Dp_avg*DSQRT(TH/PI))
         Mu_s = EPS*RO_s(M)*(MU_2_T_Kin + Mu_2_Col)
      ELSEIF(AHMADI) THEN
         IF(EPS < (ONE-ep_star_avg)) THEN
            Tmp_Ahmadi_Const = &
            ONE/(ONE+ Tau_1_avg/Tau_12_st * (ONE-EPS/(ONE-ep_star_avg))**3)
         ELSE
            Tmp_Ahmadi_Const = ONE
         ENDIF
         Mu_s = Tmp_Ahmadi_Const &
             *0.1045D0*(ONE/g0_loc+3.2D0*EPS+12.1824D0*g0_loc*EPS*EPS)  &
             *Dp_avg*RO_s(M)* DSQRT(TH)
      ENDIF
 
! Calculating frictional terms
      IF ((ONE-EPG)<= EPS_f_min) THEN
         Pf = ZERO
         Chi = ZERO
         ZETA = 1d0
      ELSE
         IF (SAVAGE.EQ.1) THEN    !form of Savage
            ZETA = ((48d0*Eta*(1d0-Eta)*RO_s(M)*EPS*EPS*g0_loc*&
                    (TH**1.5d0))/&
                    (SQRT_Pi*Dp_avg*2d0*Mu_s))**0.5d0 
         ELSEIF (SAVAGE.EQ.0)  THEN !S:S form
            ZETA = DSQRT(S_S) 
         ELSE
            ZETA = DSQRT(S_S + (TH/(Dp_avg*Dp_avg)))
         ENDIF
  
         IF (EPG < ep_star_avg) THEN
            dpc_dphi = (to_SI*Fr)*((delta**5)*(2d0*(ONE-ep_star_avg-delta) - &
               2d0*eps_f_min)+((ONE-ep_star_avg-delta)-eps_f_min)&
               *(5*delta**4))/(delta**10)

            Pc = (to_SI*Fr)*(((ONE-ep_star_avg-delta) - EPS_f_min)**N_Pc)/(delta**D_Pc)
!            Pc=  1d25*(((ONE-EPG)-(ONE-ep_star_avg))**10d0)  ! this is old Pc
         ELSE
            Pc = Fr*(((ONE-EPG) - EPS_f_min)**N_Pc)/ &
               (((ONE-ep_star_avg) - (ONE-EPG) + SMALL_NUMBER)**D_Pc)
         ENDIF
 
         IF (DEL_U .GE. ZERO) THEN
            N_Pff = DSQRT(3d0)/(2d0*Sin_Phi) !dilatation
         ELSE
            N_Pff = N_Pf !compaction
         ENDIF
  
         IF ((DEL_U/(ZETA*N_Pff*DSQRT(2d0)*Sin_Phi)) .GT. 1d0) THEN
            Pf = ZERO
         ELSEIF( DEL_U == ZERO ) THEN
            Pf = Pc
         ELSE
            Pf = Pc*(1d0 - (DEL_U/(ZETA*N_Pff*DSQRT(2d0)*Sin_Phi)))**&
               (N_Pff-1d0)
         ENDIF
 
         Chi =  DSQRT(2d0)*Pf*Sin_Phi*(N_Pff - (N_Pff-1d0)*&
            ((Pf/Pc)**(1d0/(N_Pff-1d0))))
  
         IF (Chi< ZERO) THEN
            Pf = Pc*((N_Pff/(N_Pff-1d0))**(N_Pff-1d0))
            Chi = ZERO
         ENDIF

! by writing Pf & Chi in the following form, we ensure distribution 
! of stresses amoung all solids phases (sof, Oct 24 2005)

         Pf = Pf * EPS/(ONE-EPG)
         Chi = Chi * EPS/(ONE-EPG)
 
      ENDIF
 
! Calculating gw, hw, cw
 
      Gw = (MU_s + Chi/(2d0*ZETA))*DABS(VEL - W_VEL)
 
      Hw = F_2*DABS(VEL - W_VEL) + Pf*tan_Phi_w
      IF(ZETA .NE. ZERO) Hw = Hw - Chi*S_dd*tan_Phi_w/ZETA
      
      Cw = hw * W_VEL
 
      RETURN
      END SUBROUTINE CALC_Gw_Hw_Cw
