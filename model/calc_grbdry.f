!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_GRBDRY
!  Purpose: Main controller subroutine for calculating coefficients
!   for momentum boundary conditions using kinetic & frictional theory
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

      SUBROUTINE CALC_GRBDRY(IJK1, IJK2, FCELL, COM, M, L, Gw, Hw, Cw)

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
      DOUBLE PRECISION :: EPg_avg, Mu_g_avg, RO_g_avg
! Average void fraction at packing
      DOUBLE PRECISION :: ep_star_avg
! Average scalars modified to include all solid phases
      DOUBLE PRECISION :: EPs_avg(DIMENSION_M), DP_avg(DIMENSION_M),&
                          TH_avg(DIMENSION_M), AVGX1, AVGX2, smallTheta
! Average Simonin and Ahmadi variables (sof)
      DOUBLE PRECISION :: K_12_avg, Tau_12_avg, Tau_1_avg
! Average velocities
! values of U_sm, V_sm, W_sm at appropriate place on boundary wall
      DOUBLE PRECISION :: USCM, VSCM,WSCM
      DOUBLE PRECISION :: USCM1,USCM2,VSCM1,VSCM2,WSCM1,WSCM2
! values of U_g, V_g, W_g at appropriate place on boundary wall
      DOUBLE PRECISION :: UGC, VGC, WGC
      DOUBLE PRECISION :: WGC1, WGC2, VGC1, VGC2, UGC1, UGC2
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
! Average Radial distribution function
      DOUBLE PRECISION :: g_0AVG
! radial distribution function at contact
      DOUBLE PRECISION :: g0(DIMENSION_M)
! Sum of eps*G_0
      DOUBLE PRECISION :: g0EPs_avg 
! Error message
      CHARACTER*80     LINE
 
!-----------------------------------------------
!  Function subroutines
!-----------------------------------------------
      DOUBLE PRECISION F_HW
!----------------------------------------------- 
! Include statements functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

! Note: EP_s, MU_g, and RO_g are undefined at IJK1 (wall cell).
!       Hence IJK2 (fluid cell) is used in averages.

      smallTheta = (to_SI)**4 * ZERO_EP_S

!-----------------------------------------------      
! Calculations for U momentum equation
      IF (COM .EQ. 'U')THEN

        IJK2E = EAST_OF(IJK2)
        IPJK2 = IP_OF(IJK2)

        EPg_avg = AVG_X(EP_g(IJK2), EP_g(IJK2E), I_OF(IJK2))
        ep_star_avg = AVG_X(EP_star_array(IJK2), EP_star_array(IJK2E), I_OF(IJK2))
        Mu_g_avg = AVG_X(Mu_g(IJK2), Mu_g(IJK2E), I_OF(IJK2))
        RO_g_avg = AVG_X(RO_g(IJK2), RO_g(IJK2E), I_OF(IJK2))

        IF(SIMONIN .OR. AHMADI) THEN
           K_12_avg = AVG_X(K_12(IJK2), K_12(IJK2E), I_OF(IJK2))
           Tau_12_avg = AVG_X(Tau_12(IJK2), Tau_12(IJK2E), I_OF(IJK2))
           Tau_1_avg = AVG_X(Tau_1(IJK2), Tau_1(IJK2E), I_OF(IJK2))
        ELSE
           K_12_avg = ZERO    
           Tau_12_avg = ZERO
           Tau_1_avg = ZERO
        ENDIF

        g0EPs_avg = ZERO             
        DO MM = 1, SMAX
          g0(MM)      = G_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), M, MM)
          EPs_avg(MM) = AVG_X(EP_s(IJK2, MM), EP_s(IJK2E, MM), I_OF(IJK2))
          DP_avg(MM)  = AVG_X(D_P(IJK2,MM), D_P(IJK2E,MM), I_OF(IJK2))
          g0EPs_avg   = g0EPs_avg + G_0AVG(IJK2, IJK2E, 'X', I_OF(IJK2), M, MM) &
             * AVG_X(EP_s(IJK2, MM), EP_s(IJK2E, MM), I_OF(IJK2))
          IF (.NOT.GRANULAR_ENERGY) THEN
            TH_avg(MM) = AVG_X(THETA_M(IJK2,MM), THETA_M(IJK2E,MM), I_OF(IJK2))
            IF(TH_avg(MM) < ZERO) TH_avg(MM) = smallTheta
          ENDIF
        ENDDO

        WVELS = BC_Uw_s(L,M)


        IF(FCELL .EQ. 'N')THEN
          IPJMK2 = JM_OF(IPJK2) 

! code modified for some corner cells
          IF(GRANULAR_ENERGY) THEN
            DO MM = 1, SMAX
              AVGX1 = AVG_X(Theta_m(IJK1,MM), Theta_m(IPJMK2,MM), I_OF(IJK1))
              AVGX2 = AVG_X(Theta_m(IJK2,MM), Theta_m(IPJK2,MM), I_OF(IJK2))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Y(AVGX1, AVGX2, J_OF(IJK1))
              ENDIF
            ENDDO
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

! code modified for some corner cells          
          IF(GRANULAR_ENERGY) THEN
            DO MM = 1, SMAX                  
              AVGX1 = AVG_X(Theta_m(IJK2,MM),Theta_m(IPJK2,MM),I_OF(IJK2))
              AVGX2 = AVG_X(Theta_m(IJK1,MM),Theta_m(IPJPK2,MM),I_OF(IJK1))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Y(AVGX1, AVGX2, J_OF(IJK2))
              ENDIF
            ENDDO
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

! code modified for some corner cells          
          IF(GRANULAR_ENERGY) THEN          
            DO MM = 1, SMAX
              AVGX1 = AVG_X(Theta_m(IJK1,MM),Theta_m(IPJKM2,MM),I_OF(IJK1))
              AVGX2 = AVG_X(Theta_m(IJK2,MM),Theta_m(IPJK2,MM),I_OF(IJK2))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Z(AVGX1, AVGX2, K_OF(IJK1))
              ENDIF
            ENDDO
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

! code modified for some corner cells
          IF(GRANULAR_ENERGY) THEN
            DO MM = 1, SMAX                  
              AVGX1 = AVG_X(Theta_m(IJK1,MM), Theta_m(IPJKP2,MM),I_OF(IJK1))
              AVGX2 = AVG_X(Theta_m(IJK2,MM), Theta_m(IPJK2,MM),I_OF(IJK2))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Z(AVGX1, AVGX2, K_OF(IJK2))
              ENDIF
            ENDDO
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
          CALL WRITE_ERROR('CALC_GRBDRY', LINE, 1)
          CALL exitMPI(myPE)                    
        ENDIF
 
!-----------------------------------------------      
! Calculations for V momentum equation
      ELSEIF (COM .EQ. 'V')THEN

        IJK2N = NORTH_OF(IJK2)
        IJPK2 = JP_OF(IJK2)

        EPg_avg = AVG_Y(EP_g(IJK2), EP_g(IJK2N), J_OF(IJK2))
        ep_star_avg = AVG_Y(EP_star_array(IJK2), EP_star_array(IJK2N), J_OF(IJK2))
        Mu_g_avg = AVG_Y(Mu_g(IJK2), Mu_g(IJK2N), J_OF(IJK2))
        RO_g_avg = AVG_Y(RO_g(IJK2), RO_g(IJK2N), J_OF(IJK2))

        IF(SIMONIN .OR. AHMADI) THEN
           K_12_avg = AVG_Y(K_12(IJK2), K_12(IJK2N), J_OF(IJK2))
           Tau_12_avg = AVG_Y(Tau_12(IJK2), Tau_12(IJK2N), J_OF(IJK2))
           Tau_1_avg = AVG_Y(Tau_1(IJK2), Tau_1(IJK2N), J_OF(IJK2))
        ELSE
           K_12_avg = ZERO    
           Tau_12_avg = ZERO
           Tau_1_avg = ZERO
        ENDIF

        g0EPs_avg = ZERO
        DO MM = 1, SMAX
           g0(MM)      = G_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), M, MM)
           EPs_avg(MM) = AVG_Y(EP_s(IJK2, MM), EP_s(IJK2N, MM), J_OF(IJK2))
           DP_avg(MM)  = AVG_Y(D_p(IJK2,MM), D_p(IJK2N,MM), J_OF(IJK2))
           g0EPs_avg   = g0EPs_avg + G_0AVG(IJK2, IJK2N, 'Y', J_OF(IJK2), M, MM) &
                        * AVG_Y(EP_s(IJK2, MM), EP_s(IJK2N, MM), J_OF(IJK2))
          IF(.NOT.GRANULAR_ENERGY) THEN                        
            TH_avg(MM) = AVG_Y(THETA_M(IJK2,MM), THETA_M(IJK2N,MM), J_OF(IJK2))
            IF(TH_avg(MM) < ZERO) TH_avg(MM) = smallTheta
          ENDIF
        ENDDO

        WVELS = BC_Vw_s(L, M)

        IF(FCELL .EQ. 'T')THEN
          IJPKM2 = KM_OF(IJPK2)

          IF(GRANULAR_ENERGY) THEN
            DO MM = 1, SMAX                    
              AVGX1 = AVG_Z(Theta_m(IJK1,MM), Theta_m(IJK2,MM), K_OF(IJK1))
              AVGX2 = AVG_Z(Theta_m(IJPKM2,MM), Theta_m(IJPK2,MM), K_OF(IJPKM2))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Y(AVGX1, AVGX2, J_OF(IJK1))
              ENDIF
            ENDDO
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
            DO MM = 1, SMAX
              AVGX1 = AVG_Z(Theta_m(IJK2,MM), Theta_m(IJK1,MM), K_OF(IJK2))
              AVGX2 = AVG_Z(Theta_m(IJPK2,MM), Theta_m(IJPKP2,MM), K_OF(IJPK2))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Y(AVGX1, AVGX2, J_OF(IJK2))
              ENDIF
            ENDDO
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
          VELS = VSCM 
 
        ELSEIF(FCELL .EQ. 'E')THEN
          IMJPK2= IM_OF(IJPK2)

          IF(GRANULAR_ENERGY) THEN
            DO MM = 1, SMAX
              AVGX1 = AVG_X(Theta_m(IJK1,MM),Theta_m(IJK2,MM),I_OF(IJK1))
              AVGX2 = AVG_X(Theta_m(IMJPK2,MM),Theta_m(IJPK2,MM),I_OF(IMJPK2))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Y(AVGX1, AVGX2, J_OF(IJK1))
              ENDIF
            ENDDO
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
            DO MM = 1, SMAX
              AVGX1 = AVG_X(Theta_m(IJK2,MM),Theta_m(IJK1,MM),I_OF(IJK2))
              AVGX2 = AVG_X(Theta_m(IJPK2,MM),Theta_m(IPJPK2,MM),I_OF(IJPK2))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Y(AVGX1, AVGX2, J_OF(IJK2))
              ENDIF
            ENDDO
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
          CALL WRITE_ERROR('CALC_GRBDRY', LINE, 1)
          CALL exitMPI(myPE)                    
        ENDIF


!-----------------------------------------------      
! Calculations for W momentum equation
      ELSEIF (COM .EQ. 'W')THEN
        IJK2T = TOP_OF(IJK2)
        IJKP2 = KP_OF(IJK2)

        EPg_avg = AVG_Z(EP_g(IJK2), EP_g(IJK2T), K_OF(IJK2))
        Mu_g_avg = AVG_Z(Mu_g(IJK2), Mu_g(IJK2T), K_OF(IJK2))
        RO_g_avg = AVG_Z(RO_g(IJK2), RO_g(IJK2T), K_OF(IJK2))
        ep_star_avg = AVG_Z(EP_star_array(IJK2), EP_star_array(IJK2T), K_OF(IJK2))

        IF(SIMONIN .OR. AHMADI) THEN
           K_12_avg = AVG_Z(K_12(IJK2), K_12(IJK2T), K_OF(IJK2))
           Tau_12_avg = AVG_Z(Tau_12(IJK2), Tau_12(IJK2T), K_OF(IJK2))
           Tau_1_avg = AVG_Z(Tau_1(IJK2), Tau_1(IJK2T), K_OF(IJK2))
        ELSE
           K_12_avg = ZERO    
           Tau_12_avg = ZERO
           Tau_1_avg = ZERO
        ENDIF

        g0EPs_avg = ZERO
        DO MM = 1, SMAX
           g0(MM)      = G_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, MM)
           EPs_avg(MM) = AVG_Z(EP_s(IJK2,MM), EP_s(IJK2T,MM), K_OF(IJK2))
           DP_avg(MM)  = AVG_Z(D_p(IJK2,MM), D_p(IJK2T,MM), K_OF(IJK2))
           g0EPs_avg   = g0EPs_avg + G_0AVG(IJK2, IJK2T, 'Z', K_OF(IJK2), M, MM) &
                       * AVG_Z(EP_s(IJK2, MM), EP_s(IJK2T, MM), K_OF(IJK2))
           IF(.NOT.GRANULAR_ENERGY) THEN
              TH_avg(MM) = AVG_Z(THETA_M(IJK2,MM), THETA_M(IJK2T,MM), K_OF(IJK2))
              IF(TH_avg(MM) < ZERO) TH_avg(MM) = smallTheta
           ENDIF
        ENDDO

        WVELS = BC_Ww_s(L,M)        

        IF(FCELL .EQ. 'N')THEN
          IJMKP2 = JM_OF(IJKP2)

          IF(GRANULAR_ENERGY) THEN
            DO MM = 1, SMAX
               AVGX1 = AVG_Z(Theta_m(IJK1,MM), Theta_m(IJMKP2,MM), K_OF(IJK1))
               AVGX2 = AVG_Z(Theta_m(IJK2,MM), Theta_m(IJKP2,MM), K_OF(IJK2))
               IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
               IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
               IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                 TH_avg(MM) = smallTheta
               ELSE
                 TH_avg(MM) = AVG_Y(AVGX1, AVGX2, J_OF(IJK1))
               ENDIF
            ENDDO
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
            DO MM = 1, SMAX
              AVGX1 = AVG_Z(Theta_m(IJK2,MM), Theta_m(IJKP2,MM), K_OF(IJK2))
              AVGX2 = AVG_Z(Theta_m(IJK1,MM), Theta_m(IJPKP2,MM), K_OF(IJK1))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Y(AVGX1, AVGX2, J_OF(IJK2))
              ENDIF
            ENDDO
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
            DO MM = 1, SMAX
              AVGX1 = AVG_X(Theta_m(IJK1,MM),Theta_m(IJK2,MM),I_OF(IJK1))
              AVGX2 = AVG_X(Theta_m(IMJKP2,MM),Theta_m(IJKP2,MM),I_OF(IMJKP2))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Z(AVGX1, AVGX2, K_OF(IJK1))
              ENDIF
            ENDDO
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
            DO MM = 1, SMAX                  
              AVGX1 = AVG_X(Theta_m(IJK2,MM),Theta_m(IJK1,MM),I_OF(IJK2))
              AVGX2 = AVG_X(Theta_m(IJKP2,MM),Theta_m(IPJKP2,MM),I_OF(IJKP2))
              IF(AVGX1 < ZERO .AND. AVGX2 > ZERO) AVGX1 = AVGX2
              IF(AVGX2 < ZERO .AND. AVGX1 > ZERO) AVGX2 = AVGX1
              IF(AVGX1 < ZERO .AND. AVGX2 < ZERO) THEN
                TH_avg(MM) = smallTheta
              ELSE
                TH_avg(MM) = AVG_Z(AVGX1, AVGX2, K_OF(IJK2))
              ENDIF
            ENDDO
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
          CALL WRITE_ERROR('CALC_GRBDRY', LINE, 1)
          CALL exitMPI(myPE)
        ENDIF
 

      ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown COM'
         CALL WRITE_ERROR('CALC_GRBDRY', LINE, 1)
         CALL exitMPI(myPE)
      ENDIF

! magnitude of gas-solids relative velocity
      VREL =&
         DSQRT( (UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2 )

! slip velocity for use in Jenkins bc (sof)	  
      VSLIP= DSQRT( (USCM-BC_UW_S(L,M))**2 + (VSCM-BC_VW_S(L,M))**2 &
         + (WSCM-BC_WW_S(L,M))**2 )
 
      IF (FRICTION .AND. (ONE-EP_G(IJK2))>EPS_F_MIN) THEN         
        CALL CALC_S_DDOT_S(IJK1, IJK2, FCELL, COM, M, DEL_DOT_U,&
           S_DDOT_S, S_dd)
 
        CALL CALC_Gw_Hw_Cw(g0(M), EPs_avg(M), EPg_avg, ep_star_avg, &
           g0EPs_avg, TH_avg(M), Mu_g_avg, RO_g_avg, &
           DP_avg(M), K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP,&
           DEL_DOT_U, S_DDOT_S, S_dd, VELS, WVELS, M, gw, hw, cw)
      ELSE
         GW = 1D0               
         Hw = F_Hw(g0, EPs_avg, EPg_avg, ep_star_avg, &
            g0EPs_avg, TH_avg, Mu_g_avg, RO_g_avg, &
            DP_avg, K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP, M)
         CW = HW*WVELS            
      ENDIF


      RETURN
      END SUBROUTINE CALC_GRBDRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: F_HW

!  Purpose: Function for hw                                            C
!                                                                      C
!  Author: K. Agrawal & A. Srivastava, Princeton Univ. Date: 24-JAN-98 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
!  Modified: Sofiane Benyahia, Fluent Inc.             Date: 02-FEB-05 C
!  Purpose: Include conductivity defined by Simonin and Ahmadi         C
!           Also included Jenkins small frictional limit               C
!                                                                      C
!  Literature/Document References: See calcmu_s.f for ref. on Simonin  C
!  and Ahmadi models; for Jenkins BC: Jenkins and Louge, Phys. fluids  C
!  9 (10), 2835. See equation (2) in the paper                         C
!                                                  

!  Additional Notes:
!    The current implementation of the IA (2005) kinetic theory and 
!    the GD (1999) kinetic theory do not incorporate ahmadi or simonin
!    additions nor jenkins small frictional bc mdoel

!    The granular momentum BC is written as the normal vector dot the
!    stress tensor.  Besides the gradient in velocity of phase M, the
!    stress tensor expression may contain several additional terms that 
!    would need to be accounted for when satisfying the BC.  These
!    modifications have NOT been rigorously addressed.

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION F_HW(g0,EPS,EPG, ep_star_avg, &
                                     g0EPs_avg,TH,Mu_g_avg,RO_g_avg,&
                                     DP_avg,K_12_avg, Tau_12_avg, Tau_1_avg, &
                                     VREL, VSLIP, M)
 
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE constant
      USE physprop
      USE run
      USE fldvar
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------      
! Radial distribution function of solids phase M with each
! other solids phase 
      DOUBLE PRECISION, INTENT(IN) :: g0(DIMENSION_M) 
! Average solids volume fraction of each solids phase
      DOUBLE PRECISION, INTENT(IN) :: EPS(DIMENSION_M)
! Average solids and gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPG, ep_star_avg
! Sum of eps*G_0 
      DOUBLE PRECISION, INTENT(INOUT) :: g0EPs_avg 
! Average theta_m
      DOUBLE PRECISION, INTENT(INOUT) :: TH (DIMENSION_M)      
! Average gas viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mu_g_avg
! Average gas density
      DOUBLE PRECISION, INTENT(IN) :: RO_g_avg
! Average particle diameter of each solids phase
      DOUBLE PRECISION, INTENT(IN) :: DP_avg(DIMENSION_M)
! Average cross-correlation K_12 and lagrangian integral time-scale
      DOUBLE PRECISION, INTENT(IN) :: K_12_avg, Tau_12_avg, Tau_1_avg
! Magnitude of slip velocity between two phases
      DOUBLE PRECISION, INTENT(IN) :: VREL
! Slip velocity between wall and particles
      DOUBLE PRECISION, INTENT(IN) :: VSLIP
! Solids phase index
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local Variables      
!-----------------------------------------------
! Solids phase index
      INTEGER :: LL
! Coefficient of 2nd term
      DOUBLE PRECISION :: F_2
! Coefficient of 1st term
      DOUBLE PRECISION :: Mu_s
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
! Constants in Simonin or Ahmadi model
      DOUBLE PRECISION :: Sigma_c, Tau_2_c, Tau_12_st, Nu_t
      DOUBLE PRECISION :: Tau_2, zeta_c_2, MU_2_T_Kin, Mu_2_Col
      DOUBLE PRECISION :: Tmp_Ahmadi_Const
! Variables for Iddir & Arastoopour model
      DOUBLE PRECISION :: MU_sM_sum, MU_s_MM, MU_s_LM, MU_sM_ip, MU_common_term,&
                          MU_sM_LM
      DOUBLE PRECISION :: M_PM, M_PL, MPSUM, NU_PL, NU_PM, D_PM, D_PL, DPSUMo2
      DOUBLE PRECISION :: Ap_lm, Dp_lm, R1p_lm, Bp_lm
! Variables for Garzo and Dufty model
      DOUBLE PRECISION :: c_star, zeta0_star, nu_eta_star, press_star, &
                          gamma_star, eta_k_star, eta_star, eta0
!----------------------------------------------- 

! This is done here similar to bc_theta to avoid small negative values of
! Theta coming most probably from linear solver
      IF(TH(M) .LE. ZERO)THEN
        TH(M) = 1D-8
        IF (myPE.eq.PE_IO) THEN
           WRITE(*,*) &
             'Warning: Negative granular temp at wall set to 1e-8'
!          CALL WRITE_ERROR('THETA_HW_CW', LINE, 1)
        ENDIF
      ENDIF

! In F_2 and Mu a DSQRT(T) has been left out as it appears in both
! terms and thus cancels out upon dividing the former by the latter
! The above statement was not implemented because Simonin viscosity
! doesn't have a sqrt(th) directly available to use this simplification.


      IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN  

! Use original IA theory if SWITCH_IA is false
         IF(.NOT. SWITCH_IA) g0EPs_avg = EPS(M)*RO_S(M)

         D_PM = DP_avg(M)
         M_PM = (PI/6.d0)*(D_PM**3)*RO_S(M)
         NU_PM = (EPS(M)*RO_S(M))/M_PM
 
         F_2 = (PHIP*DSQRT(3.d0*TH(M)/M_PM)*PI*RO_s(M)*EPS(M)*g0(M))/&
            (6.d0*(ONE-ep_star_avg))

! This is from Wen-Yu correlation, you can put here your own single particle drag
         Re_g = EPG*RO_g_avg*D_PM*VREL/Mu_g_avg
         IF (Re_g .lt. 1000.d0) THEN
            C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
         ELSE
            C_d = 0.44d0
         ENDIF
         DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(D_PM*EPG**(2.65d0))
         IF(VREL == ZERO) DgA = LARGE_NUMBER
         Beta = EPS(M)*DgA !this is equivalent to F_gs(ijk,m)

         Mu = (5.d0/96.d0)*D_PM* RO_S(M)*DSQRT(PI*TH(M)/M_PM)

         IF(.NOT.SWITCH_IA .OR. RO_g_avg == ZERO)THEN
            Mu_star = Mu
         ELSEIF(TH(M) .LT. SMALL_NUMBER)THEN
            MU_star = ZERO
         ELSE
            Mu_star = Mu*EPS(M)*g0(M)/ (g0EPs_avg+ 2.0d0*DgA*Mu / &
               (RO_S(M)**2 *(TH(M)/M_PM)))
         ENDIF

         MU_s_MM = (Mu_star/g0(M))*(1.d0+(4.d0/5.d0)*(1.d0+C_E)*g0EPs_avg)**2
         Mu_sM_sum = ZERO

         DO LL = 1, SMAX
            D_PL = DP_avg(LL)
            M_PL = (PI/6.d0)*(D_PL**3.)*RO_S(LL)
            MPSUM = M_PM + M_PL
            DPSUMo2 = (D_PM+D_PL)/2.d0
            NU_PL = (EPS(LL)*RO_S(LL))/M_PL

            IF ( LL .eq. M) THEN
               Ap_lm = MPSUM/(2.d0)
               Dp_lm = M_PL*M_PM/(2.d0*MPSUM)
               R1p_lm = ONE/( Ap_lm**1.5 * Dp_lm**3 )
               MU_s_LM = DSQRT(PI)*(DPSUMo2**4 / (48.d0*5.d0))*&
                  g0(LL)*(M_PL*M_PM/MPSUM)*(M_PL*M_PM/&
                  MPSUM)*((M_PL*M_PM)**1.5)*NU_PM*NU_PL*&
                  (1.d0+C_E)*R1p_lm*DSQRT(TH(M))

! solids phase 'viscosity' associated with the divergence 
! of solids phase M
               MU_sM_ip = (MU_s_MM + MU_s_LM)

            ELSE
               Ap_lm = (M_PM*TH(LL)+M_PL*TH(M))/&
                  (2.d0)
               Bp_lm = (M_PM*M_PL*(TH(LL)-TH(M) ))/&
                  (2.d0*MPSUM)
               Dp_lm = (M_PL*M_PM*(M_PM*TH(M)+M_PL*TH(LL) ))/&
                  (2.d0*MPSUM*MPSUM)
               R1p_lm = ( ONE/( Ap_lm**1.5 * Dp_lm**3 ) ) + &
                  ( (9.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**4 ) )+&
                  ( (30.d0*Bp_lm**4)/( 2.d0*Ap_lm**3.5 * Dp_lm**5 ) )
               MU_common_term = DSQRT(PI)*(DPSUMo2**4 / (48.d0*5.d0))*&
                  g0(LL)*(M_PL*M_PM/MPSUM)*(M_PL*M_PM/&
                  MPSUM)*((M_PL*M_PM)**1.5)*NU_PM*NU_PL*&
                  (1.d0+C_E)*R1p_lm
               MU_sM_LM = MU_common_term*(TH(M)*TH(M)*TH(LL)*TH(LL)*TH(LL) )

! solids phase 'viscosity' associated with the divergence
! of solids phase M       
               MU_sM_ip = MU_sM_LM

            ENDIF
            MU_sM_sum = MU_sM_sum + MU_sM_ip

          ENDDO

! Find the term proportional to the gradient in velocity
! of phase M  (viscosity in the Mth solids phase)
          Mu_s = MU_sM_sum

      ELSEIF (TRIM(KT_TYPE) .EQ. 'GD_99') THEN

         D_PM = DP_avg(M)
         M_PM = (PI/6.d0)*(D_PM**3)*RO_S(M)
         NU_PM = (EPS(M)*RO_S(M))/M_PM

         F_2 = (PHIP*DSQRT(3.d0*TH(M))*PI*RO_s(M)*EPS(M)*g0(M))/&
            (6.d0*(ONE-ep_star_avg))

! This is from Wen-Yu correlation, you can put here your own single particle drag
         Re_g = EPG*RO_g_avg*D_PM*VREL/Mu_g_avg
         IF (Re_g .lt. 1000.d0) THEN
            C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
         ELSE
            C_d = 0.44d0
         ENDIF
         DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(D_PM*EPG**(2.65d0))
         IF(VREL == ZERO) DgA = LARGE_NUMBER
         Beta = EPS(M)*DgA !this is equivalent to F_gs(ijk,m)
  
! Pressure/Viscosity/Bulk Viscosity
! Note: k_boltz = M_PM
     
! Find pressure in the Mth solids phase
         press_star = 1.d0 + 2.d0*(1.d0+C_E)*EPS(M)*G0(M)
 
! find bulk and shear viscosity
         c_star = 32.0d0*(1.0d0 - C_E)*(1.d0 - 2.0d0*C_E*C_E) &
            / (81.d0 - 17.d0*C_E + 30.d0*C_E*C_E*(1.0d0-C_E))

         zeta0_star = (5.d0/12.d0)*G0(M)*(1.d0 - C_E*C_E) &
            * (1.d0 + (3.d0/32.d0)*c_star)

         nu_eta_star = G0(M)*(1.d0 - 0.25d0*(1.d0-C_E)*(1.d0-C_E)) &
            * (1.d0-(c_star/64.d0))

         gamma_star = (4.d0/5.d0)*(32.d0/PI)*EPS(M)*EPS(M) &
            * G0(M)*(1.d0+C_E) * (1.d0 - (c_star/32.d0))

         eta_k_star = (1.d0 - (2.d0/5.d0)*(1.d0+C_E)*(1.d0-3.d0*C_E) &
            * EPS(M)*G0(M) ) / (nu_eta_star - 0.5d0*zeta0_star)

         eta_star = eta_k_star*(1.d0 + (4.d0/5.d0)*EPS(M)*G0(M) &
            * (1.d0+C_E) ) + (3.d0/5.d0)*gamma_star

         eta0 = 5.0d0*M_PM*DSQRT(TH(M)/PI) / (16.d0*D_PM*D_PM)

! added Ro_g = 0 for granular flows (no gas). 
         IF(SWITCH == ZERO .OR. RO_g_avg == ZERO) THEN 
            Mu_star = eta0
         ELSEIF(TH(M) .LT. SMALL_NUMBER)THEN
            Mu_star = ZERO               
         ELSE
            Mu_star = RO_S(M)*EPS(M)*G0(M)*TH(M)*eta0 / &
               ( RO_S(M)*EPS(M)*G0(M)*TH(M) + &
               (2.d0*DgA*eta0/RO_S(M)) )     ! Note dgA is ~F_gs/ep_s
         ENDIF
   
!  shear viscosity in Mth solids phase  (add to frictional part)
         Mu_s = Mu_star * eta_star


      ELSE   ! No modifications to original mfix if 
             ! IA or GD99 theories are not used
      
!  modify F_2 if Jenkins BC is used (sof)    
         IF(JENKINS) THEN

            IF (VSLIP == ZERO) THEN
! if solids velocity field is initialized to zero, use free slip bc
               F_2 = zero

            ELSEIF(AHMADI) THEN
! Ahmadi model uses different solids pressure model
! the coefficient mu in Jenkins paper is defined as tan_Phi_w, that's how
! I understand it from soil mechanic papers, i.e., G.I. Tardos, powder
! Tech. 92 (1997), 61-74. See his equation (1). Define Phi_w in mfix.dat!
! here F_2 divided by VSLIP to use the same bc as Johnson&Jackson
                F_2 = tan_Phi_w*RO_s(M)*EPS(M)* &
                  ((ONE + 4.0D0*g0EPs_avg) + HALF*(ONE -C_e*C_e))*TH(M)/VSLIP

            ELSE
! Simonin or granular models use same solids pressure
               F_2 = tan_Phi_w*RO_s(M)*EPS(M)*(1d0+ 4.D0 * Eta *g0EPs_avg)*TH(M)/VSLIP
            ENDIF !VSLIP == ZERO

         ELSE   ! if(.not.jenkins)
 
            F_2 = (PHIP*DSQRT(3d0*TH(M))*Pi*RO_s(M)*EPS(M)*g0(M))&
               /(6d0*(ONE-ep_star_avg))

         ENDIF   ! end if(Jenkins)/else 
 
         Mu = (5d0*DSQRT(Pi*TH(M))*DP_avg(M)*RO_s(M))/96d0
         Mu_b = (256d0*Mu*EPS(M)*g0EPs_avg)/(5d0*Pi)

! This is from Wen-Yu correlation, you can put here your own single particle drag 
         Re_g = EPG*RO_g_avg*DP_avg(M)*VREL/Mu_g_avg
         IF (Re_g.lt.1000d0) THEN
            C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
         ELSE
            C_d = 0.44d0
         ENDIF
         DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(DP_avg(M)*EPG**(2.65d0))
         IF(VREL == ZERO) DgA = LARGE_NUMBER
         Beta = SWITCH*EPS(M)*DgA

! SWITCH enables us to turn on/off the modification to the
! particulate phase viscosity. If we want to simulate gas-particle
! flow then SWITCH=1 to incorporate the effect of drag on the
! particle viscosity. If we want to simulate granular flow
! without the effects of an interstitial gas, SWITCH=0.
         IF(SWITCH == ZERO .OR. Ro_g_avg == ZERO)THEN
            Mu_star = Mu
         ELSEIF(TH(M) .LT. SMALL_NUMBER)THEN
            MU_star = ZERO
         ELSE
            Mu_star = RO_S(M)*EPS(M)* g0(M)*TH(M)* Mu/ &
               (RO_S(M)*g0EPs_avg*TH(M) + 2.0d0*DgA/RO_S(M)* Mu)
         ENDIF
 
         Mu_s = ((2d0+ALPHA)/3d0)*((Mu_star/(Eta*(2d0-Eta)*&
            g0(M)))*(ONE+1.6d0*Eta*g0EPs_avg)*&
            (ONE+1.6d0*Eta*(3d0*Eta-2d0)*&
            g0EPs_avg)+(0.6d0*Mu_b*Eta))
 
! particle relaxation time
         Tau_12_st = RO_s(M)/(DgA+small_number)

         IF(SIMONIN) THEN !see calc_mu_s for explanation of these definitions
            Sigma_c = (ONE+ C_e)*(3.d0-C_e)/5.d0
            Tau_2_c = DP_avg(M)/(6.d0*EPS(M)*g0(M)*DSQRT(16.d0*(TH(M)+Small_number)/PI))
            zeta_c_2= 2.D0/5.D0*(ONE+ C_e)*(3.d0*C_e-ONE)
            Nu_t =  Tau_12_avg/Tau_12_st
            Tau_2 = ONE/(2.D0/Tau_12_st+Sigma_c/Tau_2_c)
            MU_2_T_Kin = (2.0D0/3.0D0*K_12_avg*Nu_t + TH(M) * &
                 (ONE+ zeta_c_2*EPS(M)*g0(M)))*Tau_2
            Mu_2_Col = 8.D0/5.D0*EPS(M)*g0(M)*Eta* (MU_2_T_Kin+ &
                  DP_avg(M)*DSQRT(TH(M)/PI))
            Mu_s = EPS(M)*RO_s(M)*(MU_2_T_Kin + Mu_2_Col)

         ELSEIF(AHMADI) THEN
            IF(EPS(M) < (ONE-ep_star_avg)) THEN
               Tmp_Ahmadi_Const = ONE/&
                  (ONE+ Tau_1_avg/Tau_12_st * (ONE-EPS(M)/(ONE-ep_star_avg))**3)
            ELSE
               Tmp_Ahmadi_Const = ONE
            ENDIF
            Mu_s = Tmp_Ahmadi_Const &
               *0.1045D0*(ONE/g0(M)+3.2D0*EPS(M)+12.1824D0*g0(M)*EPS(M)*EPS(M))  &
               *DP_avg(M)*RO_s(M)* DSQRT(TH(M))
         ENDIF
        
      ENDIF    ! end if for kinetic theory type
        
 
      F_HW =  F_2/Mu_s
 
      RETURN
      END FUNCTION F_HW
      

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_CALC_GRBDRY
!  Purpose: Calculate hw and cw for kinetic theory boundary conditions C
!           Cut cell version
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

      SUBROUTINE CG_CALC_GRBDRY(IJK, TYPE_OF_CELL, M, L, F_2)

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
      Use cutcell

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! IJK indices
      INTEGER, INTENT(IN) :: IJK

      CHARACTER (LEN=*) :: TYPE_OF_CELL
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Index corresponding to boundary condition
      INTEGER, INTENT(IN) ::  L      

      DOUBLE PRECISION :: F_2

!-----------------------------------------------
! Local Variables      
!-----------------------------------------------
!
!                      IJK indices 
      INTEGER          IJK1, IJK2
!
!                      Other indices
      INTEGER          IJK2E, IPJK2, IPJKM2, IPJKP2, IPJMK2, IPJPK2
      INTEGER          IJK2N, IJPK2, IJPKM2, IJPKP2, IMJPK2
      INTEGER          IJK2T, IJKP2, IJMKP2, IMJKP2
 
!
!                      Average scalars
      DOUBLE PRECISION EPg_avg, Mu_g_avg, RO_g_avg
      
!              Average void fraction at packing
      DOUBLE PRECISION ep_star_avg
!
!                      Average scalars modified to include all solid phases
      DOUBLE PRECISION EPs_avg(DIMENSION_M), DP_avg(DIMENSION_M),&
                       TH_avg(DIMENSION_M)
!
!                      Average Simonin and Ahmadi variables (sof)
      DOUBLE PRECISION K_12_avg, Tau_12_avg, Tau_1_avg
!
!                      Average velocities
      DOUBLE PRECISION WGC1, WGC2, WGCM, VGC1, VGC2, UGC1, UGC2
!
!                      Solids phase index
      INTEGER          MM
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
!
!                      Magnitude of gas-solids relative velocity
      DOUBLE PRECISION VREL
!
!              slip velocity between wall and particles for Jenkins bc (sof)
      DOUBLE PRECISION VSLIP
!
!                      radial distribution function at contact
      DOUBLE PRECISION g0(DIMENSION_M)
!
!                      Sum of eps*G_0
      DOUBLE PRECISION g0EPs_avg 
!
!                      Error message
      CHARACTER*80     LINE
!
 
!              Average Radial distribution function
      DOUBLE PRECISION g_0AVG
!----------------------------------------------- 
! Function subroutines
!----------------------------------------------- 
      DOUBLE PRECISION F_HW
!----------------------------------------------- 
! Include statement functions
!----------------------------------------------- 
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------       

!  Note:  EP_s, MU_g, and RO_g are undefined at IJK1 (wall cell).  Hence
!         IJK2 (fluid cell) is used in averages.
!

      SELECT CASE (TYPE_OF_CELL)

         CASE('U_MOMENTUM')

            IJK2 = EAST_OF(IJK)
 
            EPg_avg = (VOL(IJK)*EP_g(IJK) + VOL(IJK2)*EP_g(IJK2))/(VOL(IJK) + VOL(IJK2))
            ep_star_avg = (VOL(IJK)*EP_star_array(IJK) + VOL(IJK2)*EP_star_array(IJK2))/(VOL(IJK) + VOL(IJK2))
            Mu_g_avg = (VOL(IJK)*Mu_g(IJK) + VOL(IJK2)*Mu_g(IJK2))/(VOL(IJK) + VOL(IJK2))
            RO_g_avg = (VOL(IJK)*RO_g(IJK) + VOL(IJK2)*RO_g(IJK2))/(VOL(IJK) + VOL(IJK2))
            g0EPs_avg = ZERO

            DO MM = 1, MMAX
               g0(MM)      = G_0AVG(IJK, IJK, 'X', I_OF(IJK), M, MM)
               EPs_avg(MM) = (VOL(IJK)*EP_s(IJK, MM) + VOL(IJK2)*EP_s(IJK2, MM))/(VOL(IJK) + VOL(IJK2))
               DP_avg(MM)  = (VOL(IJK)*D_P(IJK, MM) + VOL(IJK2)*D_P(IJK2, MM))/(VOL(IJK) + VOL(IJK2))
               g0EPs_avg   = g0EPs_avg + G_0AVG(IJK, IJK, 'X', I_OF(IJK), M, MM) &
                           * (VOL(IJK)*EP_s(IJK, MM) + VOL(IJK2)*EP_s(IJK2, MM))/(VOL(IJK) + VOL(IJK2))

!               IF(GRANULAR_ENERGY) THEN
!                  TH_avg(MM) = AVG_Y(&  ! not converted to CG
!                               AVG_X(Theta_m(IJK1,MM), Theta_m(IPJMK2,MM), I_OF(IJK1)),&
!                               AVG_X(Theta_m(IJK2,MM), Theta_m(IPJK2,MM), I_OF(IJK2)),&
!                               J_OF(IJK1))
!               ELSE
                   TH_avg(MM)  = (VOL(IJK)*Theta_m(IJK,MM) + VOL(IJK2)*Theta_m(IJK2,MM))/(VOL(IJK) + VOL(IJK2))
!               ENDIF
!
            ENDDO

	    IF(SIMONIN .OR. AHMADI) THEN  ! not converted to CG
! added for Simonin and Ahmadi model (sof)
               K_12_avg = AVG_X(K_12(IJK2), K_12(IJK2E), I_OF(IJK2))	    
               Tau_12_avg = AVG_X(Tau_12(IJK2), Tau_12(IJK2E), I_OF(IJK2))	    
               Tau_1_avg = AVG_X(Tau_1(IJK2), Tau_1(IJK2E), I_OF(IJK2))
	    ELSE
	       K_12_avg = ZERO    
               Tau_12_avg = ZERO	    
               Tau_1_avg = ZERO
	    ENDIF
 
!         Calculate velocity components at i+1/2, j+1/2, k (relative to IJK1)  ! not converted to CG
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
 
            CALL GET_CG_F2(g0, EPs_avg, EPg_avg, ep_star_avg, &
	              g0EPs_avg, TH_avg, Mu_g_avg, RO_g_avg, &
	              DP_avg, K_12_avg, Tau_12_avg, Tau_1_avg, &
                      VREL, VSLIP, M,F_2)


         CASE('V_MOMENTUM')

            IJK2 = NORTH_OF(IJK)

	    EPg_avg = (VOL(IJK)*EP_g(IJK) + VOL(IJK2)*EP_g(IJK2))/(VOL(IJK) + VOL(IJK2))
            ep_star_avg = (VOL(IJK)*EP_star_array(IJK) + VOL(IJK2)*EP_star_array(IJK2))/(VOL(IJK) + VOL(IJK2))
            Mu_g_avg = (VOL(IJK)*Mu_g(IJK) + VOL(IJK2)*Mu_g(IJK2))/(VOL(IJK) + VOL(IJK2))
            RO_g_avg = (VOL(IJK)*RO_g(IJK) + VOL(IJK2)*RO_g(IJK2))/(VOL(IJK) + VOL(IJK2))
            g0EPs_avg = ZERO
  
!
	    DO MM = 1, MMAX
               g0(MM)      = G_0AVG(IJK, IJK, 'X', I_OF(IJK), M, MM)
               EPs_avg(MM) = (VOL(IJK)*EP_s(IJK, MM) + VOL(IJK2)*EP_s(IJK2, MM))/(VOL(IJK) + VOL(IJK2))
               DP_avg(MM)  = (VOL(IJK)*D_P(IJK, MM) + VOL(IJK2)*D_P(IJK2, MM))/(VOL(IJK) + VOL(IJK2))
               g0EPs_avg   = g0EPs_avg + G_0AVG(IJK, IJK, 'X', I_OF(IJK), M, MM) &
                           * (VOL(IJK)*EP_s(IJK, MM) + VOL(IJK2)*EP_s(IJK2, MM))/(VOL(IJK) + VOL(IJK2))
!
!               IF(GRANULAR_ENERGY) THEN  ! not converted to CG
!                   TH_avg(MM) = AVG_Y(&
!                                AVG_X(Theta_m(IJK1,MM), Theta_m(IPJMK2,MM), I_OF(IJK1)),&
!                                AVG_X(Theta_m(IJK2,MM), Theta_m(IPJK2,MM), I_OF(IJK2)),&
!                                J_OF(IJK1))
!               ELSE
                   TH_avg(MM)  = (VOL(IJK)*Theta_m(IJK,MM) + VOL(IJK2)*Theta_m(IJK2,MM))/(VOL(IJK) + VOL(IJK2))

!               ENDIF
!
            ENDDO

	    IF(SIMONIN .OR. AHMADI) THEN  ! not converted to CG
! added for Simonin and Ahmadi model (sof)
               K_12_avg = AVG_X(K_12(IJK2), K_12(IJK2E), I_OF(IJK2))	    
               Tau_12_avg = AVG_X(Tau_12(IJK2), Tau_12(IJK2E), I_OF(IJK2))	    
               Tau_1_avg = AVG_X(Tau_1(IJK2), Tau_1(IJK2E), I_OF(IJK2))
	    ELSE
	       K_12_avg = ZERO    
               Tau_12_avg = ZERO	    
               Tau_1_avg = ZERO
	    ENDIF
 
!         Calculate velocity components at i+1/2, j+1/2, k (relative to IJK1)    ! not converted to CG
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


            CALL GET_CG_F2(g0, EPs_avg, EPg_avg, ep_star_avg, &
	              g0EPs_avg, TH_avg, Mu_g_avg, RO_g_avg, &
	              DP_avg, K_12_avg, Tau_12_avg, Tau_1_avg, &
                      VREL, VSLIP, M,F_2)


         CASE('W_MOMENTUM')

            IJK2 = TOP_OF(IJK)

	    EPg_avg = (VOL(IJK)*EP_g(IJK) + VOL(IJK2)*EP_g(IJK2))/(VOL(IJK) + VOL(IJK2))
            ep_star_avg = (VOL(IJK)*EP_star_array(IJK) + VOL(IJK2)*EP_star_array(IJK2))/(VOL(IJK) + VOL(IJK2))
            Mu_g_avg = (VOL(IJK)*Mu_g(IJK) + VOL(IJK2)*Mu_g(IJK2))/(VOL(IJK) + VOL(IJK2))
            RO_g_avg = (VOL(IJK)*RO_g(IJK) + VOL(IJK2)*RO_g(IJK2))/(VOL(IJK) + VOL(IJK2))
            g0EPs_avg = ZERO
  
!
	    DO MM = 1, MMAX
               g0(MM)      = G_0AVG(IJK, IJK, 'X', I_OF(IJK), M, MM)
               EPs_avg(MM) = (VOL(IJK)*EP_s(IJK, MM) + VOL(IJK2)*EP_s(IJK2, MM))/(VOL(IJK) + VOL(IJK2))
               DP_avg(MM)  = (VOL(IJK)*D_P(IJK, MM) + VOL(IJK2)*D_P(IJK2, MM))/(VOL(IJK) + VOL(IJK2))
               g0EPs_avg   = g0EPs_avg + G_0AVG(IJK, IJK, 'X', I_OF(IJK), M, MM) &
                           * (VOL(IJK)*EP_s(IJK, MM) + VOL(IJK2)*EP_s(IJK2, MM))/(VOL(IJK) + VOL(IJK2))
!
!               IF(GRANULAR_ENERGY) THEN  ! not converted to CG
!                   TH_avg(MM) = AVG_Y(&
!                                AVG_X(Theta_m(IJK1,MM), Theta_m(IPJMK2,MM), I_OF(IJK1)),&
!                                AVG_X(Theta_m(IJK2,MM), Theta_m(IPJK2,MM), I_OF(IJK2)),&
!                                J_OF(IJK1))
!               ELSE
                   TH_avg(MM)  = (VOL(IJK)*Theta_m(IJK,MM) + VOL(IJK2)*Theta_m(IJK2,MM))/(VOL(IJK) + VOL(IJK2))

!               ENDIF
!
            ENDDO

	    IF(SIMONIN .OR. AHMADI) THEN  ! not converted to CG
! added for Simonin and Ahmadi model (sof)
               K_12_avg = AVG_X(K_12(IJK2), K_12(IJK2E), I_OF(IJK2))	    
               Tau_12_avg = AVG_X(Tau_12(IJK2), Tau_12(IJK2E), I_OF(IJK2))	    
               Tau_1_avg = AVG_X(Tau_1(IJK2), Tau_1(IJK2E), I_OF(IJK2))
	    ELSE
	       K_12_avg = ZERO    
               Tau_12_avg = ZERO	    
               Tau_1_avg = ZERO
	    ENDIF
 
!         Calculate velocity components at i+1/2, j+1/2, k (relative to IJK1)    ! not converted to CG
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


            CALL GET_CG_F2(g0, EPs_avg, EPg_avg, ep_star_avg, &
	              g0EPs_avg, TH_avg, Mu_g_avg, RO_g_avg, &
	              DP_avg, K_12_avg, Tau_12_avg, Tau_1_avg, &
                      VREL, VSLIP, M,F_2)



         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: GET_INTERPOLATION_TERMS_S'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:' 
            WRITE(*,*)'U_MOMENTUM' 
            WRITE(*,*)'V_MOMENTUM' 
            WRITE(*,*)'W_MOMENTUM' 
            CALL MFIX_EXIT(myPE)
          END SELECT


      RETURN

      END SUBROUTINE CG_CALC_GRBDRY
 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_CG_F2                                            C
!  Purpose: Compute F_2 for cut cell version                           C
!                                                                      C
!  Author: K. Agrawal & A. Srivastava, Princeton Univ. Date: 24-JAN-98 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
!  Modified: Sofiane Benyahia, Fluent Inc.             Date: 02-FEB-05 C
!  Purpose: Include conductivity defined by Simonin and Ahmadi         C
!           Also included Jenkins small frictional limit               C
!                                                                      C
!  Literature/Document References: See calcmu_s.f for ref. on Simonin  C
!  and Ahmadi models; for Jenkins BC: Jenkins and Louge, Phys. fluids  C
!  9 (10), 2835. See equation (2) in the paper                         C
!                                                                      C

!  Additional Notes:
!    The current implementation of the IA (2005) kinetic theory and 
!    the GD (1999) kinetic theory do not incorporate ahmadi or simonin
!    additions nor jenkins small frictional bc mdoel

!    The granular momentum BC is written as the normal vector dot the
!    stress tensor.  Besides the gradient in velocity of phase M, the
!    stress tensor expression may contain several additional terms that 
!    would need to be accounted for when satisfying the BC.  These
!    modifications have NOT been rigorously addressed.

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GET_CG_F2(g0,EPS,EPG, ep_star_avg, &
                                     g0EPs_avg,TH,Mu_g_avg,RO_g_avg,&
                                     DP_avg,K_12_avg, Tau_12_avg, Tau_1_avg, &
                                     VREL, VSLIP, M,F_2)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE constant
      USE physprop
      USE run
      USE fldvar
      USE mpi_utility
      Use cutcell      
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------      
! Radial distribution function of solids phase M with each
! other solids phase 
      DOUBLE PRECISION, INTENT(IN) :: g0(DIMENSION_M) 
! Average solids volume fraction of each solids phase
      DOUBLE PRECISION, INTENT(IN) :: EPS(DIMENSION_M)
! Average solids and gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPG, ep_star_avg
! Sum of eps*G_0 
      DOUBLE PRECISION, INTENT(INOUT) :: g0EPs_avg 
! Average theta_m
      DOUBLE PRECISION, INTENT(INOUT) :: TH (DIMENSION_M)      
! Average gas viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mu_g_avg
! Average gas density
      DOUBLE PRECISION, INTENT(IN) :: RO_g_avg
! Average particle diameter of each solids phase
      DOUBLE PRECISION, INTENT(IN) :: DP_avg(DIMENSION_M)
! Average cross-correlation K_12 and lagrangian integral time-scale
      DOUBLE PRECISION, INTENT(IN) :: K_12_avg, Tau_12_avg, Tau_1_avg
! Magnitude of slip velocity between two phases
      DOUBLE PRECISION, INTENT(IN) :: VREL
! Slip velocity between wall and particles
      DOUBLE PRECISION, INTENT(IN) :: VSLIP
! Solids phase index
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local Variables      
!-----------------------------------------------
! Solids phase index
      INTEGER :: LL
! Coefficient of 2nd term
      DOUBLE PRECISION :: F_2
! Coefficient of 1st term
      DOUBLE PRECISION :: Mu_s
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
! Constants in Simonin or Ahmadi model
      DOUBLE PRECISION :: Sigma_c, Tau_2_c, Tau_12_st, Nu_t
      DOUBLE PRECISION :: Tau_2, zeta_c_2, MU_2_T_Kin, Mu_2_Col
      DOUBLE PRECISION :: Tmp_Ahmadi_Const
! Variables for Iddir & Arastoopour model
      DOUBLE PRECISION :: MU_sM_sum, MU_s_MM, MU_s_LM, MU_sM_ip, MU_common_term,&
                          MU_sM_LM
      DOUBLE PRECISION :: M_PM, M_PL, MPSUM, NU_PL, NU_PM, D_PM, D_PL, DPSUMo2
      DOUBLE PRECISION :: Ap_lm, Dp_lm, R1p_lm, Bp_lm
! Variables for Garzo and Dufty model
      DOUBLE PRECISION :: c_star, zeta0_star, nu_eta_star, press_star, &
                          gamma_star, eta_k_star, eta_star, eta0
!----------------------------------------------- 

! This is done here similar to bc_theta to avoid small negative values of
! Theta coming most probably from linear solver
      IF(TH(M) .LE. ZERO)THEN
        TH(M) = 1D-8
        IF (myPE.eq.PE_IO) THEN
           WRITE(*,*) &
              'Warning: Negative granular temp at wall set to 1e-8'
!          CALL WRITE_ERROR('THETA_HW_CW', LINE, 1)
        ENDIF
      ENDIF

! In F_2 and Mu a DSQRT(T) has been left out as it appears in both
! terms and thus cancels out upon dividing the former by the latter
! The above statement was not implemented because Simonin viscosity
! doesn't have a sqrt(th) directly available to use this simplification.


      IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN  

! Use original IA theory if SWITCH_IA is false
         IF(.NOT. SWITCH_IA) g0EPs_avg = EPS(M)*RO_S(M)

         D_PM = DP_avg(M)
         M_PM = (PI/6.d0)*(D_PM**3)*RO_S(M)
         NU_PM = (EPS(M)*RO_S(M))/M_PM
 
         F_2 = (PHIP*DSQRT(3.d0*TH(M)/M_PM)*PI*RO_s(M)*EPS(M)*g0(M))/&
            (6.d0*(ONE-ep_star_avg))

! This is from Wen-Yu correlation, you can put here your own single particle drag
         Re_g = EPG*RO_g_avg*D_PM*VREL/Mu_g_avg
         IF (Re_g .lt. 1000.d0) THEN
            C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
         ELSE
            C_d = 0.44d0
         ENDIF
         DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(D_PM*EPG**(2.65d0))
         IF(VREL == ZERO) DgA = LARGE_NUMBER
         Beta = EPS(M)*DgA !this is equivalent to F_gs(ijk,m)

         Mu = (5.d0/96.d0)*D_PM* RO_S(M)*DSQRT(PI*TH(M)/M_PM)

         IF(.NOT.SWITCH_IA .OR. RO_g_avg == ZERO)THEN
            Mu_star = Mu
         ELSEIF(TH(M) .LT. SMALL_NUMBER)THEN
            MU_star = ZERO
         ELSE
            Mu_star = Mu*EPS(M)*g0(M)/ (g0EPs_avg+ 2.0d0*DgA*Mu / &
               (RO_S(M)**2 *(TH(M)/M_PM)))
         ENDIF

         MU_s_MM = (Mu_star/g0(M))*(1.d0+(4.d0/5.d0)*(1.d0+C_E)*g0EPs_avg)**2
         Mu_sM_sum = ZERO

         DO LL = 1, SMAX
            D_PL = DP_avg(LL)
            M_PL = (PI/6.d0)*(D_PL**3.)*RO_S(LL)
            MPSUM = M_PM + M_PL
            DPSUMo2 = (D_PM+D_PL)/2.d0
            NU_PL = (EPS(LL)*RO_S(LL))/M_PL

            IF ( LL .eq. M) THEN
               Ap_lm = MPSUM/(2.d0)
               Dp_lm = M_PL*M_PM/(2.d0*MPSUM)
               R1p_lm = ONE/( Ap_lm**1.5 * Dp_lm**3 )
               MU_s_LM = DSQRT(PI)*(DPSUMo2**4 / (48.d0*5.d0))*&
                  g0(LL)*(M_PL*M_PM/MPSUM)*(M_PL*M_PM/&
                  MPSUM)*((M_PL*M_PM)**1.5)*NU_PM*NU_PL*&
                  (1.d0+C_E)*R1p_lm*DSQRT(TH(M))

! solids phase 'viscosity' associated with the divergence 
! of solids phase M
               MU_sM_ip = (MU_s_MM + MU_s_LM)

            ELSE
               Ap_lm = (M_PM*TH(LL)+M_PL*TH(M))/&
                  (2.d0)
               Bp_lm = (M_PM*M_PL*(TH(LL)-TH(M) ))/&
                  (2.d0*MPSUM)
               Dp_lm = (M_PL*M_PM*(M_PM*TH(M)+M_PL*TH(LL) ))/&
                  (2.d0*MPSUM*MPSUM)
               R1p_lm = ( ONE/( Ap_lm**1.5 * Dp_lm**3 ) ) + &
                  ( (9.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**4 ) )+&
                  ( (30.d0*Bp_lm**4)/( 2.d0*Ap_lm**3.5 * Dp_lm**5 ) )
               MU_common_term = DSQRT(PI)*(DPSUMo2**4 / (48.d0*5.d0))*&
                  g0(LL)*(M_PL*M_PM/MPSUM)*(M_PL*M_PM/&
                  MPSUM)*((M_PL*M_PM)**1.5)*NU_PM*NU_PL*&
                  (1.d0+C_E)*R1p_lm
               MU_sM_LM = MU_common_term*(TH(M)*TH(M)*TH(LL)*TH(LL)*TH(LL) )

! solids phase 'viscosity' associated with the divergence
! of solids phase M       
               MU_sM_ip = MU_sM_LM

            ENDIF
            MU_sM_sum = MU_sM_sum + MU_sM_ip

          ENDDO

! Find the term proportional to the gradient in velocity
! of phase M  (viscosity in the Mth solids phase)
          Mu_s = MU_sM_sum

      ELSEIF (TRIM(KT_TYPE) .EQ. 'GD_99') THEN

         D_PM = DP_avg(M)
         M_PM = (PI/6.d0)*(D_PM**3)*RO_S(M)
         NU_PM = (EPS(M)*RO_S(M))/M_PM

         F_2 = (PHIP*DSQRT(3.d0*TH(M))*PI*RO_s(M)*EPS(M)*g0(M))/&
            (6.d0*(ONE-ep_star_avg))

! This is from Wen-Yu correlation, you can put here your own single particle drag
         Re_g = EPG*RO_g_avg*D_PM*VREL/Mu_g_avg
         IF (Re_g .lt. 1000.d0) THEN
            C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
         ELSE
            C_d = 0.44d0
         ENDIF
         DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(D_PM*EPG**(2.65d0))
         IF(VREL == ZERO) DgA = LARGE_NUMBER
         Beta = EPS(M)*DgA !this is equivalent to F_gs(ijk,m)
  
! Pressure/Viscosity/Bulk Viscosity
! Note: k_boltz = M_PM
     
! Find pressure in the Mth solids phase
         press_star = 1.d0 + 2.d0*(1.d0+C_E)*EPS(M)*G0(M)
 
! find bulk and shear viscosity
         c_star = 32.0d0*(1.0d0 - C_E)*(1.d0 - 2.0d0*C_E*C_E) &
            / (81.d0 - 17.d0*C_E + 30.d0*C_E*C_E*(1.0d0-C_E))

         zeta0_star = (5.d0/12.d0)*G0(M)*(1.d0 - C_E*C_E) &
            * (1.d0 + (3.d0/32.d0)*c_star)

         nu_eta_star = G0(M)*(1.d0 - 0.25d0*(1.d0-C_E)*(1.d0-C_E)) &
            * (1.d0-(c_star/64.d0))

         gamma_star = (4.d0/5.d0)*(32.d0/PI)*EPS(M)*EPS(M) &
            * G0(M)*(1.d0+C_E) * (1.d0 - (c_star/32.d0))

         eta_k_star = (1.d0 - (2.d0/5.d0)*(1.d0+C_E)*(1.d0-3.d0*C_E) &
            * EPS(M)*G0(M) ) / (nu_eta_star - 0.5d0*zeta0_star)

         eta_star = eta_k_star*(1.d0 + (4.d0/5.d0)*EPS(M)*G0(M) &
            * (1.d0+C_E) ) + (3.d0/5.d0)*gamma_star

         eta0 = 5.0d0*M_PM*DSQRT(TH(M)/PI) / (16.d0*D_PM*D_PM)

! added Ro_g = 0 for granular flows (no gas). 
         IF(SWITCH == ZERO .OR. RO_g_avg == ZERO) THEN 
            Mu_star = eta0
         ELSEIF(TH(M) .LT. SMALL_NUMBER)THEN
            Mu_star = ZERO               
         ELSE
            Mu_star = RO_S(M)*EPS(M)*G0(M)*TH(M)*eta0 / &
               ( RO_S(M)*EPS(M)*G0(M)*TH(M) + &
               (2.d0*DgA*eta0/RO_S(M)) )     ! Note dgA is ~F_gs/ep_s
         ENDIF
   
!  shear viscosity in Mth solids phase  (add to frictional part)
         Mu_s = Mu_star * eta_star


      ELSE   ! No modifications to original mfix if 
             ! IA or GD99 theories are not used
      
!  modify F_2 if Jenkins BC is used (sof)    
         IF(JENKINS) THEN

            IF (VSLIP == ZERO) THEN
! if solids velocity field is initialized to zero, use free slip bc
               F_2 = zero

            ELSEIF(AHMADI) THEN
! Ahmadi model uses different solids pressure model
! the coefficient mu in Jenkins paper is defined as tan_Phi_w, that's how
! I understand it from soil mechanic papers, i.e., G.I. Tardos, powder
! Tech. 92 (1997), 61-74. See his equation (1). Define Phi_w in mfix.dat!
! here F_2 divided by VSLIP to use the same bc as Johnson&Jackson
               F_2 = tan_Phi_w*RO_s(M)*EPS(M)* &
                  ((ONE + 4.0D0*g0EPs_avg) + HALF*(ONE -C_e*C_e))*TH(M)/VSLIP

            ELSE
! Simonin or granular models use same solids pressure
               F_2 = tan_Phi_w*RO_s(M)*EPS(M)*(1d0+ 4.D0 * Eta *g0EPs_avg)*TH(M)/VSLIP
            ENDIF !VSLIP == ZERO

         ELSE   ! if(.not.jenkins)
 
            F_2 = (PHIP*DSQRT(3d0*TH(M))*Pi*RO_s(M)*EPS(M)*g0(M))&
               /(6d0*(ONE-ep_star_avg))

         ENDIF   ! end if(Jenkins)/else 
 
         Mu = (5d0*DSQRT(Pi*TH(M))*DP_avg(M)*RO_s(M))/96d0
         Mu_b = (256d0*Mu*EPS(M)*g0EPs_avg)/(5d0*Pi)

! This is from Wen-Yu correlation, you can put here your own single particle drag 
         Re_g = EPG*RO_g_avg*DP_avg(M)*VREL/Mu_g_avg
         IF (Re_g.lt.1000d0) THEN
            C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
         ELSE
            C_d = 0.44d0
         ENDIF
         DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(DP_avg(M)*EPG**(2.65d0))
         IF(VREL == ZERO) DgA = LARGE_NUMBER
         Beta = SWITCH*EPS(M)*DgA

! SWITCH enables us to turn on/off the modification to the
! particulate phase viscosity. If we want to simulate gas-particle
! flow then SWITCH=1 to incorporate the effect of drag on the
! particle viscosity. If we want to simulate granular flow
! without the effects of an interstitial gas, SWITCH=0.
         IF(SWITCH == ZERO .OR. Ro_g_avg == ZERO)THEN
            Mu_star = Mu
         ELSEIF(TH(M) .LT. SMALL_NUMBER)THEN
            MU_star = ZERO
         ELSE
            Mu_star = RO_S(M)*EPS(M)* g0(M)*TH(M)* Mu/ &
               (RO_S(M)*g0EPs_avg*TH(M) + 2.0d0*DgA/RO_S(M)* Mu)
         ENDIF
 
         Mu_s = ((2d0+ALPHA)/3d0)*((Mu_star/(Eta*(2d0-Eta)*&
            g0(M)))*(ONE+1.6d0*Eta*g0EPs_avg)*&
            (ONE+1.6d0*Eta*(3d0*Eta-2d0)*&
            g0EPs_avg)+(0.6d0*Mu_b*Eta))
 
! particle relaxation time
         Tau_12_st = RO_s(M)/(DgA+small_number)

         IF(SIMONIN) THEN !see calc_mu_s for explanation of these definitions
            Sigma_c = (ONE+ C_e)*(3.d0-C_e)/5.d0
            Tau_2_c = DP_avg(M)/(6.d0*EPS(M)*g0(M)*DSQRT(16.d0*(TH(M)+Small_number)/PI))
            zeta_c_2= 2.D0/5.D0*(ONE+ C_e)*(3.d0*C_e-ONE)
            Nu_t =  Tau_12_avg/Tau_12_st
            Tau_2 = ONE/(2.D0/Tau_12_st+Sigma_c/Tau_2_c)
            MU_2_T_Kin = (2.0D0/3.0D0*K_12_avg*Nu_t + TH(M) * &
                 (ONE+ zeta_c_2*EPS(M)*g0(M)))*Tau_2
            Mu_2_Col = 8.D0/5.D0*EPS(M)*g0(M)*Eta* (MU_2_T_Kin+ &
                  DP_avg(M)*DSQRT(TH(M)/PI))
            Mu_s = EPS(M)*RO_s(M)*(MU_2_T_Kin + Mu_2_Col)

         ELSEIF(AHMADI) THEN
            IF(EPS(M) < (ONE-ep_star_avg)) THEN
               Tmp_Ahmadi_Const = ONE/&
                  (ONE+ Tau_1_avg/Tau_12_st * (ONE-EPS(M)/(ONE-ep_star_avg))**3)
            ELSE
               Tmp_Ahmadi_Const = ONE
            ENDIF
            Mu_s = Tmp_Ahmadi_Const &
               *0.1045D0*(ONE/g0(M)+3.2D0*EPS(M)+12.1824D0*g0(M)*EPS(M)*EPS(M))  &
               *DP_avg(M)*RO_s(M)* DSQRT(TH(M))
         ENDIF
        
      ENDIF    ! end if for kinetic theory type
        
 
!      F_HW =  F_2/Mu_s  ! Only F_2 is actually needed 

      RETURN
      END SUBROUTINE GET_CG_F2

