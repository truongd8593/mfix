!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BC_THETA( M, A_m, B_m, IER)                            C
!                                                                      C
!  Purpose: Implementation of Johnson & Jackson boundary conditions    C
!  for the pseudo-thermal temperature.                                 C

!  Notes:
!    If BC_JJ_PS = 3, then set the wall temperature equal to the 
!    adjacent fluid cell temperature (i.e. zero flux) 
!    IF BC_JJ_PS = 1, then set the wall temperature based on 
!    Johnson and Jackson Pseudo-thermal temp B.C. 
!    IF BC_JJ_PS !=1, !=3 and >0, then set the wall temperature
!    based on Johnson and Jackson Pseudo-thermal temp B.C. without
!    the generation term

!                                                                      C
!  Author: Kapil Agrawal, Princeton University        Date: 14-MAR-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE BC_THETA(M, A_m, B_m, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------      
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
      USE compar         
      USE mpi_utility    
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Solids phase index
      INTEGER :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER :: IER
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Boundary condition
      INTEGER :: L
! Indices
      INTEGER :: I,  J, K, I1, I2, J1, J2, K1, K2, IJK,&
                 IM, JM, KM
! Granular energy coefficient
      DOUBLE PRECISION :: Gw, Hw, Cw
!----------------------------------------------- 
! Include statements functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!----------------------------------------------- 


! Setup Johnson and Jackson Pseudo-thermal temp B.C.
      DO L = 1, DIMENSION_BC
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
              DO K = K1, K2
              DO J = J1, J2
              DO I = I1, I2

                IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                IJK   = FUNIJK(I, J, K)      
                IF (FLOW_AT(IJK)) CYCLE !checks for pressure outlets
                IM    = Im1(I)
                JM    = Jm1(J)
                KM    = Km1(K)
! Setting the temperature to zero - recall the center coefficient
! and source vector are negative. The off-diagonal coefficients 
! are positive.                         
                A_m(IJK, e, M) =  ZERO
                A_m(IJK, w, M) =  ZERO
                A_m(IJK, n, M) =  ZERO
                A_m(IJK, s, M) =  ZERO
                A_m(IJK, t, M) =  ZERO
                A_m(IJK, b, M) =  ZERO
                A_m(IJK, 0, M) = -ONE
                b_m(IJK, M)    =  ZERO

! Checking if the west wall                
                IF(FLUID_AT(EAST_OF(IJK)))THEN
                   IF(EP_s(EAST_OF(IJK),M).LE.DIL_EP_s)THEN
                      A_m(IJK, e, M) = ONE
                   ELSE
                      IF (BC_JJ_PS(L).EQ.3) THEN
! Setting the wall temperature equal to the adjacent fluid cell 
! temperature (i.e. zero flux)                              
                         Gw = 1d0
                         Hw = 0d0
                         Cw = 0d0 
                      ELSE
! Setting the wall temperature based on Johnson and Jackson 
! Pseudo-thermal temp B.C.                               
                         CALL CALC_THETA_BC(IJK,EAST_OF(IJK),'E',&
                            M,L,gw,hw,cw)
! Overwriting calculated value of cw                            
                        IF (BC_JJ_PS(L).EQ.2) cw=0d0 
                      ENDIF
 
                      A_m(IJK, e, M) = -(HALF*hw - oDX_E(I)*gw)
                      A_m(IJK, 0, M) = -(HALF*hw + oDX_E(I)*gw)
 
                      IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M) = -cw
                      ELSE
                        b_m(IJK,M) = ZERO
                      ENDIF
                   ENDIF
 
                ELSEIF(FLUID_AT(WEST_OF(IJK)))THEN
                   IF(EP_s(WEST_OF(IJK),M).LE.DIL_EP_s)THEN
                      A_m(IJK, w, M) = ONE
                   ELSE
                     IF (BC_JJ_PS(L).EQ.3) THEN
                        Gw = 1d0
                        Hw = 0d0
                        Cw = 0d0
                     ELSE
                        CALL CALC_THETA_BC(IJK,WEST_OF(IJK),'W',&
                           M,L,gw,hw,cw)
                        IF (BC_JJ_PS(L).EQ.2) cw=0d0
                     ENDIF
 
                     A_m(IJK, w, M) = -(HALF*hw - oDX_E(IM)*gw)
                     A_m(IJK, 0, M) = -(HALF*hw + oDX_E(IM)*gw)
 
                     IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M) = -cw
                     ELSE
                        b_m(IJK,M) = ZERO
                     ENDIF
                   ENDIF
 
                ELSEIF(FLUID_AT(NORTH_OF(IJK)))THEN
                   IF(EP_s(NORTH_OF(IJK),M).LE.DIL_EP_s)THEN
                      A_m(IJK, n, M) = ONE
                   ELSE
                     IF (BC_JJ_PS(L).EQ.3) THEN
                        Gw = 1d0
                        Hw = 0d0
                        Cw = 0d0
                     ELSE 
                        CALL CALC_THETA_BC(IJK,NORTH_OF(IJK),'N',&
                           M,L,gw,hw,cw)
                        IF (BC_JJ_PS(L).EQ.2) cw=0d0
                     ENDIF
 
                     A_m(IJK, n, M) = -(HALF*hw - oDY_N(J)*gw)
                     A_m(IJK, 0, M) = -(HALF*hw + oDY_N(J)*gw)
 
                     IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M) = -cw
                     ELSE
                        b_m(IJK,M) = ZERO
                     ENDIF
                   ENDIF
 
                ELSEIF(FLUID_AT(SOUTH_OF(IJK)))THEN
                   IF(EP_s(SOUTH_OF(IJK),M).LE.DIL_EP_s)THEN
                      A_m(IJK, s, M) = ONE
                   ELSE
                     IF (BC_JJ_PS(L).EQ.3) THEN
                        Gw = 1d0
                        Hw = 0d0
                        Cw = 0d0
                     ELSE
                        CALL CALC_THETA_BC(IJK,SOUTH_OF(IJK),'S',&
                           M,L,gw,hw,cw)
                        IF (BC_JJ_PS(L).EQ.2) cw=0d0 
                     ENDIF
 
                     A_m(IJK, s, M) = -(HALF*hw - oDY_N(JM)*gw)
                     A_m(IJK, 0, M) = -(HALF*hw + oDY_N(JM)*gw)
 
                     IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M) = -cw
                     ELSE
                        b_m(IJK,M) = ZERO
                     ENDIF
                   ENDIF
 
                ELSEIF(FLUID_AT(TOP_OF(IJK)))THEN
                   IF(EP_s(TOP_OF(IJK),M).LE.DIL_EP_s)THEN
                      A_m(IJK, t, M) = ONE
                   ELSE
                     IF (BC_JJ_PS(L).EQ.3) THEN
                        Gw = 1d0
                        Hw = 0d0
                        Cw = 0d0
                     ELSE
                        CALL CALC_THETA_BC(IJK,TOP_OF(IJK),'T',&
                           M,L,gw,hw,cw)
                        IF (BC_JJ_PS(L).EQ.2) cw=0d0
                     ENDIF
 
                     A_m(IJK, t, M) = -(HALF*hw - oX(I)*oDZ_T(K)*gw)
                     A_m(IJK, 0, M) = -(HALF*hw + oX(I)*oDZ_T(K)*gw)
 
                     IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M) = -cw
                     ELSE
                        b_m(IJK,M) = ZERO
                     ENDIF
                   ENDIF
 
                ELSEIF(FLUID_AT(BOTTOM_OF(IJK)))THEN
                   IF(EP_s(BOTTOM_OF(IJK),M).LE.DIL_EP_s)THEN
                      A_m(IJK, b , M) = ONE
                   ELSE
                     IF (BC_JJ_PS(L).EQ.3) THEN
                        Gw = 1d0
                        Hw = 0d0
                        Cw = 0d0
                     ELSE
                        CALL CALC_THETA_BC(IJK,BOTTOM_OF(IJK),'B',&
                           M,L,gw,hw,cw)
                        IF (BC_JJ_PS(L).EQ.2) cw=0d0
                     ENDIF
 
                     A_m(IJK, b, M) = -(HALF*hw - oX(I)*oDZ_T(KM)*gw)
                     A_m(IJK, 0, M) = -(HALF*hw + oX(I)*oDZ_T(KM)*gw)
 
                     IF (BC_JJ_PS(L) .EQ. 1) THEN
                        b_m(IJK, M) = -cw
                     ELSE
                        b_m(IJK,M) = ZERO
                     ENDIF
                   ENDIF
                ENDIF

              ENDDO   ! end do over i
              ENDDO   ! end do over j
              ENDDO   ! end do over k

            ENDIF   ! end if (bc_jj_ps>0)
          ENDIF  ! end if nsw, fsw or psw
        ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc

      RETURN
      END SUBROUTINE BC_THETA
 
 
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
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_THETA_BC(IJK1,IJK2,FCELL,M,L,Gw,Hw,Cw)

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
      USE rxns            
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! IJK indices for wall cell (ijk1) and fluid cell (ijk2)
      INTEGER, INTENT(IN) :: IJK1, IJK2
! The location (e,w,n...) of fluid cell
      CHARACTER, INTENT(IN) :: FCELL
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Index corresponding to boundary condition
      INTEGER, INTENT(IN) ::  L      
! Granular energy coefficient
      DOUBLE PRECISION, INTENT(INOUT) :: Gw, Hw, Cw
! Variable specularity coefficient  
     DOUBLE PRECISION PHIP_JJ 
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! IJK indices for cells
      INTEGER ::  IJK, IJK3
! Solids phase index
      INTEGER :: MM
! Average scalars
      DOUBLE PRECISION :: EPg_avg, Mu_g_avg, RO_g_avg, smallTheta
! void fraction at packing
      DOUBLE PRECISION :: ep_star_avg
! Average scalars modified to include all solid phases
      DOUBLE PRECISION :: EPs_avg(DIMENSION_M), DP_avg(DIMENSION_M), &
                          TH_avg(DIMENSION_M)
! Average Simonin variables
      DOUBLE PRECISION :: K_12_avg, Tau_12_avg
! values of U_sm, V_sm, W_sm at appropriate place on boundary wall
      DOUBLE PRECISION :: USCM, VSCM, WSCM
! values of U_g, V_g, W_g at appropriate place on boundary wall
      DOUBLE PRECISION :: UGC, VGC, WGC
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION :: VREL
! Square of slip velocity between wall and particles
      DOUBLE PRECISION :: VSLIPSQ
! Wall velocity dot outward normal for all solids phases
      DOUBLE PRECISION :: VWDOTN (DIMENSION_M)
! Gradient in number density at the wall dot with outward normal for
! all solids phases
      DOUBLE PRECISION :: GNUWDOTN (DIMENSION_M)
! Gradient in granular temperature at the wall dot with outward normal
! for all solids phases
      DOUBLE PRECISION :: GTWDOTN (DIMENSION_M)
! Error message
      CHARACTER*80     LINE(1)
! Radial distribution function
      DOUBLE PRECISION :: g0(DIMENSION_M)
! Sum of eps*G_0
      DOUBLE PRECISION :: g0EPs_avg, G_0
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

! set location independent quantities      
      EPg_avg = EP_g(IJK2)
      ep_star_avg = EP_star_array(IJK2)
      Mu_g_avg = Mu_g(IJK2)
      RO_g_avg = RO_g(IJK2)
      g0EPs_avg = ZERO

! added for Simonin model (sof)         
      IF(SIMONIN .OR. AHMADI) THEN
         K_12_avg = K_12(IJK2)
         Tau_12_avg = Tau_12(IJK2)
      ELSE
         K_12_avg = ZERO    
         Tau_12_avg = ZERO
      ENDIF

      DO MM = 1, SMAX
         g0(MM)      = G_0(IJK2, M, MM)
         EPs_avg(MM) = EP_s(IJK2,MM)
         DP_avg(MM)  = D_P(IJK2,MM)
         g0EPs_avg   = g0EPs_avg + G_0(IJK2, M, MM)*EP_s(IJK2,MM)
      ENDDO


      IF(FCELL .EQ. 'N')THEN 
         DO MM = 1, SMAX
            TH_avg(MM) = AVG_Y(Theta_m(IJK1,MM),Theta_m(IJK2,MM),J_OF(IJK1))
            IF(TH_avg(MM) < ZERO) TH_avg(MM) = smallTheta ! for some corner cells

! added for IA (2005) theory:
! include -1 since normal vector is pointing south (-y)
! velocity at wall (i,j+1/2,k relative to ijk1) dot with the normal 
! vector at the wall
            VWDOTN(MM) = -1.d0*V_S(IJK1,MM)

! gradient in number density at wall (i,j+1/2,k relative to ijk1) dot
! with the normal vector at the wall. since nu is undefined at ijk1,
! approximate gradient with interior points: ijk2 and i,j+1,k relative
! to ijk2
            IJK3 = NORTH_OF(IJK2)
            GNUWDOTN(MM) = -1.d0*(6.d0/(PI*DP_avg(MM)))*&
                 (EP_s(IJK3,MM)-EP_s(IJK2,MM))*oDY_N(J_OF(IJK2))

! gradient in granular temperature at wall (i,j+1/2,k relative to ijk1) dot
! with the normal vector at the wall.
            GTWDOTN(MM) = -1.d0*(Theta_m(IJK3,MM)-Theta_m(IJK2,MM))*&
                 oDY_N(J_OF(IJK2))
         ENDDO

! Calculate velocity components at i, j+1/2, k (relative to IJK1)
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

! slip velocity at the wall 
         VSLIPSQ=(WSCM-BC_Ww_s(L, M))**2 + (USCM-BC_Uw_s(L, M))**2
  

      ELSEIF(FCELL .EQ. 'S')THEN 
         DO MM = 1, SMAX
            TH_avg(MM) = AVG_Y(Theta_m(IJK2, MM),Theta_m(IJK1, MM),J_OF(IJK2))
            IF(TH_avg(MM) < ZERO) TH_avg(MM) = smallTheta ! for some corner cells

! added for IA (2005) theory:
! include 1 since normal vector is pointing north (+y)
! velocity at wall (i,j+1/2,k relative to ijk2) dot with the normal 
! vector at the wall. 
            VWDOTN(MM) = 1.d0*V_S(IJK2,MM)

! gradient in number density at wall (i,j+1/2,k relative to ijk2) dot
! with the normal vector at the wall. since nu is undefined at ijk1,
! approximate gradient with interior points: ijk2 and i,j-1,k relative
! to ijk2
            IJK3 = SOUTH_OF(IJK2)
            GNUWDOTN(MM) = 1.d0*(6.d0/(PI*DP_avg(MM)))*&
                 (EP_s(IJK2,MM)-EP_s(IJK3,MM))*oDY_N(J_OF(IJK3)) 

! gradient in granular temperature at wall (i,j+1/2,k relative to ijk2)
! dot with the normal vector at the wall.
            GTWDOTN(MM) = 1.d0*(Theta_m(IJK2,MM)-Theta_m(IJK3,MM))*&
                 oDY_N(J_OF(IJK3)) 
         ENDDO

! Calculate velocity components at i, j+1/2, k (relative to IJK2)
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

! slip velocity at the wall 
         VSLIPSQ=(WSCM-BC_Ww_s(L, M))**2 + (USCM-BC_Uw_s(L, M))**2


      ELSEIF(FCELL== 'E')THEN
          DO MM = 1, SMAX
            TH_avg(MM) = AVG_X(Theta_m(IJK1,MM),Theta_m(IJK2,MM),I_OF(IJK1))
            IF(TH_avg(MM) < ZERO) TH_avg(MM) = smallTheta ! for some corner cells

! added for IA (2005) theory:
! include -1 since normal vector is pointing west (-x)
! velocity at wall (i+1/2,j,k relative to ijk1) dot with the normal 
! vector at the wall. 
            VWDOTN(MM) = -1.d0*U_S(IJK1,MM)    

! gradient in number density at wall (i+1/2,j,k relative to ijk1) dot
! with the normal vector at the wall. since nu is undefined at ijk1,
! approximate gradient with interior points: ijk2 and i,i+1,k relative
! to ijk2
            IJK3 = EAST_OF(IJK2)
            GNUWDOTN(MM) = -1.d0*(6.d0/(PI*DP_avg(MM)))*&
                 (EP_s(IJK3,MM)-EP_s(IJK2,MM))*oDX_E(I_OF(IJK2)) 

! gradient in granular temperature at wall (i+1/2,j,k relative to ijk1) 
! dot with the normal vector at the wall.
            GTWDOTN(MM) = -1.d0*(Theta_m(IJK3,MM)-Theta_m(IJK2,MM))*&
                 oDX_E(I_OF(IJK2)) 
         ENDDO

! Calculate velocity components at i+1/2, j, k (relative to IJK1)
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

! slip velocity at the wall 
         VSLIPSQ=(WSCM-BC_Ww_s(L, M))**2 + (VSCM-BC_Vw_s(L, M))**2
 
 
      ELSEIF(FCELL== 'W')THEN
          DO MM = 1, SMAX
            TH_avg(MM) = AVG_X(Theta_m(IJK2,MM),Theta_m(IJK1,MM),I_OF(IJK2))
            IF(TH_avg(MM) < ZERO) TH_avg(MM) = smallTheta ! for some corner cells

! added for IA (2005) theory:
! include 1 since normal vector is pointing west (+x)
! velocity at wall (i+1/2,j,k relative to ijk2) dot with the normal 
! vector at the wall.  
            VWDOTN(MM) = 1.d0*U_S(IJK2,MM)

! gradient in number density at wall (i+1/2,j,k relative to ijk2) dot
! with the normal vector at the wall. since nu is undefined at ijk1,
! approximate gradient with interior points: ijk2 and i,i-1,k relative
! to ijk2
            IJK3 = WEST_OF(IJK2)
            GNUWDOTN(MM) = 1.d0*(6.d0/(PI*DP_avg(MM)))*&
                 (EP_s(IJK2,MM)-EP_s(IJK3,MM))*oDX_E(I_OF(IJK3)) 

! gradient in granular temperature at wall (i+1/2,j,k relative to ijk2) 
! dot with the normal vector at the wall.
            GTWDOTN(MM) = 1.d0*(Theta_m(IJK2,MM)-Theta_m(IJK3,MM))*&
                 oDX_E(I_OF(IJK3))  
         ENDDO

! Calculate velocity components at i+1/2, j, k (relative to IJK2)
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

! slip velocity at the wall 
         VSLIPSQ=(WSCM-BC_Ww_s(L, M))**2 + (VSCM-BC_Vw_s(L, M))**2
 
 
      ELSEIF(FCELL== 'T')THEN
         DO MM = 1, SMAX
            TH_avg(MM) = AVG_Z(Theta_m(IJK1,MM),Theta_m(IJK2,MM),K_OF(IJK1))
            IF(TH_avg(MM) < ZERO) TH_avg(MM) = smallTheta ! for some corner cells

! added for IA (2005) theory:
! include -1 since normal vector is pointing in bottom dir (-z)
! velocity at wall (i,j,k+1/2 relative to ijk1) dot with the normal 
! vector at the wall.  
            VWDOTN(MM) = -1.d0*W_S(IJK1,MM)

! gradient in number density at wall (i,j,k+1/2 relative to ijk1) dot
! with the normal vector at the wall. since nu is undefined at ijk1,
! approximate gradient with interior points: ijk2 and i,j,k+1 relative
! to ijk2
            IJK3 = TOP_OF(IJK2)
            GNUWDOTN(MM) = -1.d0*(6.d0/(PI*DP_avg(MM)))*&
                 (EP_s(IJK3,MM)-EP_s(IJK2,MM))*oX(I_of(IJK2))*oDZ_T(K_OF(IJK2))     

! gradient in granular temperature at wall (i,j,k+1/2 relative to ijk1) 
! dot with the normal vector at the wall.
            GNUWDOTN(MM) = -1.d0*(Theta_m(IJK3,MM)-Theta_m(IJK2,MM))*&
                 oX(I_of(IJK2))*oDZ_T(K_OF(IJK2))  
         ENDDO

! Calculate velocity components at i, j, k+1/2 (relative to IJK1)
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

! slip velocity at the wall          
         VSLIPSQ=(VSCM-BC_Vw_s(L, M))**2 + (USCM-BC_Uw_s(L, M))**2
 
 
      ELSEIF(FCELL== 'B')THEN
         DO MM = 1, SMAX
            TH_avg(MM) = AVG_Z(Theta_m(IJK2,MM),Theta_m(IJK1,MM),K_OF(IJK2))
            IF(TH_avg(MM) < ZERO) TH_avg(MM) = smallTheta ! for some corner cells

! added for IA (2005) theory:
! include 1 since normal vector is pointing in bottom dir (+z)
! velocity at wall (i,j,k+1/2 relative to ijk2) dot with the normal 
! vector at the wall 
            VWDOTN(MM) = 1.d0*W_S(IJK2,MM)

! gradient in number density at wall (i,j,k+1/2 relative to ijk2) dot
! with the normal vector at the wall. since nu is undefined at ijk1,
! approximate gradient with interior points: ijk2 and i,j,k-1 relative
! to ijk2
            IJK3 = BOTTOM_OF(IJK2)
            GNUWDOTN(MM) = 1.d0*(6.d0/(PI*DP_avg(MM)))*&
                 (EP_s(IJK2,MM)-EP_s(IJK3,MM))*oX(I_of(IJK3))*oDZ_T(K_OF(IJK3))     

! gradient in granular temperature at wall (i,j,k+1/2 relative to ijk2) 
! dot with the normal vector at the wall.
            GTWDOTN(MM) = 1.d0*(Theta_m(IJK2,MM)-Theta_m(IJK3,MM))*&
                 oX(I_of(IJK3))*oDZ_T(K_OF(IJK3))     
         ENDDO

! Calculate velocity components at i, j, k+1/2 (relative to IJK2)
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

! slip velocity at the wall 
         VSLIPSQ=(VSCM-BC_Vw_s(L, M))**2 + (USCM-BC_Uw_s(L, M))**2

      ELSE
         WRITE(LINE,'(A, A)') 'Error: Unknown FCELL'
         CALL WRITE_ERROR('CALC_THETA_BC', LINE, 1)
         call exitMPI(myPE)          
      ENDIF

! magnitude of gas-solids relative velocity
      VREL = DSQRT( (UGC-USCM)**2 + (VGC-VSCM)**2 + (WGC-WSCM)**2 )

      CALL THETA_Hw_Cw(g0, EPs_avg, EPg_avg, ep_star_avg, &
                       g0EPs_avg, TH_avg,Mu_g_avg,RO_g_avg, &
                       DP_avg, K_12_avg,Tau_12_avg,VREL,VSLIPSQ,VWDOTN,&
                       GNUWDOTN,GTWDOTN,M,Gw,Hw,Cw,L)
!
!     Output the variable specularity coefficient 
!     This only work for one solid phase
      if(BC_JJ_M .and. PHIP_OUT_JJ)then
      	if(PHIP_OUT_ITER.eq.1)then
      		ReactionRates(IJK1,1)= PHIP_JJ(dsqrt(VSLIPSQ),TH_avg(1))
      	endif
      endif
      RETURN
      END SUBROUTINE CALC_THETA_BC
 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SUBROUTINE THETA_HW_CW
!  Purpose: Subroutine for gw, hw and cw                               C
!                                                                      C
!  Author: Kapil Agrawal, Princeton University         Date: 15-MAR-98 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
!                                                                      C
!  Modified: Sofiane Benyahia, Fluent Inc.             Date: 02-FEB-05 C
!  Purpose: Include conductivity defined by Simonin and Ahmadi         C
!           Also included Jenkins small frictional limit               C
!                                                                      C
!  Literature/Document References: See calcmu_s.f for ref. on Simonin  C
!  and Ahmadi models; for Jenkins BC: Jenkins and Louge, Phys. fluids  C
!  9 (10), 2835. See equation (3) in the paper                         C
!           
!  Additional Notes:
!    The current implementation of the IA (2005) kinetic theory and 
!    the GD (1999) kinetic theory do not incorporate ahmadi or simonin
!    additions nor jenkins small frictional bc mdoel

!    The granular energy BC is written as the normal vector dot the
!    heat flux.  In granular kinetic theory the expression for heat 
!    flux generally contains terms beyond the gradient in temperature.
!    Thus to rigorously implement the granualr energy BC would require
!    accounting for the additional terms.  Such modifications were
!    not done for the kinetic theories here.
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE THETA_HW_CW(g0,EPS,EPG, ep_star_avg, &
                             g0EPs_avg,TH,Mu_g_avg,RO_g_avg, DP_avg, &
                             K_12_avg,Tau_12_avg,VREL,VSLIPSQ,VWDOTN,&
                             GNUWDOTN,GTWDOTN,M,GW,HW,CW,L)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE run
      USE constant
      USE fldvar
      USE toleranc 
      USE bc
      USE compar
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
      DOUBLE PRECISION, INTENT(IN) :: K_12_avg, Tau_12_avg
! Magnitude of slip velocity between two phases
      DOUBLE PRECISION, INTENT(IN) :: VREL
! Square of slip velocity between wall and particles
      DOUBLE PRECISION, INTENT(IN) :: VSLIPSQ
! Wall velocity dot outward normal for all solids phases
      DOUBLE PRECISION, INTENT(IN) :: VWDOTN (DIMENSION_M)
! Gradient in number density at the wall dot with outward
! normal for all solids phases
      DOUBLE PRECISION, INTENT(IN) :: GNUWDOTN (DIMENSION_M)
! Gradient in temperature at the wall dot with outward
! normal for all solids phases
      DOUBLE PRECISION, INTENT(IN) :: GTWDOTN (DIMENSION_M)
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Coefficients in boundary conditions
      DOUBLE PRECISION, INTENT(INOUT) :: GW, HW, CW 
! Index corresponding to boundary condition
      INTEGER, INTENT(IN) :: L
!-----------------------------------------------
! Local Variables      
!-----------------------------------------------
! Solids phase index
      INTEGER :: LL
! Coefficient of 1st term
      DOUBLE PRECISION :: K_1
!      DOUBLE PRECISION :: K_10
! Conductivity
      DOUBLE PRECISION :: Kgran
!      DOUBLE PRECISION :: Lambda0
! Conductivity corrected for interstitial fluid effects
      DOUBLE PRECISION :: Kgran_star
! Reynolds number based on slip velocity
      DOUBLE PRECISION :: Re_g
! Friction Factor in drag coefficient
      DOUBLE PRECISION :: C_d
! Drag Coefficient
      DOUBLE PRECISION :: Beta, DgA
! Local parameters for Simonin model
      DOUBLE PRECISION :: Zeta_c, Omega_c, Tau_2_c, Kappa_kin, &
                          Kappa_Col, Tau_12_st
! Variables for Iddir & Arastoopour model
      DOUBLE PRECISION :: M_PM, M_PL, MPSUM, NU_PL, NU_PM, D_PM, D_PL, DPSUMo2
      DOUBLE PRECISION :: Ap_lm, Dp_lm, R0p_lm, R1p_lm, R8p_lm, R9p_lm, &
                          Bp_lm, R5p_lm, R6p_lm, R7p_lm
      DOUBLE PRECISION :: K_s_sum, K_s_MM, K_s_LM, K_sM_LM
      DOUBLE PRECISION :: Kvel_s_sum, Kvel_s_LM, Kth_sL_LM, Kth_s_sum
      DOUBLE PRECISION :: Knu_s_sum, Knu_s_LM, K_common_term, K
! Variables for Garzo & Dufty model
      DOUBLE PRECISION :: c_star, zeta0_star, press_star, eta0, &
                          kappa0, nu_kappa_star, kappa_k_star, &
                          kappa_star
! Error message
      CHARACTER*80     LINE
!----------------------------------------------- 
! Function subroutines
!----------------------------------------------- 
      DOUBLE PRECISION PHIP_JJ 
!-----------------------------------------------
 
      IF(TH(M) .LE. ZERO)THEN
         TH(M) = 1D-8
         IF (myPE.eq.PE_IO) THEN
            WRITE(*,*)  &
               'Warning: Negative granular temp at wall set to 1e-8'
!            CALL WRITE_ERROR('THETA_HW_CW', LINE, 1)
         ENDIF
      ENDIF

      IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN

! Use original IA theory if SWITCH_IA is false
         IF(.NOT. SWITCH_IA) g0EPs_avg = EPS(M)*RO_S(M)

         D_PM = DP_avg(M)        
         M_PM = (PI/6.d0)*(D_PM**3.)*RO_S(M)
         NU_PM = (EPS(M)*RO_S(M))/M_PM

         K_s_sum = ZERO
         Kvel_s_sum = ZERO
         Knu_s_sum = ZERO
         Kth_s_sum = ZERO

         Kgran = (75.d0/384.d0)*RO_s(M)*D_PM*DSQRT(Pi*TH(M)/M_PM)

! This is from Wen-Yu correlation, you can put here your own single particle drag 
         Re_g = EPG*RO_g_avg*D_PM*VREL/Mu_g_avg
         IF (Re_g .lt. 1000.d0) THEN
            C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
         ELSE
            C_d = 0.44d0
         ENDIF
         DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(DP_avg(M)*EPG**(2.65d0))
         IF(VREL == ZERO) DgA = LARGE_NUMBER
         Beta = EPS(M)*DgA
 
         IF(.NOT.SWITCH_IA .OR. RO_g_avg == ZERO)THEN
            Kgran_star = Kgran
         ELSEIF(TH(M)/M_PM .LT. SMALL_NUMBER)THEN
            Kgran_star = ZERO
         ELSE
            Kgran_star = Kgran*g0(M)*EPS(M)/ &
               (g0EPs_avg+ 1.2d0*DgA*Kgran / (RO_S(M)**2 *(TH(M)/M_PM)))
         ENDIF

         K_s_MM = (Kgran_star/(M_PM*g0(M)))*&  ! Kth doesn't include the mass.
            (1.d0+(3.d0/5.d0)*(1.d0+C_E)*(1.d0+C_E)*g0EPs_avg)**2

         DO LL = 1, SMAX
            D_PL = DP_avg(LL)
            M_PL = (PI/6.d0)*(D_PL**3.)*RO_S(LL)
            MPSUM = M_PM + M_PL
            DPSUMo2 = (D_PM+D_PL)/2.d0
            NU_PL = (EPS(LL)*RO_S(LL))/M_PL

            IF ( LL .eq. M) THEN
               K_s_sum = K_s_sum + K_s_MM
! Kth_sL_LM is zero when LL=M because it cancels with terms from K_s_LM 
! Kvel_s_LM should be zero when LL=M (same solids phase)
! Knu_s_LM should be zero when LL=M (same solids phase)
            ELSE
               Ap_lm = (M_PM*TH(LL)+M_PL*TH(M))/(2.d0)
               Bp_lm = (M_PM*M_PL*(TH(LL)-TH(M) ))/(2.d0*MPSUM)
               Dp_lm = (M_PL*M_PM*(M_PM*TH(M)+M_PL*TH(LL) ))/&
                  (2.d0*MPSUM*MPSUM)
               R0p_lm = ( 1.d0/( Ap_lm**1.5 * Dp_lm**2.5 ) )+ &
                  ( (15.d0*Bp_lm*Bp_lm)/( 2.d0* Ap_lm**2.5 * Dp_lm**3.5 ) )+&
                  ( (175.d0*(Bp_lm**4))/( 8.d0*Ap_lm**3.5 * Dp_lm**4.5 ) )
               R1p_lm = ( 1.d0/( (Ap_lm**1.5)*(Dp_lm**3) ) )+ &
                  ( (9.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**4 ) )+&
                  ( (30.d0*Bp_lm**4) /( 2.d0*Ap_lm**3.5 * Dp_lm**5 ) )
               R5p_lm = ( 1.d0/( Ap_lm**2.5 * Dp_lm**3 ) )+ &
                  ( (5.d0*Bp_lm*Bp_lm)/( Ap_lm**3.5 * Dp_lm**4 ) )+&
                  ( (14.d0*Bp_lm**4)/( Ap_lm**4.5 * Dp_lm**5 ) )
               R6p_lm = ( 1.d0/( Ap_lm**3.5 * Dp_lm**3 ) )+ &
                  ( (7.d0*Bp_lm*Bp_lm)/( Ap_lm**4.5 * Dp_lm**4 ) )+&
                  ( (126.d0*Bp_lm**4)/( 5.d0*Ap_lm**5.5 * Dp_lm**5 ) )
               R7p_lm = ( 3.d0/( 2.d0*Ap_lm**2.5 * Dp_lm**4 ) )+ &
                  ( (10.d0*Bp_lm*Bp_lm)/( Ap_lm**3.5 * Dp_lm**5 ) )+&
                  ( (35.d0*Bp_lm**4)/( Ap_lm**4.5 * Dp_lm**6 ) )
               R8p_lm = ( 1.d0/( 2.d0*Ap_lm**1.5 * Dp_lm**4 ) )+ &
                  ( (6.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**5 ) )+&
                  ( (25.d0*Bp_lm**4)/( Ap_lm**3.5 * Dp_lm**6 ) )
               R9p_lm = ( 1.d0/( Ap_lm**2.5 * Dp_lm**3 ) )+ &
                  ( (15.d0*Bp_lm*Bp_lm)/( Ap_lm**3.5 * Dp_lm**4 ) )+&
                  ( (70.d0*Bp_lm**4)/( Ap_lm**4.5 * Dp_lm**5 ) )
               K_common_term = DPSUMo2**3 * M_PL*M_PM/(2.d0*MPSUM)*&
                  (1.d0+C_E)*g0(LL) * (M_PM*M_PL)**1.5

! solids phase 'conductivity' associated with the 
! gradient in granular temperature of species M
               K_sM_LM = - K_common_term*NU_PM*NU_PL*(&
                  ((DPSUMo2*DSQRT(PI)/16.d0)*(3.d0/2.d0)*Bp_lm*R5p_lm)+&
                  ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*PI/6.d0)*&
                  (3.d0/2.d0)*R1p_lm)-(&
                  ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PM/8.d0)*Bp_lm*R6p_lm)+&
                  ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
                  8.d0)*M_PM*R9p_lm)+&
                  ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PL*M_PM/(MPSUM*MPSUM))*&
                  M_PL*Bp_lm*R7p_lm)+&
                  ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
                  2.d0)*(M_PM/MPSUM)**2 * M_PL*R8p_lm)+&
                  ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PM*M_PL/(2.d0*MPSUM))*&
                  R9p_lm)-&
                  ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*DPSUMo2*DSQRT(PI)*&
                  (M_PM*M_PL/MPSUM)*Bp_lm*R7p_lm) )*TH(LL) )*&
                  (TH(M)**2 * TH(LL)**3)

! These lines were commented because they are not currently used (sof)
! solids phase 'conductivity' associated with the gradient in granular
! temperature of species L 
!               Kth_sL_LM = K_common_term*NU_PM*NU_PL*(&
!                  (-(DPSUMo2*DSQRT(PI)/16.d0)*(3.d0/2.d0)*Bp_lm*R5p_lm)+&
!                  (-(M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*PI/6.d0)*&
!                  (3.d0/2.d0)*R1p_lm)+(&
!                  ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PL/8.d0)*Bp_lm*R6p_lm)+&
!                  ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
!                  8.d0)*M_PL*R9p_lm)+&
!                  ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PL*M_PM/(MPSUM*MPSUM))*&
!                  M_PM*Bp_lm*R7p_lm)+&
!                  ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
!                  2.d0)*(M_PM*M_PM/(MPSUM*MPSUM))*M_PM*R8p_lm)+&
!                  ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PM*M_PL/(2.d0*MPSUM))*&
!                  R9p_lm)+&
!                  ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*DPSUMo2*DSQRT(PI)*&
!                  (M_PM*M_PL/MPSUM)*Bp_lm*R7p_lm) )*TH(M) )*&
!                  (TH(LL)*TH(LL)*(TH(M)**3.))*(GTWDOTN(LL))

! solids phase 'conductivity' associated with the difference in velocities
!               Kvel_s_LM = K_common_term*NU_PM*NU_PL*&
!                  (M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(3.d0*PI/10.d0)*&
!                  R0p_lm*( (TH(M)*TH(LL))**(5./2.) )*&
!                  (VWDOTN(M)-VWDOTN(LL))

! solids phase 'conductivity' associated with the difference in the 
! gradient in number densities
!               Knu_s_LM = K_common_term*(((DPSUMo2*DSQRT(PI)/16.d0)*&
!                  Bp_lm*R5p_lm)+((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*&
!                  (DPSUMo2*PI/6.d0)*R1p_lm) )*( (TH(M)*TH(LL))**(3.) )*&
!                  (NU_PM*GNUWDOTN(LL)-NU_PL*GNUWDOTN(M))
!

               K_s_sum = K_s_sum + K_sM_LM
!               Kth_s_sum = Kth_s_sum + Kth_sL_LM
!               Kvel_s_sum = Kvel_s_sum + Kvel_s_LM
!               Knu_s_sum = Knu_s_sum + Knu_s_LM

            ENDIF   ! end if( LL .eq. M)/else
         ENDDO  ! end do LL = 1, SMAX
               
         K_1 = K_s_sum

         
! setting the coefficients for JJ BC
         GW = M_PM * K_1 !sof modified

! Note that IA (2005) theory defines theta in terms of mass so
! the boundary conditions must be adjusted for this definition
! of granular temperature (JJ BC do not have mass in definition
! of granular temperature)
         HW = (PI*DSQRT(3.d0)/(4.d0*(ONE-ep_star_avg)))*(1.d0-e_w*e_w)*&
            RO_s(M)*EPS(M)*g0(M)*DSQRT(TH(M)/M_PM)
          if(.NOT. BC_JJ_M)then                
          CW = (PI*DSQRT(3.d0)/(6.d0*(ONE-ep_star_avg)))*PHIP*RO_s(M)*&
               EPS(M)*g0(M)*DSQRT(TH(M)*M_PM)*VSLIPSQ
          else
          CW = (PI*DSQRT(3.d0)/(6.d0*(ONE-ep_star_avg)))*PHIP_JJ(dsqrt(vslipsq),th(m))*RO_s(M)*&
               EPS(M)*g0(M)*DSQRT(TH(M)*M_PM)*VSLIPSQ
          endif

! Note that the velocity term is not included here because it should
! become zero when dotted with the outward normal (i.e. no solids 
! flux through the wall; assuming that the solids velocity at the
! wall in the normal direction is zero)
!         CW = CW + Kvel_s_sum

! Note that the gradient in number density at the wall must be
! approximated with interior points since there is no value of
! number density associated with the wall location; the ghost
! cell values are undefined for volume fraction.  
!         CW = CW + Knu_s_sum

! The gradient in temperature of phase LL at the wall
!         CW = CW + Kth_s_sum


      ELSEIF (TRIM(KT_TYPE) .EQ. 'GD_99') THEN
         D_PM = DP_avg(M)        
         M_PM = (PI/6.d0)*(D_PM**3.)*RO_S(M)
         NU_PM = (EPS(M)*RO_S(M))/M_PM

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

! Conductivity
! Note: k_boltz = M_PM

! Find pressure in the Mth solids phase
         press_star = 1.d0 + 2.d0*(1.d0+C_E)*EPS(M)*G0(M)
 
! find conductivity  
         eta0 = 5.0d0*M_PM*DSQRT(TH(M)/PI) / (16.d0*D_PM*D_PM)
 
         c_star = 32.0d0*(1.0d0 - C_E)*(1.d0 - 2.0d0*C_E*C_E) &
            / (81.d0 - 17.d0*C_E + 30.d0*C_E*C_E*(1.0d0-C_E))

         zeta0_star = (5.d0/12.d0)*G0(M)*(1.d0 - C_E*C_E) &
            * (1.d0 + (3.d0/32.d0)*c_star)

         kappa0 = (15.d0/4.d0)*eta0

         nu_kappa_star = (G0(M)/3.d0)*(1.d0+C_E) * ( 1.d0 + &
            (33.d0/16.d0)*(1.d0-C_E) + ((19.d0-3.d0*C_E)/1024.d0)*c_star)
!         nu_mu_star = nu_kappa_star

         kappa_k_star = (2.d0/3.d0)*(1.d0 + 0.5d0*(1.d0+press_star)*c_star + &
            (3.d0/5.d0)*EPS(M)*G0(M)*(1.d0+C_E)*(1.d0+C_E) * &
            (2.d0*C_E - 1.d0 + ( 0.5d0*(1.d0+C_E) - 5.d0/(3*(1.d0+C_E))) * &
            c_star ) ) / (nu_kappa_star - 2.d0*zeta0_star)

         kappa_star = kappa_k_star * (1.d0 + (6.d0/5.d0)*EPS(M)* &
            G0(M)*(1.d0+C_E) ) + (256.d0/25.d0)*(EPS(M)* &
            EPS(M)/PI)*G0(M)*(1.d0+C_E)*(1.d0+(7.d0/32.d0)* &
            c_star)

         IF(SWITCH == ZERO .OR. Ro_g_avg == ZERO)THEN
            Kgran_star = kappa0
         ELSEIF(TH(M) .LT. SMALL_NUMBER)THEN
            Kgran_star = ZERO
         ELSE
            Kgran_star = RO_S(M)*EPS(M)* G0(M)*TH(M)* kappa0/ &
               (RO_S(M)*G0(M)*EPS(M)*TH(M) + &
               1.2d0*DgA*kappa0/RO_S(M))     ! Note dgA is ~F_gs/ep_s
         ENDIF
  
! setting the coefficients for JJ BC         
! granular conductivity in Mth solids phase
         GW = Kgran_star * kappa_star

         HW = (Pi*DSQRT(3d0)/(4.D0*(ONE-ep_star_avg)))*(1d0-e_w*e_w)*&
            RO_s(M)*EPS(M)*g0(M)*DSQRT(TH(M))
 
          if(.NOT. BC_JJ_M)then 
          CW = (Pi*DSQRT(3d0)/(6.D0*(ONE-ep_star_avg)))*PHIP*RO_s(M)*&
               EPS(M)*g0(M)*DSQRT(TH(M))*VSLIPSQ
          else
          CW = (Pi*DSQRT(3d0)/(6.D0*(ONE-ep_star_avg)))*PHIP_JJ(dsqrt(vslipsq),th(m))*RO_s(M)*&
               EPS(M)*g0(M)*DSQRT(TH(M))*VSLIPSQ
	  endif

! transport coefficient of the Mth solids phase associated
! with gradient in volume fraction in heat flux
!         qmu_k_star = 2.d0*( (1.d0+EPS(M)*DG_0DNU(EPS(M)))* &
!            zeta0_star*kappa_k_star + ( (press_star/3.d0) + (2.d0/3.d0)* &
!            EPS(M)*(1.d0+C_E) * (G0(M)+EPS(M)* &
!            DG_0DNU(EPS(M))) )*c_star - (4.d0/5.d0)*EPS(M)* &
!            G0(M)* (1.d0+(EPS(M)/2.d0)*DG_0DNU(EPS(M)))* &
!            (1.d0+C_E) * ( C_E*(1.d0-C_E)+0.25d0*((4.d0/3.d0)+C_E* &
!            (1.d0-C_E))*c_star ) ) / (2.d0*nu_kappa_star-3.d0*zeta0_star)
!         qmu_star = qmu_k_star*(1.d0+(6.d0/5.d0)*EPS(M)*G0(M)*&
!            (1.d0+C_E) )
!         Kphis = (TH(M)*Kgran_star/NU_PM)*qmu_star


      ELSE   ! No modifications to original mfix if 
             ! IA or GD99 theories are not used
 
         Kgran = 75d0*RO_s(M)*DP_avg(M)*DSQRT(Pi*TH(M))/(48*Eta*(41d0-33d0*Eta))
 
         Re_g = EPG*RO_g_avg*DP_avg(M)*VREL/Mu_g_avg
         IF (Re_g.lt.1000d0) THEN
            C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
         ELSE
            C_d = 0.44d0
         ENDIF
         DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(DP_avg(M)*EPG**(2.65d0))
         IF(VREL == ZERO) DgA = LARGE_NUMBER
         Beta = SWITCH*EPS(M)*DgA

! particle relaxation time
         Tau_12_st = RO_s(M)/(DgA+small_number)
 
! SWITCH enables us to turn on/off the modification to the
! particulate phase viscosity. If we want to simulate gas-particle
! flow then SWITCH=1 to incorporate the effect of drag on the
! particle viscosity. If we want to simulate granular flow
! without the effects of an interstitial gas, SWITCH=0.
         IF(SWITCH == ZERO .OR. Ro_g_avg == ZERO)THEN
            Kgran_star = Kgran
         ELSEIF(TH(M) .LT. SMALL_NUMBER)THEN
            Kgran_star = ZERO
         ELSE
            Kgran_star = RO_S(M)*EPS(M)* g0(M)*TH(M)* Kgran/ &
               (RO_S(M)*g0EPs_avg*TH(M) + &
               1.2d0*SWITCH*DgA/RO_S(M)* Kgran)
         ENDIF
 
         K_1 = Kgran_star/g0(M)*(&
            ( ONE + (12d0/5.d0)*Eta*g0EPs_avg )&
            * ( ONE + (12d0/5.d0)*Eta*Eta*(4d0*Eta-3d0)&
            *g0EPs_avg )&
            + (64d0/(25d0*Pi)) * (41d0-33d0*Eta) *&
            (Eta*g0EPs_avg)**2 &
            )
 
         IF(SIMONIN) THEN
            Zeta_c  = (ONE+ C_e)*(49.d0-33.d0*C_e)/100.d0

            Omega_c = 3.d0*(ONE+ C_e)**2 *(2.d0*C_e-ONE)/5.d0 

            Tau_2_c = DP_avg(M)/(6.d0*EPS(M)*g0(M) &
               *DSQRT(16.d0*(TH(M)+Small_number)/PI))

! Defining Simonin's Solids Turbulent Kinetic diffusivity: Kappa
            Kappa_kin = (9.d0/10.d0*K_12_avg*(Tau_12_avg/Tau_12_st) &
               + 3.0D0/2.0D0 * TH(M)*(ONE+ Omega_c*EPS(M)*g0(M)))/     &
               (9.d0/(5.d0*Tau_12_st) + zeta_c/Tau_2_c)

            Kappa_Col = 18.d0/5.d0*EPS(M)*g0(M)*Eta* (Kappa_kin+ &
               5.d0/9.d0*DP_avg(M)*DSQRT(TH(M)/PI))

            K_1 =  EPS(M)*RO_s(M)*(Kappa_kin + Kappa_Col)
 
         ELSEIF(AHMADI) THEN
            K_1 =  0.1306D0*RO_s(M)*DP_avg(M)*(ONE+C_e**2)* (  &
               ONE/g0(M)+4.8D0*EPS(M)+12.1184D0 *EPS(M)*EPS(M)*g0(M) )*&
               DSQRT(TH(M))

         ENDIF !for simonin or ahmadi models
      
! setting the coefficients for JJ BC
         GW = K_1
      
! modify HW and CW if Jenkins BC is used (sof)    
         IF(JENKINS) THEN

            IF(AHMADI) THEN
! Ahmadi model uses different solids pressure model
               HW = 3.D0/8.D0*DSQRT(3.D0*TH(M))*((1d0-e_w))*&
                  RO_s(M)*EPS(M)*((ONE + 4.0D0*g0EPs_avg) +&
                  HALF*(ONE -C_e*C_e))

! the coefficient mu in Jenkins paper is defined as tan_Phi_w, that's how
! I understand it from soil mechanic papers, i.e., G.I. Tardos, powder
! Tech. 92 (1997), 61-74. See his equation (1). Define Phi_w in mfix.dat!
               CW = tan_Phi_w*tan_Phi_w*(ONE+e_w)*21.d0/16.d0*&
                  DSQRT(3.D0*TH(M)) *RO_s(M)*EPS(M)*&
                  ((ONE + 4.0D0*g0EPs_avg) + HALF*(ONE -C_e*C_e))*TH(M)

            ELSE  
! Simonin or granular models use same solids pressure
               HW = 3.D0/8.D0*DSQRT(3.*TH(M))*((1d0-e_w))*&
                  RO_s(M)*EPS(M)*(1d0+ 4.D0*Eta*g0EPs_avg)
               CW = tan_Phi_w*tan_Phi_w*(ONE+e_w)*21.D0/16.D0*&
                  DSQRT(3.D0*TH(M))*RO_s(M)*EPS(M)*&
                  (1d0+ 4.D0*Eta*g0EPs_avg)*TH(M)
            ENDIF   ! end if for Ahmadi
 
         ELSE   ! if(.not.jenkins) branch

            HW = (Pi*DSQRT(3d0)/(4.D0*(ONE-ep_star_avg)))*(1d0-e_w*e_w)*&
               RO_s(M)*EPS(M)*g0(M)*DSQRT(TH(M))
               if(.NOT. BC_JJ_M)then 
               CW = (Pi*DSQRT(3d0)/(6.D0*(ONE-ep_star_avg)))*PHIP*RO_s(M)*&
                    EPS(M)*g0(M)*DSQRT(TH(M))*VSLIPSQ
               else
               CW = (Pi*DSQRT(3d0)/(6.D0*(ONE-ep_star_avg)))*PHIP_JJ(dsqrt(vslipsq),th(m))*RO_s(M)*&
                    EPS(M)*g0(M)*DSQRT(TH(M))*VSLIPSQ  
               endif	

         ENDIF   ! end if(Jenkins)/else

      ENDIF   ! end if for kinetic theory type
 
 
      RETURN
      END SUBROUTINE THETA_HW_CW

!----------------------------------------------------------------------------
!Module name phip_jj
!calculate the specularity coefficient for JJ boundary conditions
!Author: Tingwen Li
!Reference:
!T. Li and S. Benyahia, Revisiting Johnson and Jackson boundary conditions 
!for granular flows, AIChE Journal, DOI: 10.1002/aic.12728
!----------------------------------------------------------------------------

	double precision function phip_jj(uslip,g_theta)
	use constant  !e_w, phi_w, PI
	implicit none
	double precision uslip,g_theta,r4phi
	
!	k4phi=7.d0/2.d0*mu4phi*(1.d0+e_w)
	r4phi=uslip/dsqrt(3.d0*g_theta)
	
	if(r4phi .le. 4.d0*k4phi/7.d0/dsqrt(6.d0*PI)/phip0)then
		phip_jj=-7.d0*dsqrt(6.d0*PI)*phip0**2/8.d0/k4phi*r4phi+phip0
        else
        	phip_jj=2.d0/7.d0*k4phi/r4phi/dsqrt(6.d0*PI)
	endif

	return      	
	end
	      
