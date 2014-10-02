!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_OUTFLOW                                             C
!  Purpose: Set specified pressure outflow bc for a specified range    C
!           of cells. This routine is also called for mass_outlow or   C
!           outflow bcs                                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables modified:                                                 C
!     T_g, X_g, T_s, X_s, Scalar, Theta_m,                             C
!     K_Turb_g & E_Turb_g (if K_epsilon),                              C
!     P_star, P_s, MW_MIX_g, P_g (if not PO)                           C
!     RO_g (if ro_g0 == undefined), ROP_s, EP_g, ROP_g,                C
!     Flux_gE, Flux_gN, Flux_gT, Flux_sE, Flux_sN, Flux_sT;            C
!     Flux_gE, Flux_gN, Flux_gT Flux_sE, Flux_sN, &                    C
!     Flux_sT (if ADDED_MASS); Flux_nE, Flux_nN, Flux_nT (if GHD);     C
!     U_g, V_g, W_g, U_s, V_s, W_s                                     C
!                                                                      C
!  Comments:                                                           C
!     If the outflow boundary is on the W, S or B side of the domain   C
!     and the component of velocity through the plane is defined then  C
!     this routine will NOT modify it (i.e., its value from the        C
!     momentum solver is maintained).                                  C
!                                                                      C
!  Local variables: IJK, I, J, K, M, FIJK                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_OUTFLOW(BCV, I1, I2, J1, J2, K1, K2)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE bc
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE scalars
      USE run
      USE compar
      USE mflux
      USE discretelement
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Boundary condition number
      INTEGER, INTENT(IN) :: BCV
! Starting and ending I index
      INTEGER, INTENT(IN) :: I1, I2
! Starting and ending J index
      INTEGER, INTENT(IN) :: J1, J2
! Starting and ending K index
      INTEGER, INTENT(IN) :: K1, K2
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K, M
! index for boundary cell
      INTEGER :: IJK
! index for a fluid cell adjacent to the boundary cell
      INTEGER :: FIJK
! alias for gas and solids velocity
      DOUBLE PRECISION :: RVEL_G, RVEL_S(DIMENSION_M)
!-----------------------------------------------

      DO K = K1, K2
         DO J = J1, J2
            DO I = I1, I2
! Check if current i,j,k resides on this PE
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF(DEAD_CELL_AT(I,J,K)) CYCLE
               IJK = FUNIJK(I,J,K)

! Fluid cell at West
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(IM_OF(IJK))) THEN
                  FIJK = IM_OF(IJK)
                  RVEL_G = U_G(FIJK)
                  DO M = 1,MMAX
                     RVEL_S(M) = U_S(FIJK,M)
                  ENDDO

                  CALL SET_OUTFLOW2(BCV, IJK, FIJK, RVEL_G, RVEL_S)

! set boundary cell values of velocity according to adjacent fluid cell
! values
                  IF (ROP_G(IJK) > ZERO) THEN
! scale boundary velocity to adjacent fluid velocity based on the
! concentration ratio of fluid cell to boundary cell. This ratio is most
! likely 1 except for compressible cases with a PO boundary where P_g
! of the boundary is set and may differ from the value of the adjacent
! fluid cell
                     U_G(IJK) = ROP_G(FIJK)*U_G(FIJK)/ROP_G(IJK)
                  ELSE
                     U_G(IJK) = ZERO
                  ENDIF
                  V_G(IJK) = V_G(FIJK)
                  W_G(IJK) = W_G(FIJK)

                  IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
                     IF (MMAX > 0) THEN
                        WHERE (ROP_S(IJK,:MMAX) > ZERO)
! scale boundary velocity to adjacent fluid velocity based on the
! concentration ratio of fluid cell to boundary cell. Since solids are
! incompressible this is likely unnecessary. However, differences may
! arise if bc_rop_s is set.
                           U_S(IJK,:MMAX) = ROP_S(FIJK,:MMAX)*&
                              U_S(FIJK,:MMAX)/ROP_S(IJK,:MMAX)
                        ELSEWHERE
                           U_S(IJK,:MMAX) = ZERO
                        END WHERE
                        V_S(IJK,:MMAX) = V_S(FIJK,:MMAX)
                        W_S(IJK,:MMAX) = W_S(FIJK,:MMAX)
                     ENDIF
                  ENDIF   ! end if (.not.discrete_element .or.
                          !         des_continuum_hybrid)
               ENDIF   ! end if (fluid_at(im_of(ijk)))
! ----------------------------------------------------------------<<<


! Fluid cell at East
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(IP_OF(IJK))) THEN
                  FIJK = IP_OF(IJK)
                  RVEL_G = -U_G(IJK)
                  DO M = 1,MMAX
                     RVEL_S(M) = -U_S(FIJK,M)
                  ENDDO

                  CALL SET_OUTFLOW2(BCV, IJK, FIJK, RVEL_G, RVEL_S)

! provide an initial value for the velocity component through the domain
! otherwise its present value (from solution of the corresponding
! momentum eqn) is kept. values for the velocity components in the off
! directions are modified
                  IF (U_G(IJK) == UNDEFINED) THEN
                     IF (ROP_G(IJK) > ZERO) THEN
                        U_G(IJK) = ROP_G(FIJK)*U_G(FIJK)/ROP_G(IJK)
                     ELSE
                        U_G(IJK) = ZERO
                     ENDIF
                  ENDIF
                  V_G(IJK) = V_G(FIJK)
                  W_G(IJK) = W_G(FIJK)

                  IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
                     DO M = 1, MMAX
                        IF (U_S(IJK,M) == UNDEFINED) THEN
                           IF (ROP_S(IJK,M) > ZERO) THEN
                              U_S(IJK,M) = ROP_S(FIJK,M)*&
                                 U_S(FIJK,M)/ROP_S(IJK,M)
                           ELSE
                              U_S(IJK,M) = ZERO
                           ENDIF
                        ENDIF
                        V_S(IJK,M) = V_S(FIJK,M)
                        W_S(IJK,M) = W_S(FIJK,M)
                     ENDDO
                  ENDIF   ! end if (.not.discrete_element .or.
                          !         des_continuum_hybrid)
               ENDIF   ! end if (fluid_at(ip_of(ijk)))
! ----------------------------------------------------------------<<<



! Fluid cell at South
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(JM_OF(IJK))) THEN
                  FIJK = JM_OF(IJK)
                  RVEL_G = V_G(FIJK)
                  DO M= 1,MMAX
                     RVEL_S(M) = V_S(FIJK,M)
                  ENDDO

                  CALL SET_OUTFLOW2(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  U_G(IJK) = U_G(FIJK)
                  IF (ROP_G(IJK) > ZERO) THEN
                     V_G(IJK) = ROP_G(FIJK)*V_G(FIJK)/ROP_G(IJK)
                  ELSE
                     V_G(IJK) = ZERO
                  ENDIF
                  W_G(IJK) = W_G(FIJK)

                  IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
                     IF (MMAX > 0) THEN
                         U_S(IJK,:MMAX) = U_S(FIJK,:MMAX)
                         WHERE (ROP_S(IJK,:MMAX) > ZERO)
                            V_S(IJK,:MMAX) = ROP_S(FIJK,:MMAX)*&
                               V_S(FIJK,:MMAX)/ROP_S(IJK,:MMAX)
                         ELSEWHERE
                            V_S(IJK,:MMAX) = ZERO
                         END WHERE
                         W_S(IJK,:MMAX) = W_S(FIJK,:MMAX)
                      ENDIF
                  ENDIF   ! end if (.not.discrete_element .or.
                          !         des_continuum_hybrid)
               ENDIF   ! end if (fluid_at(jm_of(ijk)))
! ----------------------------------------------------------------<<<


! Fluid cell at North
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(JP_OF(IJK))) THEN
                  FIJK = JP_OF(IJK)
                  RVEL_G = -V_G(IJK)
                  DO M = 1,MMAX
                     RVEL_S(M) = -V_S(FIJK,M)
                  ENDDO

                  CALL SET_OUTFLOW2(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  U_G(IJK) = U_G(FIJK)
                  IF (V_G(IJK) == UNDEFINED) THEN
                     IF (ROP_G(IJK) > ZERO) THEN
                        V_G(IJK) = ROP_G(FIJK)*V_G(FIJK)/ROP_G(IJK)
                     ELSE
                        V_G(IJK) = ZERO
                     ENDIF
                  ENDIF
                  W_G(IJK) = W_G(FIJK)

                  IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
                     DO M = 1, MMAX
                        U_S(IJK,M) = U_S(FIJK,M)
                        IF (V_S(IJK,M) == UNDEFINED) THEN
                           IF (ROP_S(IJK,M) > ZERO) THEN
                              V_S(IJK,M) = ROP_S(FIJK,M)*&
                                 V_S(FIJK,M)/ROP_S(IJK,M)
                           ELSE
                              V_S(IJK,M) = ZERO
                           ENDIF
                        ENDIF
                        W_S(IJK,M) = W_S(FIJK,M)
                     ENDDO
                  ENDIF   ! end if (.not.discrete_element .or.
                          !         des_continuum_hybrid)
               ENDIF   ! if (fluid_at(jp_of(ijk)))
! ----------------------------------------------------------------<<<


! Fluid cell at Bottom
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(KM_OF(IJK))) THEN
                  FIJK = KM_OF(IJK)
                  RVEL_G = W_G(IJK)
                  DO M = 1,MMAX
                     RVEL_S = W_S(IJK,M)
                  ENDDO

                  CALL SET_OUTFLOW2(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  U_G(IJK) = U_G(FIJK)
                  V_G(IJK) = V_G(FIJK)
                  IF (ROP_G(IJK) > ZERO) THEN
                     W_G(IJK) = ROP_G(FIJK)*W_G(FIJK)/ROP_G(IJK)
                  ELSE
                     W_G(IJK) = ZERO
                  ENDIF

                  IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
                     IF (MMAX > 0) THEN
                        U_S(IJK,:MMAX) = U_S(FIJK,:MMAX)
                        V_S(IJK,:MMAX) = V_S(FIJK,:MMAX)
                        WHERE (ROP_S(IJK,:MMAX) > ZERO)
                           W_S(IJK,:MMAX) = ROP_S(FIJK,:MMAX)*&
                              W_S(FIJK,:MMAX)/ROP_S(IJK,:MMAX)
                        ELSEWHERE
                           W_S(IJK,:MMAX) = ZERO
                        END WHERE
                     ENDIF
                  ENDIF   ! end if (.not.discrete_element .or.
                          !         des_continuum_hybrid)
               ENDIF   ! if (fluid_at(km_of(ijk)))
! ----------------------------------------------------------------<<<


! Fluid cell at Top
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(KP_OF(IJK))) THEN
                  FIJK = KP_OF(IJK)
                  RVEL_G = -W_G(IJK)
                  DO M = 1, MMAX
                     RVEL_S = -W_S(IJK,M)
                  ENDDO

                  CALL SET_OUTFLOW2(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  U_G(IJK) = U_G(FIJK)
                  V_G(IJK) = V_G(FIJK)
                  IF (W_G(IJK) == UNDEFINED) THEN
                     IF (ROP_G(IJK) > ZERO) THEN
                        W_G(IJK) = ROP_G(FIJK)*W_G(FIJK)/ROP_G(IJK)
                     ELSE
                        W_G(IJK) = ZERO
                     ENDIF
                  ENDIF

                  IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
                     DO M = 1, MMAX
                        U_S(IJK,M) = U_S(FIJK,M)
                        V_S(IJK,M) = V_S(FIJK,M)
                        IF (W_S(IJK,M) == UNDEFINED) THEN
                           IF (ROP_S(IJK,M) > ZERO) THEN
                              W_S(IJK,M) = ROP_S(FIJK,M)*&
                                 W_S(FIJK,M)/ROP_S(IJK,M)
                           ELSE
                              W_S(IJK,M) = ZERO
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF   ! end if (.not.discrete_element .or.
                          !         des_continuum_hybrid)
               ENDIF   ! if (fluid_at(kp_of(ijk)))
! ----------------------------------------------------------------<<<

            ENDDO   ! end do (i=i1,i2)
         ENDDO   ! end do (j=j1,j2)
      ENDDO   ! end do (k=k1,k2)

      RETURN
      END SUBROUTINE SET_OUTFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_OUTFLOW2                                            C
!  Purpose: Set various quantities in outflow (PO, MO, O) boundaries   C
!     according to their value in the adjacent fluid cell that would   C
!     not otherwise be set in these boundaries.                        C
!     For many of the scalar field variables (i.e., T_g, T_s, X_g,     C
!     X_s, k_turb_g, e_turb_g, theta_m and scalar) this routine        C
!     serves to setup initial values in outflow cells since for these  C
!     variables the respective governing equation's solver routine     C
!     will also set its own value in the outflow boundary cells -      C
!     making this routine redundant in that aspect. For the other      C
!     variables (e.g., P_star, P_s, MW_MIX_g, P_g, RO_g, ROP_s, EP_g   C
!     ROP_g, etc.) this routine does set their value in the outflow    C
!     cells since no othe routine will.                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_OUTFLOW2(BCV, IJK, FIJK, RVEL_G, RVEL_S)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE bc
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE scalars
      USE run
      USE compar
      USE mflux
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Boundary condition number
      INTEGER, INTENT(IN) :: BCV
! ijk index for boundary cell
      INTEGER, INTENT(IN) :: IJK
! ijk index for adjacent fluid cell
      INTEGER, INTENT(IN) :: FIJK
! the gas or solids velocity in the fluid cell adjacent to the boundary
! cell dot with the outward normal of that bc plane; defines the gas or
! solids velocity component normal to the bc plane as positive when it
! is flowing into the bc cell from the fluid cell. so, for example, for
! an outflow on the eastern boundary this is the u component of velocity
! while for an outflow on the western boundary this is the -u component,
! etc.
      DOUBLE PRECISION, INTENT(IN) :: RVEL_G
      DOUBLE PRECISION, INTENT(IN), DIMENSION(DIMENSION_M) :: RVEL_S
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: M, N
! solids volume fraction
      DOUBLE PRECISION :: EPs
! sum of solids phases volume fractions
      DOUBLE PRECISION :: SUM_EPs
! sum of solids phases bulk densities
      DOUBLE PRECISION :: SUM_ROPS
!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: EOSG
!-----------------------------------------------

      IF (RVEL_G >=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN
! initially ep_g may be undefined (the initial step of a new run) but
! otherwise ep_g should always be defined. so this if effectively checks
! for backflow, and if backflow occurs skip these assignments.

         IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(FIJK)
         T_G(IJK) = T_G(FIJK)
         IF (NMAX(0) > 0) &
            X_G(IJK,:NMAX(0)) = X_G(FIJK,:NMAX(0))
         MW_MIX_G(IJK) = MW_MIX_G(FIJK)
! At this point, P_g, T_g and MW_MIX_G have been defined at IJK based on
! their values at FIJK
         IF (RO_G0 == UNDEFINED) RO_G(IJK) = &
            EOSG(MW_MIX_G(IJK),P_G(IJK),T_G(IJK))
      ENDIF
      P_STAR(IJK) = P_STAR(FIJK)

! setting scalar quantities
      DO N = 1, NScalar
         M = Phase4Scalar(N)
         IF(M == 0)Then
            IF (RVEL_G>=ZERO) THEN
               Scalar(IJK, N) = Scalar(FIJK, N)
            ENDIF
         ELSE
            IF (RVEL_S(M)>=ZERO) THEN
               Scalar(IJK, N) = Scalar(FIJK, N)
            ENDIF
         ENDIF
      ENDDO

! setting turbulence quantities
      IF(K_Epsilon) THEN
         IF (RVEL_G >= ZERO) THEN
            K_Turb_G(IJK) = K_Turb_G(FIJK)
            E_Turb_G(IJK) = E_Turb_G(FIJK)
         ENDIF
      ENDIF

! initializing solids quantities
      SUM_ROPS = ZERO
      SUM_EPS = ZERO

! DEM simulations do not employ variables for continuum solids.
      IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
         DO M = 1, SMAX
            P_S(IJK,M) = P_S(FIJK,M)
! check if solids are entering boundary cell from adjacent fluid cell,
! if so set boundary cell values of scalar quantities rop_s, t_s and
! theta_m according to their values in the adjacent fluid cell.
            IF (RVEL_S(M) >= ZERO) THEN
               ROP_S(IJK,M) = ROP_S(FIJK,M)
               T_S(IJK,M) = T_S(FIJK,M)
               THETA_M(IJK,M) =  THETA_M(FIJK,M)
               IF (NMAX(M) > 0) &
                  X_S(IJK,M,:NMAX(M)) = X_S(FIJK,M,:NMAX(M))
            ELSE
               ROP_S(IJK,M) = ZERO
! why not maintain the current value of theta_m (i.e. do nothing..)
! cannot have theta_m = zero anywhere
               THETA_M(IJK,M) =  THETA_M(FIJK,M)
            ENDIF

! if bc_rop_s is defined, set value of rop_s in the ijk boundary cell
! according to user definition
            IF(BC_ROP_S(BCV,M)/=UNDEFINED) ROP_S(IJK,M)=BC_ROP_S(BCV,M)

! add to total solids phase bulk density and solids volume fraction
            SUM_ROPS = SUM_ROPS + ROP_S(IJK,M)
            SUM_EPS = SUM_EPS + EP_S(IJK,M)
         ENDDO   ! end do (m=1,smax)
      ENDIF   ! end if (.not.discrete_element .or. des_continuum_hybrid)

      IF (KT_TYPE_ENUM == GHD_2007) THEN
         P_S(IJK,MMAX) = P_S(FIJK,MMAX)
         THETA_M(IJK,MMAX) =  THETA_M(FIJK,MMAX)
         ROP_S(IJK,MMAX) = SUM_ROPS
      ENDIF

! this section must be skipped until after the initial setup of the
! discrete element portion of the simulation (set_bc1 is called once
! before the initial setup).
      IF (DISCRETE_ELEMENT .AND. ALLOCATED(DES_ROP_S)) THEN
         DO M = 1, DES_MMAX
! unlike in the two fluid model, in the discrete element model it is
! possible to actually calculate the bulk density in a flow boundary
! cell. Currently, however, such calculations are not strictly enforced.
! therefore use the bulk density of the adjacent fluid cell
            DES_ROP_S(IJK,M) = DES_ROP_S(FIJK,M)
            SUM_ROPS = SUM_ROPS + DES_ROP_S(IJK,M)
            EPS = DES_ROP_S(IJK,M)/DES_RO_S(M)
            SUM_EPS = SUM_EPS + EPS
         ENDDO
      ENDIF

! if bc_ep_g undefined, set ep_g to 1
      IF (BC_EP_G(BCV) == UNDEFINED) &
         EP_G(IJK) = ONE - SUM_EPS

! now that ep_g in the boundary cell is known, define the bulk density
! of the gas phase in the boundary cell
      ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)

! set boundary cell values of convective fluxes according to adjacent
! fluid cell values. Should these assigned be more like velocity?
      Flux_gE(IJK) = Flux_gE(FIJK)
      Flux_gN(IJK) = Flux_gN(FIJK)
      Flux_gT(IJK) = Flux_gT(FIJK)
      IF(ADDED_MASS) THEN
        Flux_gSE(IJK) = Flux_gSE(FIJK)
        Flux_gSN(IJK) = Flux_gSN(FIJK)
        Flux_gST(IJK) = Flux_gST(FIJK)
      ENDIF

      IF (MMAX >0) THEN
         Flux_sE(IJK,:MMAX) = Flux_sE(FIJK,:MMAX)
         Flux_sN(IJK,:MMAX) = Flux_sN(FIJK,:MMAX)
         Flux_sT(IJK,:MMAX) = Flux_sT(FIJK,:MMAX)
         IF(ADDED_MASS) THEN
           Flux_sSE(IJK) = Flux_sSE(FIJK)
           Flux_sSN(IJK) = Flux_sSN(FIJK)
           Flux_sST(IJK) = Flux_sST(FIJK)
         ENDIF
      ENDIF

      IF (KT_TYPE_ENUM == GHD_2007) THEN
         Flux_nE(IJK) = Flux_nE(FIJK)
         Flux_nN(IJK) = Flux_nN(FIJK)
         Flux_nT(IJK) = Flux_nT(FIJK)
      ENDIF

      RETURN
      END SUBROUTINE SET_OUTFLOW2

