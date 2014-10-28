!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_SS_DEM_NONINTERP                                   !
!  Purpose: This subroutine is called from the calc_des_force          !
!           and is called from the DISCRETE phase.                     !
!           This subroutine calculates the drag force exerted on the   !
!           particles by a continuous solids phase using cell average  !
!           quantities (i.e., non-interpolated version). The drag      !
!           coefficient (F_SS) is calculated here during the discrete  !
!           time step(s) and also during the continuum time step.      !
!           Accordingly, the drag coefficient in each call will be     !
!           based on the most currently available values (see notes)   !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_SS_DEM_NONINTERP

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE drag
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      use desmpi
      USE cutcell
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
! see the discussion for IJK_U ..... in comments
      INTEGER :: IJK_U, IJK_V, IJK_W
! average solids velocity at scalar cell center in array form
      DOUBLE PRECISION :: VELCS_ARR(3), &
                          VELDS_ARR(3)
! local drag force
      DOUBLE PRECISION :: SS_DRAG (DIMENSION_3, DES_MMAX, 3)
! Index of continuum solids phase
      INTEGER :: CM, M
! Index of discrete solids 'phase'
      INTEGER :: DM, L
! particle number index, used for looping
      INTEGER :: NP
! relative velocity between solids phase m and l
      DOUBLE PRECISION :: VREL
! particle diameters of phase M and phase L
      DOUBLE PRECISION :: D_pm, D_pl
! particle densities of phase M and phase L
      DOUBLE PRECISION :: RO_M, RO_L
! radial distribution function between phase M and L
      DOUBLE PRECISION :: G0_ML
! Solids volume fraction, and void fraction
      DOUBLE PRECISION :: EP_SM, EPS, EPG
! Sum over all phases of ratio volume fraction over particle diameter
      DOUBLE PRECISION :: EPSoDP
! one over solids volume fraction in fluid cell and one over the
! volume of fluid cell
      DOUBLE PRECISION :: OEPS
! solid-solid drag coefficient
      DOUBLE PRECISION :: lDss
! for error messages
      INTEGER :: IER
!-----------------------------------------------
! External Functions
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: G_0
!-----------------------------------------------


! initializing
      SS_DRAG(:,:,:) = ZERO

!!$omp parallel do default(shared)                                 &
!!$omp private(ijk,i,j,k,imjk,ijmk,ijkm,ijk_u,ijk_v,ijk_w,         &
!!$omp         velds_arr,velcs_arr,vrel,ss_drag,ldss,              &
!!$omp         dm,cm,m,l,epg,eps,ep_sm,epsodp,oeps,                &
!!$omp         d_pm,d_pl,ro_l,ro_m,g0_ml)                          &
!!$omp schedule (guided,50)
      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IMJK = IM_OF(IJK)
         IJMK = JM_OF(IJK)
         IJKM = KM_OF(IJK)
         IJK_U = IMJK
         IJK_V = IJMK
         IJK_W = IJKM

         IF(PINC(IJK).GT.0) THEN

            DO DM = 1, DES_MMAX
               EP_SM = DES_ROP_S(IJK,DM)/DES_RO_S(DM)

               IF(EP_SM.GT.ZERO) THEN

! calculate the cell center average continuum solids velocities.
! no manipulation is needed for the discrete solids velocities
! since these are already determined at cell centers
! ---------------------------------------------------------------->>>

! defining array form of average solids velocity
                  VELDS_ARR(1) = DES_U_S(IJK,DM)
                  VELDS_ARR(2) = DES_V_S(IJK,DM)
                  VELDS_ARR(3) = merge(DES_W_S(IJK,DM), ZERO, DO_K)

                  DO CM = 1, MMAX
                     IF(CUT_U_TREATMENT_AT(IJK_U)) THEN
                        VELCS_ARR(1) = &
                           (Theta_Ue_bar(IJK_U)*U_S(IJK_U,CM) + &
                            Theta_Ue(IJK_U)    *U_S(IP_OF(IJK_U),CM))
                     ELSE
                        VELCS_ARR(1) = &
                           AVG_X_E(U_S(IMJK,CM),U_S(IJK,CM),I)
                     ENDIF
                     IF(CUT_V_TREATMENT_AT(IJK_V)) THEN
                        VELCS_ARR(2) = &
                           (Theta_Vn_bar(IJK_V)*V_S(IJK_V,CM) + &
                            Theta_Vn(IJK_V)    *V_S(JP_OF(IJK_V),CM))
                     ELSE
                        VELCS_ARR(2) = &
                           AVG_Y_N(V_S(IJMK,CM),V_S(IJK,CM))
                     ENDIF

! calculating the relative velocity in 2D (overwrite if 3D)
                     VREL = SQRT((VELCS_ARR(1)-VELDS_ARR(1))**2+&
                                 (VELCS_ARR(2)-VELDS_ARR(2))**2)

                     IF(DO_K) THEN
                        IF(CUT_W_TREATMENT_AT(IJK_W)) THEN
                           VELCS_ARR(3) = &
                              (Theta_Wt_bar(IJK_W)*W_S(IJK_W,CM) + &
                               Theta_Wt(IJK_W)    * W_S(KP_OF(IJK_W),CM))
                        ELSE
                           VELCS_ARR(3) = &
                              AVG_Z_T(W_S(IJKM,CM),W_S(IJK,CM))
                        ENDIF
! calculating the relative velocity in 3D (overwrite 2D calculation)
                        VREL = SQRT((VELCS_ARR(1)-VELDS_ARR(1))**2+&
                                    (VELCS_ARR(2)-VELDS_ARR(2))**2+&
                                    (VELCS_ARR(3)-VELDS_ARR(3))**2)
                     ELSE
                        VELCS_ARR(3) = ZERO
                     ENDIF
! ----------------------------------------------------------------<<<


! Update the solids-solids drag coefficient
! ---------------------------------------------------------------->>>
! setting aliases for easy reference
                     D_PM = D_P(IJK,CM)
                     D_PL = DES_D_P0(DM)
                     RO_M = RO_S(IJK,CM)
                     RO_L = DES_RO_S(DM)

! evaluating g0 - taken from G_0.f subroutine (lebowitz form)
! this section is needed to account for all solids phases until g0 for
! multiple solids types (i.e. discrete & continuum) can be addressed
! more effectively.
                     EPSoDP = ZERO
                     DO M = 1, MMAX
                        EPS = EP_s(IJK,M)
                        EPSoDP = EPSoDP + EPS / D_p(IJK,M)
                     ENDDO
                     DO L = 1, DES_MMAX
                        EPS = DES_ROP_S(IJK,L)/DES_RO_S(L)
                        EPSoDP = EPSoDP + EPS / DES_D_p0(L)
                     ENDDO
                     EPg = EP_g(IJK)
                     G0_ML = ONE/EPg + 3.0d0*EPSoDP*D_pM*D_PL / &
                        (EPg*EPg *(D_pM + D_pL))

                     CALL DRAG_SS_SYAM(lDss,D_PM,D_PL,RO_M,RO_L,G0_ML,VREL)

                     F_SDS(IJK,CM,DM) = lDss*ROP_S(IJK,CM)*&
                        DES_ROP_S(IJK,DM)

! accounting for particle-particle drag due to enduring contact in a
! close-packed system
                     IF(CLOSE_PACKED(CM)) F_SDS(IJK,CM,DM) = &
                        F_SDS(IJK,CM,DM) + &
                        SEGREGATION_SLOPE_COEFFICIENT*P_star(IJK)
! ----------------------------------------------------------------<<<


! calculating the accumulated solids-solids drag force on discrete
! solids phase dm
                     SS_DRAG(IJK,DM,:) = SS_DRAG(IJK,DM,:) - &
                        F_SDS(IJK,CM,DM)*(VELDS_ARR(:)-VELCS_ARR(:))

                  ENDDO   ! end do loop (cm=1,mmax)

                  OEPS = ONE/EP_SM
                  SS_DRAG(IJK,DM,:) = SS_DRAG(IJK,DM,:)*OEPS

               ENDIF  ! end if ep_sm>0
            ENDDO   ! end do loop (dm=1,des_mmax)

         ENDIF      ! end if(pinc(ijk).gt.0)

      ENDDO         ! end do loop (ijk=ijkstart3, ijkend3)
!!$omp end parallel do


!!$omp parallel do private(np,ijk,m,ss_drag)                       &
!!$omp schedule (guided,100)
      DO NP = 1, MAX_PIP
! skipping indices that do not represent particles and ghost particles
         if(.not.pea(np,1)) cycle
         if(pea(np,4)) cycle

         IJK = PIJK(NP,4)
         M = PIJK(NP,5)
! Update the contact forces (FC) on the particle to include
! solids-solids drag force
!----------------------------------------------------------------->>>
         FC(:,NP) = FC(:,NP) + SS_DRAG(IJK,M,:)*PVOL(NP)
!-----------------------------------------------------------------<<<

      ENDDO   ! end do loop (np=1,max_pip)
!!$omp end parallel do



      RETURN
      END SUBROUTINE DRAG_SS_DEM_NONINTERP



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DRAG_SS                                             C
!  Purpose: This subroutine is only called from the CONTINUUM side.    C
!     It performs the following functions                              C
!     - If non-interpolated then calculate the solids-solids drag      C
!       force coefficient between continuum and discrete solids.       C
!     - If interpolated, then it calculates the solids-particle        C
!       drag coefficient (F_SP) based on the particle velocity and     C
!       interpolated continuum solids velocity. It then determines     C
!       the contributions of solids-particle drag to the center        C
!       coefficient of the A matrix and the b (source) vector in the   C
!       matrix equation (A*VEL_FP=b) equation for the continuum        C
!       solids phase x, y and z momentum balances using F_SP           C
!                                                                      C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DRAG_SS_TFM_NONINTERP


!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE constant
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE drag
      USE discretelement
      use desmpi
      USE cutcell
      USE fun_avg
      USE functions
      IMPLICIT NONE

!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: I, IJK, IMJK, IJMK, IJKM, IJK_U, IJK_V, IJK_W
! Index of continuum solids phase
      INTEGER :: CM, M
! Index of discrete solids 'phase'
      INTEGER :: DM, L
! average solids velocity in x, y and z directions at scalar cell center
      DOUBLE PRECISION :: USCM, VSCM, WSCM, USDM, VSDM, WSDM
! average fluid and solid velocity in array form
      DOUBLE PRECISION :: VELG_ARR(3), &
                          VELDS_ARR(DES_MMAX, 3)
! relative velocity between solids phase m and l
      DOUBLE PRECISION :: VREL
! particle diameters of phase M and phase L
      DOUBLE PRECISION :: D_pm, D_pl
! particle densities of phase M and phase L
      DOUBLE PRECISION :: RO_M, RO_L
! radial distribution function between phase M and L
      DOUBLE PRECISION :: G0_ML
! Solids volume fraction, and void fraction
      DOUBLE PRECISION :: EPS, EPG
! Sum over all phases of ratio volume fraction over particle diameter
      DOUBLE PRECISION :: EPSoDP
! solid-solid drag coefficient
      DOUBLE PRECISION :: lDss
!-----------------------------------------------
! External Functions
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: G_0
!-----------------------------------------------

! NON-INTERPOLATED solid-solid drag:
! calculate the gas solids drag coefficient (F_SDS) using the cell
! averaged particle velocity and the cell average continuum solids
! velocity
!----------------------------------------------------------------->>>
!      IF (.NOT.DES_INTERP_ON) THEN
      DO CM = 1, MMAX
         DO DM = 1, DES_MMAX

            DO IJK = ijkstart3, ijkend3

               IF (FLUID_AT(IJK) .AND. PINC(IJK) >0) THEN   ! IF(.NOT.WALL_AT(IJK)) THEN?
                  I = I_OF(IJK)
                  IMJK = IM_OF(IJK)
                  IJMK = JM_OF(IJK)
                  IJKM = KM_OF(IJK)
                  IJK_U = IMJK
                  IJK_V = IJMK
                  IJK_W = IJKM

! calculating the cell center average continuum and discrete solids
! velocities. no manipulation is needed for the discrete soldis
! velocities since these are already determined at cell centers
                  IF(CUT_U_TREATMENT_AT(IJK_U)) THEN
                     USCM = (Theta_Ue_bar(IJK_U)*U_S(IJK_U,CM) + &
                            Theta_Ue(IJK_U)     *U_S(IP_OF(IJK_U),CM))
                  ELSE
                     USCM = AVG_X_E(U_S(IMJK,CM),U_S(IJK,CM),I)
                  ENDIF
                  IF(CUT_V_TREATMENT_AT(IJK_V)) THEN
                     VSCM = (Theta_Vn_bar(IJK_V)*V_S(IJK_V,CM) + &
                            Theta_Vn(IJK_V)     *V_S(JP_OF(IJK_V),CM))
                  ELSE
                     VSCM = AVG_Y_N(V_S(IJMK,CM),V_S(IJK,CM))
                  ENDIF
                  USDM = DES_U_S(IJK,DM)
                  VSDM = DES_V_S(IJK,DM)

! calculating the relative velocity in 2D (overwrite if 3D)
                  VREL = SQRT((USCM-USDM)**2 + (VSCM-VSDM)**2)

                  IF (DO_K) THEN
                     IF(CUT_W_TREATMENT_AT(IJK_W)) THEN
                        WSCM = (Theta_Wt_bar(IJK_W)*W_S(IJK_W,CM) + &
                               Theta_Wt(IJK_W)     *W_S(KP_OF(IJK_W),CM))
                     ELSE
                        WSCM = AVG_Z_T(W_S(IJKM,CM),W_S(IJK,CM))
                     ENDIF
                     WSDM = DES_W_S(IJK,DM)
! calculating the relative velocity in 3D
                     VREL = SQRT((USCM-USDM)**2 + (VSCM-VSDM)**2 +&
                                 (WSCM-WSDM)**2)
                  ELSE
                     WSDM = ZERO
                  ENDIF

! setting aliases for easy reference
                  D_PM = D_P(IJK,CM)
                  D_PL = DES_D_P0(DM)
                  RO_M = RO_S(IJK,CM)
                  RO_L = DES_RO_S(DM)

! evaluating g0 - taken from G_0.f subroutine (lebowitz form)
! this section is needed to account for all solids phases until g0 for
! multiple solids types (i.e. discrete & continuum) can be addressed
! more effectively.
                  EPSoDP = ZERO
                  DO M = 1, MMAX
                     EPS = EP_s(IJK, M)
                     EPSoDP = EPSoDP + EPS / D_p(IJK,M)
                  ENDDO
                  DO L = 1, DES_MMAX
                     EPS = DES_ROP_S(IJK,L)/DES_RO_S(L)
                     EPSoDP = EPSoDP + EPS / DES_D_p0(L)
                  ENDDO
                  EPg = EP_g(IJK)
                  G0_ML = ONE/EPg + 3.0d0*EPSoDP*D_pM*D_PL / &
                     (EPg*EPg *(D_pM + D_pL))

                  CALL DRAG_SS_SYAM(lDss,D_PM,D_PL,RO_M,RO_L,G0_ML,VREL)

                  F_SDS(IJK,CM,DM) = lDss*ROP_S(IJK,CM)*&
                     DES_ROP_S(IJK,DM)

! accounting for particle-particle drag due to enduring contact in a
! close-packed system.
                  IF(CLOSE_PACKED(CM)) F_SDS(IJK,CM,DM) = &
                     F_SDS(IJK,CM,DM) + &
                     SEGREGATION_SLOPE_COEFFICIENT*P_star(IJK)

               ELSE   ! else branch of if(fluid_at(ijk) .and. pinc(ijk)>0)

                  F_SDS(IJK,CM,DM) = ZERO

               ENDIF   ! end if/else (fluid_at(ijk))
            ENDDO    ! end do (ijk=ijkstart3,ijkend3)
         ENDDO   ! end do (dm=1,des_mmax)
      ENDDO  ! end do (cm=1,smax)

      RETURN

!      ENDIF   ! endif(.not.des_interp_on)
!-----------------------------------------------------------------<<<


! INTERPOLATED solids-solids drag (the rest of this routine):
! Calculate the solids-solids drag coefficient using the particle
! velocity and the solids velocity interpolated to particle
! position.
!----------------------------------------------------------------->>>
! initializations
!      sdrag_am = ZERO
!      sdrag_bm = ZERO

! interpolated version will require steps to
! - interpolate mth phase continuum solids velocity and solids volume
!   fraction to particle position
! - create a routine (des_drag_sp) to calculate the mth phase continuum
!   solids-particle drag coefficient (f_sp).
! - calculate contributions to a and b matrix in mth phase solids
!   momentum equations

!-----------------------------------------------------------------<<<
 
      END SUBROUTINE DRAG_SS_TFM_NONINTERP
