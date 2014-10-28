!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_GS_DES_NONINTERP                                   !
!  Purpose: This subroutine is only called from the DISCRETE side.     !
!     It performs the following functions.                             !
!     - Calculates the fluid-solids drag force exerted on the          !
!       particles by the fluid phase using cell average quantities.    !
!       The gas solids drag coefficient (F_GDS) is calculated from     !
!       the subroutine drag_gs, which is called during the continuum   !
!       time step (from calc_drag) and during the discrete time(s)     !
!       (from here).                                                   !
!     - The total contact force on the particle is then updated to     !
!       include the gas-solids drag force and gas pressure force       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_GS_DES_NONINTERP

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE drag
      use desmpi
      USE cutcell
      USE mfix_pic
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
! average fluid and solid velocity at scalar cell center in array form
      DOUBLE PRECISION :: VELG_ARR(3), &
                          VELDS_ARR(DES_MMAX, 3)
! local drag force
      DOUBLE PRECISION :: GS_DRAG (DIMENSION_3, DES_MMAX, 3)
! index of solid phase that particle NP belongs to
      INTEGER :: M
! particle number index, used for looping
      INTEGER :: NP
! solids volume fraction of phase M in fluid cell
      DOUBLE PRECISION :: EP_SM
! one over solids volume fraction in fluid cell and one over the
! volume of fluid cell
      DOUBLE PRECISION :: OEPS
! for error messages
      INTEGER :: IER
! Flag only used when the hybrid model is invoked and notifies the
! routine that the solid phase index M refers to the indice of a
! discrete 'phase' not a continuous phase so that the appropriate
! variables are referenced.
      LOGICAL :: DISCRETE_FLAG
!-----------------------------------------------


! initializing
      GS_DRAG(:,:,:) = ZERO

! computing F_GDS (gas-solids drag coefficient) with the latest average
! fluid and solid velocity fields.  see comments below
      DISCRETE_FLAG = .TRUE.   ! only matters if des_continuum_hybrid
      DO M = 1, DES_MMAX
         IF (RO_G0/=ZERO) THEN
! this call can not be readily replaced with a call to des_drag_gp
! due to difficulties going from particle phase and ijk loops to a
! particle loop & vice versa
            CALL DRAG_GS (M, DISCRETE_FLAG, IER)
         ENDIF
      ENDDO

!$omp parallel do default(shared)                                 &
!$omp private(ijk,i,j,k,imjk,ijmk,ijkm,ijk_u,ijk_v,ijk_w,         &
!$omp         velg_arr,velds_arr,                                 &
!$omp         m,oeps,ep_sm)                               &
!$omp schedule (guided,50)
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

! UGC (and likewise for VGC and WGC) are computed at the center of
! the scalar cell. The center of the Ith scalar cell can also be
! thought of as east face of the (I-1)th U- cell.
! For cut-cell, it is easier to think in terms of U, V, W grids/cell
! as the interpolation arays (like Theta_Ue_bar) are based on the
! respective grids, i.e., theta_Ue_bar is based on U- grid. See
! conv_diff_u_g for an example of this usage.
! For U uncut grid, the average at the east face of IJK_U will be
! U_AVG(at East face of IJK_U) = HALF*(U(IJK_U) + U(IP_OF(IJK_U)))
! It can be verified that the above formula and old formula of
! U_AVG(center of IJK) = HALF*(U(IJK) + U(IM_OF(IJK))) are identical
! since IJK_U = IM_OF(IJK)

         IF(PINC(IJK).GT.0) THEN

! average fluid velocity at scalar cell center in array form
            IF(CUT_U_TREATMENT_AT(IJK_U)) THEN
               VELG_ARR(1) = (Theta_Ue_bar(IJK_U)*U_G(IJK_U) + &
                              Theta_Ue(IJK_U)    *U_G(IP_OF(IJK_U)))
            ELSE
               VELG_ARR(1) = HALF * (U_G(IJK_U) + U_G(IP_OF(IJK_U)))
            ENDIF

            IF(CUT_V_TREATMENT_AT(IJK_V)) THEN
               VELG_ARR(2) = (Theta_Vn_bar(IJK_V)*V_G(IJK_V) + &
                              Theta_Vn(IJK_V)    *V_G(JP_OF(IJK_V)))
            ELSE
               VELG_ARR(2) = HALF * (V_G(IJK_V) + V_G(JP_OF(IJK_V)))
            ENDIF

            VELDS_ARR(:,1) = DES_U_S(IJK,:)
            VELDS_ARR(:,2) = DES_V_S(IJK,:)

            IF(DO_K) THEN
               IF(CUT_W_TREATMENT_AT(IJK_W)) THEN
                  VELG_ARR(3) = (Theta_Wt_bar(IJK_W)*W_G(IJK_W) + &
                                 Theta_Wt(IJK_W)    * W_G(KP_OF(IJK_W)))
               ELSE
                  VELG_ARR(3) = HALF * (W_G(IJK_W) + W_G(KP_OF(IJK_W)))
               ENDIF
               VELDS_ARR(:,3) = DES_W_S(IJK,:)
            ELSE
               VELG_ARR(3) = ZERO
               VELDS_ARR(:,3) = ZERO
            ENDIF


            DO M = 1, DES_MMAX
! the call to drag coefficient should probably be here for cut-cell
! since cell center velocities are now known and these would be used
! in vrel within drag correlations). this would require some
! rearrangement/reworking?  alternatively these calculations should be
! moved into drag_gs?

               EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)

               IF(EP_SM.GT.ZERO) THEN
                  IF (MPPIC .AND. MPPIC_PDRAG_IMPLICIT) THEN
                     GS_DRAG(IJK,M, :) = F_GDS(IJK,M)*VELG_ARR(:)
                  ELSEIF (DES_CONTINUUM_HYBRID) THEN
                     GS_DRAG(IJK,M,:) = -F_GDS(IJK,M)*&
                        (VELDS_ARR(M,:)-VELG_ARR(:))
                  ENDIF   ! end if/else (des_continuum_hybrid)
                  OEPS = ONE/EP_SM
                  GS_DRAG(IJK,M,:) = GS_DRAG(IJK,M,:)*OEPS
               ENDIF  ! end if ep_sm>0

            ENDDO   ! end do loop (dm=1,des_mmax)
         ENDIF      ! end if(pinc(ijk).gt.0)

      ENDDO         ! end do loop (ijk=ijkstart3, ijkend3)
!$omp end parallel do


!$omp parallel do private(np,ijk,m,oeps,ep_sm)             &
!$omp schedule (guided,100)
      DO NP = 1, MAX_PIP
! skipping indices that do not represent particles and ghost particles
         if(.not.pea(np,1)) cycle
         if(pea(np,4)) cycle

         IJK = PIJK(NP,4)
         M = PIJK(NP,5)

! Update the contact forces (FC) on the particle to include
! gas pressure and gas-solids drag
!----------------------------------------------------------------->>>
         FC(:,NP) = FC(:,NP) + GS_DRAG(IJK,M,:)*PVOL(NP)

         IF(.NOT.MODEL_B) THEN
! Add the pressure gradient force
! P_force is evaluated as -dp/dx
            FC(:,NP) = FC(:,NP) + (P_FORCE(IJK,:))*PVOL(NP)
         ENDIF
!-----------------------------------------------------------------<<<

! For mppic calculate the gas-particle drag coefficient (F_GP).
!----------------------------------------------------------------->>>
         IF(MPPIC) THEN
            EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
            IF(EP_SM.GT.ZERO) THEN
               OEPS = ONE/EP_SM
            ELSE
               OEPS = ZERO
            ENDIF
            IF(MPPIC_PDRAG_IMPLICIT) THEN
               F_gp(NP) = (F_GDS(IJK,M)*PVOL(NP))*OEPS
            ELSE
               F_gp(NP) = ZERO
            ENDIF
         ENDIF   ! end if (mppic)
!-----------------------------------------------------------------<<<

      ENDDO   ! end do loop (np=1,max_pip)
!$omp end parallel do

      RETURN
      END SUBROUTINE DRAG_GS_DES_NONINTERP



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_GS_GAS_NONINTERP                                   !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side.    !
!  The code is directed to the subroutine drag_gs for the appropriate  !
!  calculations (calculation of the fluid-solids drag coefficient).    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_GS_GAS_NONINTERP

      use discretelement, only: DES_MMAX
      use physprop, only: RO_G0
      use param1, only: ZERO

      IMPLICIT NONE

! index of solid phase that particle NP belongs to
      INTEGER :: M
! Flag only used when the hybrid model is invoked and notifies the
! routine that the solid phase index M refers to the indice of a
! discrete 'phase' not a continuous phase so that the appropriate
! variables are referenced.
      LOGICAL :: DISCRETE_FLAG
! Error Flag.
      INTEGER :: IER


! Calculate the fluid solids drag coefficient (F_GDS) using the cell
! averaged particle velocity and the cell average fluid velocity
!---------------------------------------------------------------------//
      DISCRETE_FLAG = .TRUE.   ! only matters if des_continuum_hybrid
      DO M = 1, DES_MMAX
! this call can not be readily replaced with a call to des_drag_gp
! due to difficulties going from particle phase and ijk loops to a
! particle loop & vice versa
         IF (RO_G0/=ZERO) CALL DRAG_GS (M, DISCRETE_FLAG, IER)
      ENDDO


      RETURN
      END SUBROUTINE DRAG_GS_GAS_NONINTERP
