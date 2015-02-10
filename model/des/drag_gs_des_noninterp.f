!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_GS_DES_NONINTERP                                   !
!                                                                      !
!  Purpose: This routine is called from the DISCRETE side to calculate !
!  the drag force acting on each particle using the scalar cell center !
!  gas velocity. The total contact force is also updated to include    !
!  the gas pressure force.                                             !
!                                                                      !
!  Notes:                                                              !
!                                                                      !
!   F_gp is obtained from des_drag_gp subroutine is given as:          !
!    F_GP = beta*VOL_P/EP_s where VOL_P is the particle volume.        !
!                                                                      !
!  The drag force on each particle is equal to:                        !
!    D_FORCE = beta*VOL_P/EP_s*(Ug - Us) = F_GP *(Ug - Us)             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_GS_DES_NONINTERP

! Global Variables:
!---------------------------------------------------------------------//
! The count and a list of particles in IJK
      use discretelement, only: PINC, PIC
! Flags indicating the state of particle
      use discretelement, only: PEA
! Drag force on each particle
      use discretelement, only: F_GP
! Particle velocity
      use discretelement, only: DES_VEL_NEW
! Total forces acting on particle
      use discretelement, only: FC
! Gas pressure force by fluid cell
      use discretelement, only: P_FORCE
! Particle volume.
      use discretelement, only: PVOL
! Model B momentum equation
      use run, only: MODEL_B
! Flag for MPPIC runs.
      use mfix_pic, only: MPPIC
! Flag to use implicit drag for MPPIC
      use mfix_pic, only: MPPIC_PDRAG_IMPLICIT
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to deterine if a cell contains fluid.
      use functions, only: FLUID_AT

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
! local variable used for debugging
! general i, j, k indices
      INTEGER :: IJK
! particle number index, used for looping
      INTEGER :: NP, NINDX
! Scalar center gas velocity
      DOUBLE PRECISION :: VelFp(3)
! Drag force acting on each particle.
      DOUBLE PRECISION :: D_FORCE(3)
! Flag for Model A momentum equation
      LOGICAL :: MODEL_A

! Set flag for Model A momentum equation.
      MODEL_A = .NOT.MODEL_B

! Calculate the drag on each particle using scalar center gas velocity
!---------------------------------------------------------------------//
!!$omp parallel do default(none)                                        &
!!$omp shared(IJKSTART3, IJKEND3, PINC, PEA, PIC, DES_VEL_NEW, F_GP,    &
!!$omp   FC, P_FORCE, PVOL, MPPIC, MPPIC_PDRAG_IMPLICIT, MODEL_A)       &
!!$omp private(IJK, VELFP, NINDX, NP, D_FORCE)
      DO IJK = IJKSTART3,IJKEND3

         IF(.NOT.FLUID_AT(IJK)) CYCLE
         IF(PINC(IJK) == 0) CYCLE

! Calculate the average fluid velocity at scalar cell center.
         CALL CALC_NONINTERP_VELFP_GAS(IJK, VELFP)

! Calculate the drag force for each particle in the current cell.
         DO NINDX = 1,PINC(IJK)
            NP = PIC(IJK)%P(NINDX)
! skipping indices that do not represent particles and ghost particles
            IF(.NOT.PEA(NP,1)) CYCLE
            IF(PEA(NP,4)) CYCLE

! Calculate the particle centered drag coefficient (F_GP) using the
! particle velocity and the cell averaged gas velocity.
            CALL DES_DRAG_GP(NP, VelFp, DES_VEL_NEW(:,NP))

! Calculate the gas-solids drag force on the particle
            IF(MPPIC .AND. MPPIC_PDRAG_IMPLICIT) THEN
! implicit treatment of the drag term for mppic
               D_FORCE(:) = F_GP(NP)*(VelFp)
            ELSE
! default case
               D_FORCE(:) = F_GP(NP)*(VelFp-DES_VEL_NEW(:,NP))
            ENDIF

! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag
            FC(:,NP) = FC(:,NP) + D_FORCE(:)

! P_force is evaluated as -dp/dx
            IF(MODEL_A) FC(:,NP) = FC(:,NP) + P_FORCE(:,IJK)*PVOL(NP)
         ENDDO ! end do (nindx = 1,pinc(ijk))

      ENDDO ! end do (ijk=ijkstart3,ijkend3)
!!$omp end parallel do

      RETURN
      END SUBROUTINE DRAG_GS_DES_NONINTERP



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_GS_GAS_NONINTERP                                   !
!                                                                      !
!  Purpose: This routine is called from the CONTINUUM. It calculates   !
!  the scalar cell center drag force acting on the fluid using the     !
!  scalar cell center gas velocity.                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_GS_GAS_NONINTERP

! Global Variables:
!---------------------------------------------------------------------//
! The count and a list of particles in IJK
      use discretelement, only: PINC, PIC
! Flags indicating the state of particle
      use discretelement, only: PEA
! Drag force on each particle
      use discretelement, only: F_GP
! Particle velocity
      use discretelement, only: DES_VEL_NEW
! Contribution to gas momentum equation due to drag
      use discretelement, only: DRAG_BM
! Scalar cell center total drag force
      use discretelement, only: F_GDS
! Flag for MPPIC runs
      use mfix_pic, only: MPPIC
! Statical weight of each MPPIC parcel
      use mfix_pic, only: DES_STAT_WT
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to deterine if a cell contains fluid.
      use functions, only: FLUID_AT
! Volume of scalar cell.
      use geometry, only: VOL

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO, ONE

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
! Loop counter for fluid grid
      INTEGER :: IJK
! Loop counters for particles
      INTEGER :: NP, NINDX
! Scalar cell center fluid velocity
      DOUBLE PRECISION :: VelFp(3)
! One divided by fluid cell volume
      DOUBLE PRECISION :: OoVOL
! Drag force (intermediate calculation)
      DOUBLE PRECISION :: lFORCE

!......................................................................!

! Calculate the drag for each fluid cell if it contains particles.
!---------------------------------------------------------------------//
!!$omp parallel do default(none)                                        &
!!$omp shared(IJKSTART3, IJKEND3, PINC, PEA, PIC, DES_VEL_NEW, F_GP,    &
!!$omp   VOL, F_GDS, DRAG_BM, MPPIC, DES_STAT_WT)                       &
!!$omp private(IJK, VELFP, NINDX, NP, OoVol, lFORCE)
      DO IJK = IJKSTART3,IJKEND3

! Initialize fluid cell values.
         F_GDS(IJK) = ZERO
         DRAG_BM(IJK,:) = ZERO

! Skip non-fluid cells and cells without particles.
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         IF(PINC(IJK) == 0)  CYCLE

! Calculate the average fluid velocity at scalar cell center.
         CALL CALC_NONINTERP_VELFP_GAS(IJK, VELFP)

         OoVOL = ONE/VOL(IJK)

! loop through particles in the cell
         DO NINDX = 1,PINC(IJK)
            NP = PIC(IJK)%P(NINDX)
! skipping indices that do not represent particles and ghost particles
            IF(.NOT.PEA(NP,1)) CYCLE
            IF(PEA(NP,4)) CYCLE

! Calculate the particle centered drag coefficient (F_GP) using the
! particle velocity and the cell averaged gas velocity.
            CALL DES_DRAG_GP(NP, VelFp, DES_VEL_NEW(:,NP))

            lFORCE = OoVOL*F_GP(NP)
            IF(MPPIC) lFORCE = lFORCE*DES_STAT_WT(NP)

            F_GDS(IJK) = F_GDS(IJK) + lFORCE
            DRAG_BM(IJK,:) = DRAG_BM(IJK,:) + lFORCE*DES_VEL_NEW(:,NP)
         ENDDO

      ENDDO
!!$omp end parallel do


      RETURN
      END SUBROUTINE DRAG_GS_GAS_NONINTERP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_GS_EXPLICIT_NONINTERP                              !
!                                                                      !
!  Purpose: This routine is called from the CONTINUUM. It calculates   !
!  the scalar cell center drag force acting on the fluid using the     !
!  scalar cell center gas velocity.                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_GS_EXPLICIT_NONINTERP

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase volume fraction
      use fldvar, only: EP_G
! Max number of particles on current process
      use discretelement, only: MAX_PIP
! Flags indicating the state of particle
      use discretelement, only: PEA, PIJK
! Drag force on each particle
      use discretelement, only: F_GP
! Particle drag force
      use discretelement, only: DRAG_FC
! Particle velocity
      use discretelement, only: DES_VEL_NEW
! Contribution to gas momentum equation due to drag
      use discretelement, only: DRAG_BM
! Scalar cell center total drag force
      use discretelement, only: F_GDS
! Flag for MPPIC runs
      use mfix_pic, only: MPPIC
! Statical weight of each MPPIC parcel
      use mfix_pic, only: DES_STAT_WT
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to deterine if a cell contains fluid.
      use functions, only: FLUID_AT
! Volume of scalar cell.
      use geometry, only: VOL
! Cell-center gas velocities.
      use tmp_array, only: UGC => ARRAY1
      use tmp_array, only: VGC => ARRAY2
      use tmp_array, only: WGC => ARRAY3
! Lock/Unlock the temp arrays to prevent double usage.
      use tmp_array, only: LOCK_TMP_ARRAY
      use tmp_array, only: UNLOCK_TMP_ARRAY
! MPI wrapper for halo exchange.
      use sendrecv, only: SEND_RECV


! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO, ONE

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
! Loop counter for fluid grid
      INTEGER :: IJK
! Loop counters for particles
      INTEGER :: NP, NINDX
! Scalar cell center fluid velocity
      DOUBLE PRECISION :: VelFp(3)
! Drag force (intermediate calculation)
      DOUBLE PRECISION :: lFORCE

!......................................................................!

! Initialize fluid cell values.
      F_GDS = ZERO
      DRAG_BM = ZERO

! Lock the temp arrays.
      CALL LOCK_TMP_ARRAY

! Calculate the cell center gas velocities.
      CALL CALC_CELL_CENTER_GAS_VEL

! Calculate the drag for each fluid cell if it contains particles.
!---------------------------------------------------------------------//
      DO NP=1,MAX_PIP
         IF(.NOT.PEA(NP,1)) CYCLE

! The drag force is not calculated on entering or exiting particles
! as their velocities are fixed and may exist in 'non fluid' cells.
         IF(any(PEA(NP,2:3))) THEN
            DRAG_FC(:,NP) = ZERO
            CYCLE
         ENDIF

! Fluid index containing the particle
         IJK = PIJK(NP,4)

! Gas phase velocity
         VELFP(1) = UGC(IJK)
         VELFP(2) = VGC(IJK)
         VELFP(3) = WGC(IJK)

         CALL DES_DRAG_GP_NEW(NP, DES_VEL_NEW(:,NP), VELFP, EP_g(IJK))

! Evaluate the drag force acting on the particle.
         DRAG_FC(:,NP) = F_GP(NP)*(VELFP - DES_VEL_NEW(:,NP))

         lFORCE = F_GP(NP)/VOL(IJK)
         IF(MPPIC) lFORCE = lFORCE*DES_STAT_WT(NP)

         F_GDS(IJK) = F_GDS(IJK) + lFORCE

         DRAG_BM(IJK,:) = DRAG_BM(IJK,:) +                          &
            lFORCE*(DES_VEL_NEW(:,NP)-VELFP)

      ENDDO

! Unlock the temp arrays.
      CALL UNLOCK_TMP_ARRAY

! Update the drag force and sources in ghost layers.
      CALL SEND_RECV(F_GDS, 2)
      CALL SEND_RECV(DRAG_BM, 2)


      RETURN
      END SUBROUTINE DRAG_GS_EXPLICIT_NONINTERP



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_NONINTER_VELFP                                     !
!  Author: J.Musser                                   Date: 07-NOV-14  !
!                                                                      !
!  Purpose: Calculate the scalar cell center gas velocity. This code   !
!  is common to the DEM and GAS calls for non-interpolated drag        !
!  routines.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_NONINTERP_VELFP_GAS(IJK, lVELFP)

! Global Variables:
!---------------------------------------------------------------------//
! Fluid velocities.
      use fldvar, only: U_G, V_G, W_G
! Functions to average momentum to scalar cell center.
      use fun_avg, only: AVG_X_E, AVG_Y_N, AVG_Z_T
! Flags and correction factors for cut momentum cells.
      use cutcell, only: CUT_U_TREATMENT_AT, THETA_UE, THETA_UE_BAR
      use cutcell, only: CUT_V_TREATMENT_AT, THETA_VN, THETA_VN_BAR
      use cutcell, only: CUT_W_TREATMENT_AT, THETA_WT, THETA_WT_BAR
! Functions to lookup adjacent cells by index.
      use functions, only: I_OF, IM_OF, JM_OF, KM_OF
! Flag for 3D simulatoins.
      use geometry, only: DO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision parameters.
      use param1, only: ZERO

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Fluid cell index
      INTEGER, INTENT(IN) :: IJK
! Fluid velocity vector at IJK cell center.
      DOUBLE PRECISION, INTENT(OUT) :: lVELFP(3)

! Local variables:
!---------------------------------------------------------------------//
! Indices of adjacent cells
      INTEGER :: IMJK, IJMK, IJKM


! Calculate the average fluid velocity at scalar cell center.
      IMJK = IM_OF(IJK)
      IF(CUT_U_TREATMENT_AT(IMJK)) THEN
         lVELFP(1) = (THETA_UE_BAR(IMJK)*U_G(IMJK) +                 &
            THETA_UE(IMJK)*U_G(IJK))
      ELSE
         lVELFP(1) = AVG_X_E(U_G(IMJK),U_G(IJK),I_OF(IJK))
      ENDIF

      IJMK = JM_OF(IJK)
      IF(CUT_V_TREATMENT_AT(IJMK)) THEN
         lVELFP(2) = (THETA_VN_BAR(IJMK)*V_G(IJMK) +                 &
            THETA_VN(IJMK)*V_G(IJK))
      ELSE
         lVELFP(2) = AVG_Y_N(V_G(IJMK),V_G(IJK))
      ENDIF

      IF(DO_K) THEN
         IJKM = KM_OF(IJK)
         IF(CUT_W_TREATMENT_AT(IJKM)) THEN
            lVELFP(3) = (THETA_WT_BAR(IJKM)*W_G(IJKM) +            &
               THETA_WT(IJKM)* W_G(IJK))
         ELSE
            lVELFP(3) = AVG_Z_T(W_G(IJKM),W_G(IJK))
         ENDIF
      ELSE
         lVELFP(3) = ZERO
      ENDIF

      RETURN
      END SUBROUTINE CALC_NONINTERP_VELFP_GAS
