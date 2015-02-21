!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_GS_DES_INTERP1                                     !
!  Author: J.Musser                                   Date: 21-NOV-14  !
!                                                                      !
!  Purpose: This routine is called from the DISCRETE side to calculate !
!  the gas-based forces acting on each particle using interpolated     !
!  values for gas velocity, gas pressure, and gas volume fraction.     !
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
      SUBROUTINE DRAG_GS_DES_INTERP1

! Gas phase volume fraction
      use fldvar, only: EP_G
! Gas phase velocities
      use fldvar, only: U_G, V_G, W_G

      use discretelement, only: MAX_PIP
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_WEIGHT

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

! Cell-center gas velocities.
      use tmp_array, only: UGC => ARRAY1
      use tmp_array, only: VGC => ARRAY2
      use tmp_array, only: WGC => ARRAY3
! Lock/Unlock the temp arrays to prevent double usage.
      use tmp_array, only: LOCK_TMP_ARRAY
      use tmp_array, only: UNLOCK_TMP_ARRAY


! Flag for MPPIC runs.
      use mfix_pic, only: MPPIC
! Flag to use implicit drag for MPPIC
      use mfix_pic, only: MPPIC_PDRAG_IMPLICIT
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Flag for 3D simulatoins.
      use geometry, only: DO_K
! Function to deterine if a cell contains fluid.
      use functions, only: FLUID_AT

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO


      IMPLICIT NONE

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, LC
! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
! Interpolated gas phase quanties.
      DOUBLE PRECISION :: lEPg, VELFP(3), lPF(3)
! Drag force acting on each particle.
      DOUBLE PRECISION :: D_FORCE(3)
! Flag for Model A momentum equation
      LOGICAL :: MODEL_A
! Loop bound for filter
      INTEGER :: LP_BND


! Set flag for Model A momentum equation.
      MODEL_A = .NOT.MODEL_B
! Loop bounds for interpolation.
      LP_BND = merge(27,9,DO_K)

! Lock the temp arrays.
      CALL LOCK_TMP_ARRAY

! Calculate the cell center gas velocities.
      CALL CALC_CELL_CENTER_GAS_VEL

! Calculate the gas phae forces acting on each particle.
      DO NP=1,MAX_PIP
         IF(.NOT.PEA(NP,1)) CYCLE
         IF(any(PEA(NP,2:4))) CYCLE

         lEPG = ZERO
         VELFP = ZERO
         lPF = ZERO

         DO LC=1,LP_BND

            IJK = FILTER_CELL(LC,NP)
            WEIGHT = FILTER_WEIGHT(LC,NP)
! Gas phase volume fraction.
            lEPG = lEPG + EP_G(IJK)*WEIGHT
! Gas phase velocity.
            VELFP(1) = VELFP(1) + UGC(IJK)*WEIGHT
            VELFP(2) = VELFP(2) + VGC(IJK)*WEIGHT
            VELFP(3) = VELFP(3) + WGC(IJK)*WEIGHT
! Gas pressure force.
            lPF = lPF + P_FORCE(:,IJK)*WEIGHT
         ENDDO

         CALL DES_DRAG_GP_NEW(NP, DES_VEL_NEW(:,NP), VELFP, lEPg)

! Calculate the gas-solids drag force on the particle
         IF(MPPIC .AND. MPPIC_PDRAG_IMPLICIT) THEN
! implicit treatment of the drag term for mppic
            D_FORCE = F_GP(NP)*VELFP
         ELSE
! default case
            D_FORCE = F_GP(NP)*(VELFP - DES_VEL_NEW(:,NP))
         ENDIF

! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag
         FC(:,NP) = FC(:,NP) + D_FORCE(:)

         IF(MODEL_A) FC(:,NP) = FC(:,NP) + lPF*PVOL(NP)

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

! Unlock the temp arrays.
      CALL UNLOCK_TMP_ARRAY

      RETURN
      END SUBROUTINE DRAG_GS_DES_INTERP1



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_GS_GAS_INTERP1                                     !
!  Author: J.Musser                                   Date: 21-NOV-14  !
!                                                                      !
!                                                                      !
!  Purpose: This routine is called from the CONTINUUM. It calculates   !
!  the scalar cell center drag force acting on the fluid using         !
!  interpolated values for the gas velocity and volume fraction. The   !
!  The resulting sources are interpolated back to the fluid grid.      !
!                                                                      !
!  NOTE: The loop over particles includes ghost particles so that MPI  !
!  communications are needed to distribute overlapping force between   !
!  neighboring grid cells. This is possible because only cells "owned" !
!  by the current process will have non-zero weights.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_GS_GAS_INTERP1

! Gas phase volume fraction
      use fldvar, only: EP_G
! Gas phase velocities
      use fldvar, only: U_G, V_G, W_G

      use discretelement, only: MAX_PIP
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_WEIGHT

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
! Volume of scalar cell.
      use geometry, only: VOL
! Flag for 3D simulatoins.
      use geometry, only: DO_K
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
      use param1, only: ZERO

      IMPLICIT NONE

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, LC
! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
! Interpolated gas phase quanties.
      DOUBLE PRECISION :: lEPg, VELFP(3)
! Loop bound for filter
      INTEGER :: LP_BND
! Drag force (intermediate calculation)
      DOUBLE PRECISION :: lFORCE
! Drag sources for fluid (intermediate calculation)
      DOUBLE PRECISION :: lDRAG_BM(3)


! Initialize fluid cell values.
      F_GDS = ZERO
      DRAG_BM = ZERO

! Loop bounds for interpolation.
      LP_BND = merge(27,9,DO_K)

! Lock the temp arrays.
      CALL LOCK_TMP_ARRAY

! Calculate the cell center gas velocities.
      CALL CALC_CELL_CENTER_GAS_VEL

! Calculate the gas phae forces acting on each particle.
      DO NP=1,MAX_PIP
         IF(.NOT.PEA(NP,1)) CYCLE

         lEPG = ZERO
         VELFP = ZERO

         DO LC=1,LP_BND
            IJK = FILTER_CELL(LC,NP)
            WEIGHT = FILTER_WEIGHT(LC,NP)
! Gas phase volume fraction.
            lEPG = lEPG + EP_G(IJK)*WEIGHT
! Gas phase velocity.
            VELFP(1) = VELFP(1) + UGC(IJK)*WEIGHT
            VELFP(2) = VELFP(2) + VGC(IJK)*WEIGHT
            VELFP(3) = VELFP(3) + WGC(IJK)*WEIGHT
         ENDDO

         CALL DES_DRAG_GP_NEW(NP, DES_VEL_NEW(:,NP), VELFP, lEPg)

         lFORCE = F_GP(NP)
         IF(MPPIC) lFORCE = lFORCE*DES_STAT_WT(NP)

         lDRAG_BM = lFORCE*DES_VEL_NEW(:,NP)

         DO LC=1,LP_BND
            IJK = FILTER_CELL(LC,NP)
            WEIGHT = FILTER_WEIGHT(LC,NP)/VOL(IJK)

            DRAG_BM(IJK,:) = DRAG_BM(IJK,:) + lDRAG_BM*WEIGHT
            F_GDS(IJK) = F_GDS(IJK) + lFORCE*WEIGHT
         ENDDO

      ENDDO

! Unlock the temp arrays.
      CALL UNLOCK_TMP_ARRAY

! Update the drag force and sources in ghost layers.
      CALL SEND_RECV(F_GDS, 2)
      CALL SEND_RECV(DRAG_BM, 2)


      RETURN
      END SUBROUTINE DRAG_GS_GAS_INTERP1

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_GS_EXPLICIT_INTERP1                                !
!  Author: J.Musser                                   Date: 21-NOV-14  !
!                                                                      !
!                                                                      !
!  Purpose: This routine is called from the CONTINUUM. It calculates   !
!  the scalar cell center drag force acting on the fluid using         !
!  interpolated values for the gas velocity and volume fraction. The   !
!  The resulting sources are interpolated back to the fluid grid.      !
!                                                                      !
!  NOTE: The loop over particles includes ghost particles so that MPI  !
!  communications are needed to distribute overlapping force between   !
!  neighboring grid cells. This is possible because only cells "owned" !
!  by the current process will have non-zero weights.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_GS_EXPLICIT_INTERP1

! Gas phase volume fraction
      use fldvar, only: EP_G
! Gas phase velocities
      use fldvar, only: U_G, V_G, W_G

      use discretelement, only: MAX_PIP
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_WEIGHT

! Flags indicating the state of particle
      use discretelement, only: PEA
! Drag force on each particle
      use discretelement, only: F_GP
! Particle velocity
      use discretelement, only: DES_VEL_NEW
! Particle drag force
      use discretelement, only: DRAG_FC
! Contribution to gas momentum equation due to drag
      use discretelement, only: DRAG_BM
! Scalar cell center total drag force
      use discretelement, only: F_GDS
! Flag for MPPIC runs
      use mfix_pic, only: MPPIC
! Statical weight of each MPPIC parcel
      use mfix_pic, only: DES_STAT_WT
! Volume of scalar cell.
      use geometry, only: VOL
! Flag for 3D simulatoins.
      use geometry, only: DO_K
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
      use param1, only: ZERO

      IMPLICIT NONE

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, LC
! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
! Interpolated gas phase quanties.
      DOUBLE PRECISION :: lEPg, VELFP(3)
! Loop bound for
      INTEGER :: LP_BND
! Drag force (intermediate calculation)
      DOUBLE PRECISION :: lFORCE
! Drag source for fluid (intermediate calculation)
      DOUBLE PRECISION :: lDRAG_BM(3)


! Initialize fluid cell values.
      F_GDS = ZERO
      DRAG_BM = ZERO

! Loop bounds for interpolation.
      LP_BND = merge(27,9,DO_K)

! Lock the temp arrays.
      CALL LOCK_TMP_ARRAY

! Calculate the cell center gas velocities.
      CALL CALC_CELL_CENTER_GAS_VEL

! Calculate the gas phase forces acting on each particle.
!---------------------------------------------------------------------//
!$omp parallel default(none)                                           &
!$omp private(np, lepg, velfp, lc, ijk, weight, lforce, ldrag_bm)      &
!$omp shared(max_pip, pea, lp_bnd, drag_fc, f_gp, des_vel_new, ugc,    &
!$omp   vgc, wgc, mppic, drag_bm, vol, f_gds, ep_g, filter_weight,     &
!$omp   des_stat_wt, filter_cell)
!$omp do
      DO NP=1,MAX_PIP
         IF(.NOT.PEA(NP,1)) CYCLE

! The drag force is not calculated on entering or exiting particles
! as their velocities are fixed and may exist in 'non fluid' cells.
        IF(any(PEA(NP,2:3))) THEN
            DRAG_FC(:,NP) = ZERO
            CYCLE
         ENDIF

         lEPG = ZERO
         VELFP = ZERO

         DO LC=1,LP_BND
            IJK = FILTER_CELL(LC,NP)
            WEIGHT = FILTER_WEIGHT(LC,NP)
! Gas phase volume fraction.
            lEPG = lEPG + EP_G(IJK)*WEIGHT
! Gas phase velocity.
            VELFP(1) = VELFP(1) + UGC(IJK)*WEIGHT
            VELFP(2) = VELFP(2) + VGC(IJK)*WEIGHT
            VELFP(3) = VELFP(3) + WGC(IJK)*WEIGHT
         ENDDO

         CALL DES_DRAG_GP_NEW(NP, DES_VEL_NEW(:,NP), VELFP, lEPg)

! Evaluate the drag force acting on the particle.
         DRAG_FC(:,NP) = F_GP(NP)*(VELFP - DES_VEL_NEW(:,NP))

! Calculate the force on the fluid.
         lFORCE = F_GP(NP)
         IF(MPPIC) lFORCE = lFORCE*DES_STAT_WT(NP)

         lDRAG_BM(:) = lFORCE*(DES_VEL_NEW(:,NP) - VELFP)

         DO LC=1,LP_BND
            IJK = FILTER_CELL(LC,NP)
            WEIGHT = FILTER_WEIGHT(LC,NP)/VOL(IJK)

!$omp atomic
            DRAG_BM(IJK,1) = DRAG_BM(IJK,1) + lDRAG_BM(1)*WEIGHT
!$omp atomic
            DRAG_BM(IJK,2) = DRAG_BM(IJK,2) + lDRAG_BM(2)*WEIGHT
!$omp atomic
            DRAG_BM(IJK,3) = DRAG_BM(IJK,3) + lDRAG_BM(3)*WEIGHT
!$omp atomic
            F_GDS(IJK) = F_GDS(IJK) + lFORCE*WEIGHT
         ENDDO

      ENDDO
!$omp end do
!$omp end parallel

! Unlock the temp arrays.
      CALL UNLOCK_TMP_ARRAY

! Update the drag force and sources in ghost layers.
      CALL SEND_RECV(F_GDS, 2)
      CALL SEND_RECV(DRAG_BM, 2)


      RETURN
      END SUBROUTINE DRAG_GS_EXPLICIT_INTERP1



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_CELL_CENTER_GAS_VEL                                !
!  Author: J.Musser                                   Date: 07-NOV-14  !
!                                                                      !
!  Purpose: Calculate the scalar cell center gas velocity. This code   !
!  is common to the DEM and GAS calls for non-interpolated drag        !
!  routines.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_CELL_CENTER_GAS_VEL

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
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Flag for 3D simulatoins.
      use geometry, only: DO_K
! Function to deterine if a cell contains fluid.
      use functions, only: FLUID_AT

      use tmp_array, only: UGC => ARRAY1
      use tmp_array, only: VGC => ARRAY2
      use tmp_array, only: WGC => ARRAY3

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision parameters.
      use param1, only: ZERO

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
! Indices of adjacent cells
      INTEGER :: IJK, IMJK, IJMK, IJKM

! Calculate the cell center gas velocity components.
      DO IJK=IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN
            IMJK = IM_OF(IJK)
            IF(CUT_U_TREATMENT_AT(IMJK)) THEN
               UGC(IJK) = (THETA_UE_BAR(IMJK)*U_G(IMJK) +              &
                  THETA_UE(IMJK)*U_G(IJK))
            ELSE
               UGC(IJK) = AVG_X_E(U_G(IMJK),U_G(IJK),I_OF(IJK))
            ENDIF

            IJMK = JM_OF(IJK)
            IF(CUT_V_TREATMENT_AT(IJMK)) THEN
               VGC(IJK) = (THETA_VN_BAR(IJMK)*V_G(IJMK) +              &
                  THETA_VN(IJMK)*V_G(IJK))
            ELSE
               VGC(IJK) = AVG_Y_N(V_G(IJMK),V_G(IJK))
            ENDIF

            IF(DO_K) THEN
               IJKM = KM_OF(IJK)
               IF(CUT_W_TREATMENT_AT(IJKM)) THEN
                  WGC(IJK) = (THETA_WT_BAR(IJKM)*W_G(IJKM) +           &
                     THETA_WT(IJKM)* W_G(IJK))
               ELSE
                  WGC(IJK) = AVG_Z_T(W_G(IJKM),W_G(IJK))
               ENDIF
            ELSE
               WGC(IJK) = ZERO
            ENDIF
         ELSE
            UGC(IJK) = ZERO
            VGC(IJK) = ZERO
            WGC(IJK) = ZERO
         ENDIF
      ENDDO


      RETURN
      END SUBROUTINE CALC_CELL_CENTER_GAS_VEL

