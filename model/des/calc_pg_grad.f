!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_PG_GRAD                                            !
!  Purpose: Calculate cell centered pressure force exerted on the      !
!           particles in the cell by the gas/fluid phase               !
!           Note that P_force is evaluated as -dp/dx                   !
!                                                                      !
!  Notes: This pressure force only needs to be calculated once during  !
!         the DEM loop (at the beginning) since the gas/fluid phase    !
!         is essentially static at that point (i.e., gas field is not  !
!         updated during DEM loop                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PG_GRAD

      use cutcell, only: CARTESIAN_GRID

! Model B momentum equation
      use run, only: MODEL_B

      use discretelement, only: MAX_PIP, PIJK, DES_EXPLICITLY_COUPLED
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_WEIGHT
      use particle_filter, only: DES_INTERP_ON

! Particle volume.
      use discretelement, only: PVOL
! Gas pressure force by fluid cell
      use discretelement, only: P_FORCE
! Particle drag force
      use discretelement, only: DRAG_FC
! Flag for 3D simulatoins.
      use geometry, only: DO_K

      use discretelement, only: IS_NONEXISTENT, IS_NORMAL, IS_ENTERING, IS_EXITING, IS_ENTERING_GHOST, IS_EXITING_GHOST

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO

      implicit none

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, LC
! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
! Interpolated gas phase quanties.
      DOUBLE PRECISION :: lPF(3)
! Loop bound for
      INTEGER :: LP_BND

      IF(CARTESIAN_GRID) THEN
         CALL CALC_PG_GRAD_CG
      ELSE
         CALL CALC_PG_GRAD_STD
      ENDIF

      IF(DES_EXPLICITLY_COUPLED .AND. .NOT.MODEL_B) THEN

! Loop bounds for interpolation.
         LP_BND = merge(27,9,DO_K)

! Calculate the gas phase forces acting on each particle.

!$omp  parallel do default(none) &
!$omp              private(NP,lPF,lc,ijk,weight) &
!$omp              shared(MAX_PIP,DES_INTERP_ON,LP_BND,p_force,drag_fc,filter_cell,filter_weight,pijk,pvol)
         DO NP=1,MAX_PIP
            IF(IS_NONEXISTENT(NP).or.IS_ENTERING(NP).or.IS_EXITING(NP).or.IS_ENTERING_GHOST(NP).or.IS_EXITING_GHOST(NP)) CYCLE

            IF(DES_INTERP_ON) THEN
               lPF = ZERO
               DO LC=1,LP_BND
                  IJK = FILTER_CELL(LC,NP)
                  WEIGHT = FILTER_WEIGHT(LC,NP)
                  lPF = lPF + P_FORCE(:,IJK)*WEIGHT
               ENDDO
            ELSE
               lPF = P_FORCE(:,PIJK(NP,4))
            ENDIF

! Include gas pressure and gas-solids drag
            DRAG_FC(:,NP) = DRAG_FC(:,NP) + lPF*PVOL(NP)

         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CALC_PG_GRAD

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_PG_GRAD_STD                                        !
!  Purpose: Calculate cell centered pressure force exerted on the      !
!           particles in the cell by the gas/fluid phase               !
!           Note that P_force is evaluated as -dp/dx                   !
!                                                                      !
!  Notes: This pressure force only needs to be calculated once during  !
!         the DEM loop (at the beginning) since the gas/fluid phase    !
!         is essentially static at that point (i.e., gas field is not  !
!         updated during DEM loop                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PG_GRAD_STD

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE cutcell
      USE fun_avg
      USE functions
      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
! temporary variables used to calculate pressure at scalar cell edge
      DOUBLE PRECISION :: TEMP1, TEMP2
! mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION :: MPG_CYCLIC(3)
!-----------------------------------------------


      MPG_CYCLIC(1:3) = ZERO

      IF(CYCLIC_X_PD) MPG_CYCLIC(1) = DELP_X/XLENGTH
      IF(CYCLIC_Y_PD) MPG_CYCLIC(2) = DELP_Y/YLENGTH
      IF(CYCLIC_Z_PD.AND.DO_K) MPG_CYCLIC(3) = DELP_Z/ZLENGTH

      DO IJK = IJKSTART3, IJKEND3
         P_FORCE(:,IJK) = ZERO
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IMJK = IM_OF(IJK)
         IPJK = IP_OF(IJK)
         IJMK = JM_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKM = KM_OF(IJK)
         IJKP = KP_OF(IJK)

         IF(IMIN1.EQ.IMAX1) THEN
            P_FORCE(1,IJK) = MPG_CYCLIC(1)
         ELSEIF(I.EQ.IMIN1) THEN
            TEMP2 = AVG_X(P_G(IJK), P_G(IPJK), I)
            TEMP1 = P_G(IJK)
            P_FORCE(1,IJK) = 2.d0*(TEMP1-TEMP2)/DX(I)  +  MPG_CYCLIC(1)
         ELSEIF(I.EQ.IMAX1) THEN
            TEMP2 = AVG_X(P_G(IMJK), P_G(IJK), I-1)
            TEMP1 = P_G(IJK)
            P_FORCE(1,IJK) = 2.d0*(TEMP2 - TEMP1)/DX(I) +  MPG_CYCLIC(1)
         ELSEIF((I.GT.IMIN1).AND.(I.LT.IMAX1)) THEN
            TEMP2 = AVG_X(P_G(IJK),  P_G(IPJK), I)
            TEMP1 = AVG_X(P_G(IMJK), P_G(IJK),  I-1)
            P_FORCE(1,IJK) = (TEMP1 - TEMP2)/DX(I) + MPG_CYCLIC(1)
         ENDIF

         IF(JMIN1.EQ.JMAX1) THEN
            P_FORCE(2,IJK) = MPG_CYCLIC(2)
         ELSEIF(J.EQ.JMIN1) THEN
            TEMP2 = AVG_Y(P_G(IJK), P_G(IJPK), J)
            TEMP1 = P_G(IJK)
            P_FORCE(2,IJK) = 2.d0*(TEMP1 - TEMP2)/DY(J)  + MPG_CYCLIC(2)
         ELSEIF(J.EQ.JMAX1) THEN
            TEMP2 = AVG_Y(P_G(IJMK), P_G(IJK), J-1)
            TEMP1 = P_G(IJK)
            P_FORCE(2,IJK) = 2.d0*(TEMP2 - TEMP1)/DY(J) + MPG_CYCLIC(2)
         ELSEIF((J.GT.JMIN1).AND.(J.LT.JMAX1)) THEN
            TEMP2 = AVG_Y(P_G(IJK),  P_G(IJPK), J)
            TEMP1 = AVG_Y(P_G(IJMK), P_G(IJK),  J-1)
            P_FORCE(2,IJK) = (TEMP1 - TEMP2)/DY(J) +  MPG_CYCLIC(2)
         ENDIF

         IF(DO_K) THEN
            IF(KMIN1.EQ.KMAX1) THEN
               P_FORCE(3,IJK) = MPG_CYCLIC(3)
            ELSEIF(K.EQ.KMIN1) THEN
               TEMP2 = AVG_Z(P_G(IJK), P_G(IJKP), K)
               TEMP1 = P_G(IJK)
               P_FORCE(3,IJK) = 2.d0*(TEMP1 - TEMP2)/DZ(K) +  MPG_CYCLIC(3)
            ELSEIF(K.EQ.KMAX1) THEN
               TEMP2 = AVG_Z(P_G(IJKM), P_G(IJK), K-1)
               TEMP1 = P_G(IJK)
               P_FORCE(3,IJK) = 2.d0*(TEMP2 - TEMP1)/DZ(K) +  MPG_CYCLIC(3)
            ELSEIF((K.GT.KMIN1).AND.(K.LT.KMAX1)) THEN
               TEMP2 = AVG_Z(P_G(IJK),  P_G(IJKP), K)
               TEMP1 = AVG_Z(P_G(IJKM), P_G(IJK),  K-1)
               P_FORCE(3,IJK) = (TEMP1 - TEMP2)/DZ(K) +  MPG_CYCLIC(3)
            ENDIF
         ENDIF
      ENDDO         ! end do loop over ijk

      RETURN
      END SUBROUTINE CALC_PG_GRAD_STD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_PG_GRAD_CG                                         !
!  Purpose: Calculate cell centered pressure force exerted on the      !
!           particles in the cell by the gas/fluid phase               !
!           (cut-cell version)                                         !
!                                                                      !
!  Notes: This pressure force needs to be calculated once in a DEM     !
!         time step (at the beginning) since the gas/fluid phase is    !
!         not updated (is static) during the DEM portion of the        !
!         simulation.                                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PG_GRAD_CG

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE fun_avg
      USE functions
      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
! counters
      INTEGER :: X_COUNT, Y_COUNT, Z_COUNT
! mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION :: MPG_CYCLIC(3)
!-----------------------------------------------
      MPG_CYCLIC(1:3) = ZERO

      IF(CYCLIC_X_PD) MPG_CYCLIC(1) = DELP_X/XLENGTH
      IF(CYCLIC_Y_PD) MPG_CYCLIC(2) = DELP_Y/YLENGTH
      IF(CYCLIC_Z_PD.AND.DO_K) MPG_CYCLIC(3) = DELP_Z/ZLENGTH

      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         P_FORCE(:,IJK) = ZERO

         IF(.NOT.FLUID_AT(IJK).OR..NOT.IS_ON_myPE_owns(I, J, K)) CYCLE

         X_COUNT = 0
         Y_COUNT = 0
         Z_COUNT = 0

         IJKE = IP_OF(IJK)
         IJKW = IM_OF(IJK)
         IJKN = JP_OF(IJK)
         IJKS = JM_OF(IJK)
         IJKT = KP_OF(IJK)
         IJKB = KM_OF(IJK)

         IF(FLUID_AT(IJKE)) THEN
            X_COUNT = X_COUNT + 1
            P_FORCE(1,IJK) = P_FORCE(1,IJK) +                          &
               2.d0*(P_G(IJKE) - P_G(IJK))/(DX(I) + DX(I_OF(IJKE)))
         ENDIF
         IF(FLUID_AT(IJKW)) THEN
            X_COUNT = X_COUNT + 1
            P_FORCE(1,IJK) = P_FORCE(1,IJK) +                          &
               2.d0*(P_G(IJK) - P_G(IJKW))/(DX(I) + DX(I_OF(IJKW)))
         ENDIF
         X_COUNT = MAX(1, X_COUNT) !to prevent division from zero
! P_FORCE (by convention) is stored as -dp/dx. MPG_CYCLIC is already -dp/dx.
! therefore, P_force is multiplied by "-" for consistency
         P_FORCE(1,IJK) = MPG_CYCLIC(1) - P_FORCE(1,IJK)/REAL(X_COUNT)

         IF(FLUID_AT(IJKN)) THEN
            Y_COUNT = Y_COUNT + 1
            P_FORCE(2,IJK) = P_FORCE(2,IJK) +                          &
               2.d0*(P_G(IJKN) - P_G(IJK))/(DY(J) + DY(J_OF(IJKN)))
         ENDIF

         IF(FLUID_AT(IJKS)) THEN
            Y_COUNT = Y_COUNT + 1
            P_FORCE(2,IJK) = P_FORCE(2,IJK) +                          &
               2.d0*(P_G(IJK) - P_G(IJKS))/(DY(J) + DY(J_OF(IJKS)))
         ENDIF
         Y_COUNT = MAX(1, Y_COUNT) !to prevent division from zero
         P_FORCE(2,IJK) = MPG_CYCLIC(2) - P_FORCE(2,IJK)/REAL(Y_COUNT)

         IF(DO_K) THEN
            IF(FLUID_AT(IJKT)) THEN
               Z_COUNT = Z_COUNT + 1
               P_FORCE(3,IJK) = P_FORCE(3,IJK) +                       &
                  2.d0*(P_G(IJKT) - P_G(IJK))/(DZ(K) + DZ(K_OF(IJKT)))
            ENDIF
            IF(FLUID_AT(IJKB)) THEN
               Z_COUNT = Z_COUNT + 1
               P_FORCE(3,IJK) = P_FORCE(3,IJK) +                       &
                  2.d0*(P_G(IJK) - P_G(IJKB))/(DZ(K) + DZ(K_OF(IJKB)))
            ENDIF
            Z_COUNT = MAX(1, Z_COUNT) !to prevent division from zero
            P_FORCE(3,IJK) = MPG_CYCLIC(3)-P_FORCE(3,IJK)/REAL(Z_COUNT)
         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE CALC_PG_GRAD_CG
