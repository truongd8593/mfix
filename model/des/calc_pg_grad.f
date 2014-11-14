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

      implicit none

      IF(CARTESIAN_GRID) THEN
         CALL CALC_PG_GRAD_CG
      ELSE
         CALL CALC_PG_GRAD_STD
      ENDIF

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
         P_FORCE(IJK, :) = ZERO
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
            P_FORCE(IJK,1) = MPG_CYCLIC(1)
         ELSEIF(I.EQ.IMIN1) THEN
            TEMP2 = AVG_X(P_G(IJK), P_G(IPJK), I)
            TEMP1 = P_G(IJK)
            P_FORCE(IJK,1) = 2.d0*(TEMP1-TEMP2)/DX(I)  +  MPG_CYCLIC(1)
         ELSEIF(I.EQ.IMAX1) THEN
            TEMP2 = AVG_X(P_G(IMJK), P_G(IJK), I-1)
            TEMP1 = P_G(IJK)
            P_FORCE(IJK,1) = 2.d0*(TEMP2 - TEMP1)/DX(I) +  MPG_CYCLIC(1)
         ELSEIF((I.GT.IMIN1).AND.(I.LT.IMAX1)) THEN
            TEMP2 = AVG_X(P_G(IJK),  P_G(IPJK), I)
            TEMP1 = AVG_X(P_G(IMJK), P_G(IJK),  I-1)
            P_FORCE(IJK,1) = (TEMP1 - TEMP2)/DX(I) + MPG_CYCLIC(1)
         ENDIF

         IF(JMIN1.EQ.JMAX1) THEN
            P_FORCE(IJK,2) = MPG_CYCLIC(2)
         ELSEIF(J.EQ.JMIN1) THEN
            TEMP2 = AVG_Y(P_G(IJK), P_G(IJPK), J)
            TEMP1 = P_G(IJK)
            P_FORCE(IJK,2) = 2.d0*(TEMP1 - TEMP2)/DY(J)  + MPG_CYCLIC(2)
         ELSEIF(J.EQ.JMAX1) THEN
            TEMP2 = AVG_Y(P_G(IJMK), P_G(IJK), J-1)
            TEMP1 = P_G(IJK)
            P_FORCE(IJK,2) = 2.d0*(TEMP2 - TEMP1)/DY(J) + MPG_CYCLIC(2)
         ELSEIF((J.GT.JMIN1).AND.(J.LT.JMAX1)) THEN
            TEMP2 = AVG_Y(P_G(IJK),  P_G(IJPK), J)
            TEMP1 = AVG_Y(P_G(IJMK), P_G(IJK),  J-1)
            P_FORCE(IJK,2) = (TEMP1 - TEMP2)/DY(J) +  MPG_CYCLIC(2)
         ENDIF

         IF(DO_K) THEN
            IF(KMIN1.EQ.KMAX1) THEN
               P_FORCE(IJK,3) = MPG_CYCLIC(3)
            ELSEIF(K.EQ.KMIN1) THEN
               TEMP2 = AVG_Z(P_G(IJK), P_G(IJKP), K)
               TEMP1 = P_G(IJK)
               P_FORCE(IJK,3) = 2.d0*(TEMP1 - TEMP2)/DZ(K) +  MPG_CYCLIC(3)
            ELSEIF(K.EQ.KMAX1) THEN
               TEMP2 = AVG_Z(P_G(IJKM), P_G(IJK), K-1)
               TEMP1 = P_G(IJK)
               P_FORCE(IJK,3) = 2.d0*(TEMP2 - TEMP1)/DZ(K) +  MPG_CYCLIC(3)
            ELSEIF((K.GT.KMIN1).AND.(K.LT.KMAX1)) THEN
               TEMP2 = AVG_Z(P_G(IJK),  P_G(IJKP), K)
               TEMP1 = AVG_Z(P_G(IJKM), P_G(IJK),  K-1)
               P_FORCE(IJK,3) = (TEMP1 - TEMP2)/DZ(K) +  MPG_CYCLIC(3)
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
!         time step (at the beggining) since the gas/fluid phase is    !
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
         P_FORCE(IJK, :) = ZERO

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
            P_FORCE(IJK,1) = P_FORCE(IJK,1) +                          &
               2.d0*(P_G(IJKE) - P_G(IJK))/(DX(I) + DX(I_OF(IJKE)))
         ENDIF
         IF(FLUID_AT(IJKW)) THEN
            X_COUNT = X_COUNT + 1
            P_FORCE(IJK, 1) = P_FORCE(IJK, 1) +                        &
               2.d0*(P_G(IJK) - P_G(IJKW))/(DX(I) + DX(I_OF(IJKW)))
         ENDIF
         X_COUNT = MAX(1, X_COUNT) !to prevent division from zero
! P_FORCE (by convention) is stored as -dp/dx. MPG_CYCLIC is already -dp/dx.
! therefore, P_force is multiplied by "-" for consistency
         P_FORCE(IJK, 1) = MPG_CYCLIC(1) - P_FORCE(IJK,1)/REAL(X_COUNT)

         IF(FLUID_AT(IJKN)) THEN
            Y_COUNT = Y_COUNT + 1
            P_FORCE(IJK, 2) = P_FORCE(IJK, 2) +                        &
               2.d0*(P_G(IJKN) - P_G(IJK))/(DY(J) + DY(J_OF(IJKN)))
         ENDIF

         IF(FLUID_AT(IJKS)) THEN
            Y_COUNT = Y_COUNT + 1
            P_FORCE(IJK, 2) = P_FORCE(IJK, 2) +                        &
               2.d0*(P_G(IJK) - P_G(IJKS))/(DY(J) + DY(J_OF(IJKS)))
         ENDIF
         Y_COUNT = MAX(1, Y_COUNT) !to prevent division from zero
         P_FORCE(IJK, 2) = MPG_CYCLIC(2) - P_FORCE(IJK,2)/REAL(Y_COUNT)

         IF(DO_K) THEN
            IF(FLUID_AT(IJKT)) THEN
               Z_COUNT = Z_COUNT + 1
               P_FORCE(IJK, 3) = P_FORCE(IJK, 3) +                     &
                  2.d0*(P_G(IJKT) - P_G(IJK))/(DZ(K) + DZ(K_OF(IJKT)))
            ENDIF
            IF(FLUID_AT(IJKB)) THEN
               Z_COUNT = Z_COUNT + 1
               P_FORCE(IJK, 3) = P_FORCE(IJK, 3) +                     &
                  2.d0*(P_G(IJK) - P_G(IJKB))/(DZ(K) + DZ(K_OF(IJKB)))
            ENDIF
            Z_COUNT = MAX(1, Z_COUNT) !to prevent division from zero
            P_FORCE(IJK, 3) = MPG_CYCLIC(3)-P_FORCE(IJK,3)/REAL(Z_COUNT)
         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE CALC_PG_GRAD_CG
