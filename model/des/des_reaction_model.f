!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_REACTION_MODEL                                     !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_REACTION_MODEL

      USE compar
      Use constant
      Use des_rxns
      Use des_thermo
      Use derived_types, only: pic
      Use discretelement
      USE geometry
      USE indices
      Use param, only: dimension_n_s
      use run, only: ANY_SPECIES_EQ, SPECIES_EQ
      use physprop, only: NMAX
      use run, only: DT
      use run, only: SOLVE_ROs
      use functions

      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! None

! Local variables
!-----------------------------------------------
! Loop counter
      INTEGER :: NN
! total rate of consumption/production of species (g/sec)
      DOUBLE PRECISION  :: SUM_DES_Rs(1:MAX_PIP)

      DOUBLE PRECISION :: PIx4o3
      DOUBLE PRECISION :: o3 = 1.0d0/3.0d0

      DOUBLE PRECISION :: lDT
! Logical for Adams-Bashfort integration.
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

!---------------------------------------------------------------------//

      IF(.NOT.ANY_SPECIES_EQ) RETURN

      PIx4o3 = Pi*4.0d0/3.0d0

      lDT = merge(DT, DTSOLID, DES_EXPLICITLY_COUPLED)

! First-order method: Euler
      IF(INTG_EULER) THEN
         WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)            &
            SUM_DES_Rs(:MAX_PIP) = sum(DES_R_s(:MAX_PIP,:),DIM=2)

         WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)            &
            PMASS(:MAX_PIP) = PMASS(:MAX_PIP) + lDT*SUM_DES_Rs(:MAX_PIP)

         FORALL(NN=1:DIMENSION_N_S)
            WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)         &
               DES_X_s(:MAX_PIP,NN) = max(DES_X_s(:MAX_PIP,NN) + lDT*  &
               (DES_R_s(:MAX_PIP,NN) - DES_X_s(:MAX_PIP,NN)*           &
               SUM_DES_Rs(:MAX_PIP))/PMASS(:MAX_PIP), ZERO)
         END FORALL

      ELSE
         IF(FIRST_PASS) THEN
         ENDIF
      ENDIF

      DO NN=1,MAX_PIP
         IF(IS_NORMAL(NN)) THEN
            IF(SOLVE_ROs(PIJK(NN,5))) THEN
               RO_Sol(NN)= PMASS(NN)/PVOL(NN)
            ELSE
               DES_RADIUS(NN) = (PMASS(NN)/(Pix4o3*RO_SOL(NN)))**o3
               PVOL(NN) = PMASS(NN)/RO_SOL(NN)
            ENDIF
         ENDIF
      ENDDO

! Clear the necessary variables.
      DES_R_s = ZERO

! Flag that the first pass is over
      FIRST_PASS = .FALSE.

      RETURN
      END SUBROUTINE DES_REACTION_MODEL
