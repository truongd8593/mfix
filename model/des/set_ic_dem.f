
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE Name: DES_SET_IC                                         !
!                                                                      !
!  Purpose: Assign initial conditions to particles basded upon thier   !
!  location within the domain.                                         !
!                                                                      !
!  Author: J.Musser                                   Date: 15-Feb-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_IC_DEM


      use run, only: ENERGY_EQ, SPECIES_EQ
      use run, only: RUN_TYPE

      use ic

      use des_thermo, only: DES_T_s_NEW

      use discretelement, only: PINC, PIC
      use discretelement, only: PIJK

      USE des_rxns, only: DES_X_s

      use physprop, only: NMAX

      USE compar
      use indices
      use geometry


      IMPLICIT NONE

! Dummy indices
      INTEGER :: ICV
      INTEGER :: I, J, K, IJK
      INTEGER :: M, N
      INTEGER :: NP
      INTEGER :: NINDX

      LOGICAL , EXTERNAL :: COMPARE 


      include "../function.inc"


      IF(RUN_TYPE /= 'NEW') RETURN

      DO ICV = 1, DIMENSION_IC
         IF(.NOT.IC_DEFINED(ICV)) CYCLE

         DO K = IC_K_B(ICV), IC_K_T(ICV)
         DO J = IC_J_S(ICV), IC_J_N(ICV)
         DO I = IC_I_W(ICV), IC_I_E(ICV)

! Set the initial conditions for particles in cells that are
! not dead and that this rank owns.
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IF (DEAD_CELL_AT(I,J,K)) CYCLE

            IJK = FUNIJK(I,J,K)

! Loop through particles in cell IJK.
            DO NINDX = 1,PINC(IJK)
               NP = PIC(IJK)%P(NINDX)
               M = PIJK(NP,5)

! Set the initial particle temperature.
               IF(ENERGY_EQ) THEN
                  DES_T_s_NEW(NP) = IC_T_s(ICV,M)
               ENDIF

! Set the initial species composition.
               IF(SPECIES_EQ(M)) THEN
                  DES_X_s(NP,:) = ZERO
                  DO N = 1, NMAX(M)
                     write(*,*)'IC_X_s:', IC_X_s(ICV,M,N)
                     DES_X_s(NP,N) = IC_X_s(ICV,M,N)
                  ENDDO
               ENDIF
            ENDDO

         ENDDO
         ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE SET_IC_DEM
