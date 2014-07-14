!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_THERMO_DES                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_THERMO_DES

      use physprop, only: SMAX
      use physprop, only: K_s0

      USE compar
      Use des_rxns
      Use des_thermo
      Use discretelement
      Use fldvar
      USE geometry
      USE indices
      Use interpolation
      Use param1
      Use run

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! Index of neighbor particle of particle I such that I < J
      INTEGER I, J, K, IJK
! Loop index for particles.
      INTEGER NP, lNP
! Phase index for particle NP
      INTEGER M
! Identifies that the indicated particle is of interest for debugging
      LOGICAL FOCUS
! Variables needed for calculating new interpolation quantities for
! species and energy equations
      INTEGER INTERP_IJK(2**3)
      DOUBLE PRECISION INTERP_WEIGHTS(2**3)
! Flag to calculate gas/particle convective heat transfer.
      LOGICAL :: CALC_CONV
! Flag to calculate radiative heat transfer
      LOGICAL :: CALC_RADT(DIM_M)
! Flag to calculate conductive heat transfer
      LOGICAL :: CALC_COND(DIM_M)

! Functions
!---------------------------------------------------------------------//
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      INCLUDE '../function.inc'


! This is a quick work-around to keep the thermo routines from causes
! issues while the "check_data" routines are rewritten. Moving forward
! this routine should be split apart to avoid the particle loops for
! cold-flow, non-reacting cases.
      IF(.NOT.ENERGY_EQ .AND. .NOT.ANY_SPECIES_EQ) RETURN

      CALC_CONV = .FALSE.
      CALC_RADT = .FALSE.
      CALC_COND = .FALSE.

! Set flags for energy equations:
      IF(ENERGY_EQ) THEN
! Flag to calculate convection.
         CALC_CONV = DES_CONTINUUM_COUPLED
         DO M=1, SMAX + DES_MMAX
! Flag to calculate radiation.
            IF(DES_Em(M) > ZERO) CALC_RADT(M) = .TRUE.
! Flag to calculate conduction.
            IF(K_s0(M) > ZERO) CALC_COND(M) = .TRUE.
         ENDDO
      ENDIF

! Loop over fluid cells.
!---------------------------------------------------------------------//
      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE IJK_LP
         IF(PINC(IJK) == 0) CYCLE IJK_LP

! Interpolation: Removed J.Musser 11/8/2012
!---------------------------------------------------------------------//
!     IF(DES_INTERP_ON .AND. (ANY_SPECIES_EQ .OR. DES_CONV_EQ)) THEN
!         INTERP_IJK(:) = -1
!         INTERP_WEIGHTS(:) = ZERO
!         CALL INTERPOLATE_CC(NP, INTERP_IJK, INTERP_WEIGHTS, FOCUS)
!      ENDIF

! Preform user-defined calculations from fluid grid.
         IF(CALL_USR) CALL USR4_DES(IJK)

! Loop over all particles in cell IJK.
!---------------------------------------------------------------------//
         lNP_LP: DO lNP = 1, PINC(IJK)
            NP = PIC(IJK)%p(lNP)

! Skip indices that do not represent particles
            IF(.NOT.PEA(NP,1)) CYCLE lNP_LP

! Skip indices that represent ghost particles
            IF(PEA(NP,4)) CYCLE lNP_LP

! Reset the debug flag
            FOCUS = .FALSE.

! Calculate time dependent physical properties
            CALL DES_PHYSICAL_PROP(NP, FOCUS)

! Identify the solid phases of each particle
            M = PIJK(NP,5)

! calculate heat transfer via convection
            IF(CALC_CONV) CALL DES_CONVECTION(NP, M, IJK, &
               INTERP_IJK, INTERP_WEIGHTS, FOCUS)

! calculate heat transfer via radiation
            IF(CALC_RADT(M)) CALL DES_RADIATION(NP, M, IJK, FOCUS)

! Loop over thermodynamic neighbor for conduction and radiation
            IF(CALC_COND(M)) CALL DES_CONDUCTION(NP, M, IJK, FOCUS)

! Calculate reaction rates and interphase mass transfer
            IF(ANY_SPECIES_EQ) CALL DES_RRATES0(NP, M, IJK, &
               INTERP_IJK, INTERP_WEIGHTS, FOCUS)

         ENDDO lNP_LP ! End loop over all particles
      ENDDO IJK_LP ! End loop over fluid cells

      END SUBROUTINE CALC_THERMO_DES
