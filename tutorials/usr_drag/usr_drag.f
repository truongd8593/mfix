!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_DRAG                                               !
!                                                                      !
!  Purpose: Provide a hook for user defined drag law implementation.   !
!                                                                      !
!  This routine is called from inside a fluid-loop and passes the      !
!  fluid cell index, and the index of the phase for which the drag     !
!  force is being calculated.                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_USR(IJK, MNP, lDgA, EPg, Mug, ROg, VREL, DPM, ROs)

      use constant, only: PSI_s => C
      use constant, only: GRAVITY
! Fluid cell I, J, K, IJK containing particle and phase
      use discretelement, only: PIJK

      use error_manager

      IMPLICIT NONE

! Index of fluid cell:
      INTEGER, INTENT(IN) :: IJK
! TFM SOLIDS --> Index of phase (M)
! DES SOLIDS --> Index of particle (NP); M = PIJK(NP,4)
      INTEGER, INTENT(IN) :: MNP

! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
! particle density of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: ROs

! Drag correlation.
      DOUBLE PRECISION :: C_d

! Phase index
      INTEGER :: M

! Get the particle phase index.
      M = PIJK(MNP,4)

      C_d = (0.69d0*GRAVITY*(DPM**3)*ROg*1.33*(ROs-ROg))/    &
         (Mug**2*((GRAVITY*(PSI_s(M)**1.6)*(DPM**3)*ROg*(ROs-ROg)) / &
         Mug**2)**1.0412)

      lDgA = 0.75 * C_d * VREL * (EPg * ROg) * EPg**(-2.65) / DPM

      RETURN
      END SUBROUTINE DRAG_USR
