!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_DRAG                                               !
!                                                                      !
!  Purpose: Provide a hook for user defined drag law implementation.   !
!                                                                      !
!  This routine is called from inside fluid (TFM) and particle (DES)   !
!  loops. The fluid cell index (IJK) and phase (TFM) or particle index !
!  (DES) is passed.                                                    !
!                                                                      !
!  ***************************   WARNING   **************************  !
!  *----------------------------------------------------------------*  !
!  * The dummy arguments changed in the 2015-1 MFIX Release.        *  !
!  *                                                                *  !
!  *   1) Phase index (M) is now particle index (NP) for DES. This  *  !
!  *      is reflected in the name change M --> M_NP.               *  !
!  *                                                                *  !
!  *   2) The fluid velocity was added as a dummy argument. This    *  !
!  *      provides access to the interpolated gas velocity for      *  !
!  *      coupled DES simulations.                                  *  !
!  *                                                                *  !
!  ******************************************************************  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_USR(IJK, M_NP, lDgA, EPg, Mug, ROg, VREL, DPM, &
         ROs, lUg, lVg, lWg)

      use constant, only: PSI_s => C
      use constant, only: GRAVITY
! Fluid cell I, J, K, IJK containing particle and phase
      use discretelement, only: PIJK

      use error_manager

      IMPLICIT NONE

! Index of fluid cell:
      INTEGER, INTENT(IN) :: IJK
! TFM SOLIDS --> Index of phase (M)
! DES SOLIDS --> Index of particle (NP); M = PIJK(NP,5)
      INTEGER, INTENT(IN) :: M_NP

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
! fluid velocity components:
! o TFM: Averaged from faces to cell center
! o DES: Interpolated to the particle's position
      DOUBLE PRECISION, INTENT(IN) :: lUg, lVg, lWg

! Drag correlation.
      DOUBLE PRECISION :: C_d

! Phase index
      INTEGER :: M

! Get the particle phase index.
      M = PIJK(M_NP,5)

      C_d = (0.69d0*GRAVITY*(DPM**3)*ROg*1.33*(ROs-ROg))/    &
         (Mug**2*((GRAVITY*(PSI_s(M)**1.6)*(DPM**3)*ROg*(ROs-ROg)) / &
         Mug**2)**1.0412)

      lDgA = 0.75 * C_d * VREL * (EPg * ROg) * EPg**(-2.65) / DPM

      RETURN
      END SUBROUTINE DRAG_USR
