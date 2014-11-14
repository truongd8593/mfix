!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_CONSTANTS                                           C
!  Purpose: This module sets all the constants                         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 31-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: compute constants for plastic_flow stress terms            C
!  Author: M. Syamlal                                 Date: 11-FEB-93  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Add K_scale                                                C
!  Author: W. Sams                                    Date: 03-MAY-93  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: Add to_SI to change from CGS to SI in some routines        C
!  Author: S. Dartevelle                              Date: 03-MAY-02  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: UNITS                                         C
!  Variables modified: GRAVITY, GAS_CONST, K_scale, Pi, SQRT_Pi,       C
!                      ETA, tan_phi_w, sin_phi, sin2_phi, f_phi        C
!                      lam_hys                                         C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_CONSTANTS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE visc_s
      USE energy
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE run
      USE funits
      USE drag
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------

! Note that the cell flags are not set when this routine is called.

! Dimensionless constants
      K_SCALE = .08D0   ! this actually isn't used anywhere...

! Enter the value of all constants in various units (CGS or SI)
      IF (UNITS == 'SI') THEN
         IF (GRAVITY == UNDEFINED) GRAVITY = 9.80665D0 ! m/s2
         GAS_CONST = 8314.56D0                     !Pa.m3/kmol.K, or kg m2/s2 kmol K (Perry and Green, 1984)
         to_SI = 0.1D0                             !to convert dyne/cm2 to Pa, see s_pr2.inc, see calc_mu_g.f
         IF (LAM_HYS == UNDEFINED) LAM_HYS = 0.000001d0    ! m
      ELSEIF (UNITS == 'CGS') THEN
         IF (GRAVITY == UNDEFINED) GRAVITY = 980.665D0 !cm/s2
         GAS_CONST = 8314.56D4                   !g.cm2/s2.mol.K
         to_SI = ONE                             !does not do anything in CGS,  see s_pr2.inc, see calc_mu_g.f
         IF (LAM_HYS == UNDEFINED) LAM_HYS = 0.0001d0    ! cm
      ELSE
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) UNITS
         CALL MFIX_EXIT(myPE)
      ENDIF

! If the gravitational acceleration vector is undefined, then used default value in negative y-direction
! This ensured backward compatibility with the old (legacy) GRAVITY keyword.
! At this point GRAVITY is defined,either from mfix.dat or by default above
      IF(GRAVITY_X==ZERO.AND.GRAVITY_Y==ZERO.AND.GRAVITY_Z==ZERO) THEN
         GRAVITY_X = ZERO
         GRAVITY_Y = - GRAVITY
         GRAVITY_Z = ZERO
      ENDIF



      RETURN
 1000 FORMAT(/1X,70('*')//'From: SET_CONSTANTS',/&
         ' Message: Unknown UNITS: ',1A16,/1X,70('*')/)

      END SUBROUTINE SET_CONSTANTS
