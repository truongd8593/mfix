!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_CONSTANTS                                          C
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
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: UNITS                                         C
!                                                                      C
!  Variables modified: G, GAS_CONST, K_scale, Pi, SQRT_Pi, SQRT_3,     C
!                      ETA, D_p3, oD_p3, MASS_s                        C
!                                                                      C
!  Local variables: IJK
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_CONSTANTS 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
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
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJK, M 
!-----------------------------------------------
!
!
!
!   Note that the cell flags are not set when this routine is called.
!
!
!
!  Dimensionless constants
!
      PI = 4.*ATAN(ONE) 
      SQRT_PI = SQRT(PI) 
      K_SCALE = .08 
      EP_S_CP = 1. - EP_STAR 
      ETA = (1D0 + C_E)*0.5D0 
!
!  plastic regime stress
!
      if(mmax > 0) then
        TAN_PHI_W = TAN(PHI_W*PI/180.D0) 
        SIN_PHI = SIN(PHI*PI/180.D0) 
        SIN2_PHI = SIN_PHI*SIN_PHI 
        F_PHI = (3.0 - 2.0*SIN2_PHI)/3.0 
      endif
!
!  Enter the value of all constants in various units
!
      IF (UNITS == 'SI') THEN 
!                                                ! m/s^2         (Perry and Green,
         IF (GRAVITY == UNDEFINED) GRAVITY = 9.80665 
         GAS_CONST = 8314.56                     ! Pa.m^3/kmol.K (Perry and Green, 1984) 
      ELSE IF (UNITS == 'CGS') THEN 
         IF (GRAVITY == UNDEFINED) GRAVITY = 980.665 
         GAS_CONST = 8314.56E4 
      ELSE 
         WRITE (UNIT_LOG, 1000) UNITS 
         CALL MFIX_EXIT 
      ENDIF 
!
      RETURN  
 1000 FORMAT(/70('*')//'From: SET_CONSTANTS'/'Message: Unknown UNITS: ',1A16,/&
         70('*')/) 
      END SUBROUTINE SET_CONSTANTS 
