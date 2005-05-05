!TO DO
!  define ep_s_max as an input variable... Done (sof).
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
!  Revision Number: 3                                                  C
!  Purpose: Add to_SI to change from CGS to SI in some routines        C
!  Author: S. Dartevelle                              Date: 03-MAY-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: UNITS                                         C
!                                                                      C
!  Variables modified: G, GAS_CONST, K_scale, Pi, SQRT_Pi, SQRT_3,     C
!                      ETA, D_p3, oD_p3, MASS_s                        C
!                                                                      C
!  Local variables: IJK                                                C
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
      USE compar
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
!       For multiple particle types
! commented by sof (05-04-2005) ep_s_max(MMAX) will be defined in mfix.dat
!        ep_s_max(1) = MAX_SOLID_1_PACKING  ! maximum packing volume fraction for spheres, typically 0.6
!	ep_s_max(2) = MAX_SOLID_2_PACKING ! maximum packing volume fraction for solids, typically 0.6

	if (d_p(2) .GT. d_p(1)) then
	  d_p_ratio(1,2) = (d_p(1)/d_p(2))**0.5
	else
	  d_p_ratio(1,2) = (d_p(2)/d_p(1))**0.5
	end if
	ep_s_max_ratio(1,2) = ep_s_max(1)/(ep_s_max(1)+(1.-ep_s_max(1))*ep_s_max(2))  ! refer to Syam's dissertation
	
!
!  Dimensionless constants
      PI = 4.*ATAN(ONE) 
      SQRT_PI = SQRT(PI) 
      K_SCALE = .08 
      EP_S_CP = 1. - EP_STAR 
      ETA = (1D0 + C_E)*0.5D0 
!
!  plastic regime stress
!Angle given in degree but calculated in radian within the fortran codes
      if(mmax > 0) then
        TAN_PHI_W = TAN(PHI_W*PI/180.D0) 
        SIN_PHI = SIN(PHI*PI/180.D0) 
        SIN2_PHI = SIN_PHI*SIN_PHI 
        F_PHI = (3.0 - 2.0*SIN2_PHI)/3.0 
      endif
!
!  Enter the value of all constants in various units (CGS or SI)
      IF (UNITS == 'SI') THEN
         IF (GRAVITY == UNDEFINED) GRAVITY = 9.80665 ! m/s2
         GAS_CONST = 8314.56                     !Pa.m3/kmol.K, or kg m2/s2 kmol K (Perry and Green, 1984)
         to_SI = 0.1                             !to convert dyne/cm2 to Pa, see s_pr2.inc, see calc_mu_g.f
      ELSE IF (UNITS == 'CGS') THEN
         IF (GRAVITY == UNDEFINED) GRAVITY = 980.665 !cm/s2
         GAS_CONST = 8314.56E4                   !g.cm2/s2.mol.K
         to_SI = ONE                             !does not do anything in CGS,  see s_pr2.inc, see calc_mu_g.f
      ELSE 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) UNITS
         CALL MFIX_EXIT(myPE) 
      ENDIF 
!
      RETURN  
 1000 FORMAT(/70('*')//'From: SET_CONSTANTS'/'Message: Unknown UNITS: ',1A16,/&
         70('*')/) 
      END SUBROUTINE SET_CONSTANTS 
