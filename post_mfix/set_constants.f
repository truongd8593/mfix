CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: SET_CONSTANTS                                          C
C  Purpose: This module sets all the constants                         C
C                                                                      C
C  Author: M. Syamlal                                 Date: 30-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 31-JAN-92  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: compute constants for plastic_flow stress terms            C
C  Author: M. Syamlal                                 Date: 11-FEB-93  C
C                                                                      C
C  Revision Number: 2                                                  C
C  Purpose: Add K_scale                                                C
C  Author: W. Sams                                    Date: 03-MAY-93  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: UNITS                                         C
C                                                                      C
C  Variables modified: G, GAS_CONST, K_scale, Pi, SQRT_Pi, SQRT_3,     C
C                      ETA, D_p3, oD_p3, MASS_s                        C
C                                                                      C
C  Local variables: IJK                                                 C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE SET_CONSTANTS
C
      IMPLICIT NONE
C
C             Local index
      INTEGER IJK, M
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'visc_s.inc'
      INCLUDE 'energy.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'run.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'drag.inc'
C
C   Note that the cell flags are not when this routine is called.


C
C  Dimensionless constants
C
      Pi  = 4. * ATAN(ONE)
      SQRT_Pi = SQRT(Pi)
      K_scale = .08
      EP_s_cp = 1. - EP_star
      Eta = (1d0 + C_e)*0.5d0
C
C  plastic regime stress
C
      tan_Phi_w  = TAN(PHI_w*Pi/180.D0)
      Sin_Phi  = SIN(PHI*Pi/180.D0)
      Sin2_Phi = Sin_Phi * Sin_Phi
      F_Phi    = (3.0 - 2.0 * Sin2_Phi)/ 3.0
C
C  Enter the value of all constants in various units
C
      IF( UNITS .EQ. 'SI' ) THEN
        IF(GRAVITY .EQ. UNDEFINED)
     &    GRAVITY   =  9.80665        ! m/s^2         (Perry and Green, 1984)
        GAS_CONST = 8314.56         ! Pa.m^3/kmol.K (Perry and Green, 1984)
      ELSEIF( UNITS .EQ. 'CGS' ) THEN
        IF(GRAVITY .EQ. UNDEFINED)
     &    GRAVITY   = 980.665
        GAS_CONST = 8314.56E4
      ELSE
        WRITE(UNIT_LOG, 1000) UNITS
        CALL EXIT
      ENDIF

      RETURN
1000  FORMAT(/70('*')//'From: SET_CONSTANTS'/'Message: Unknown UNITS: ',
     & 1A16,/70('*')/)
      END
