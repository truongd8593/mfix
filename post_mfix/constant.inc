CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: constants.inc                                          C
C  Purpose: Common block containing physical constants and constants   C
C           used in the numerical technique                            C
C                                                                      C
C  Author: M. Syamlal                                 Date: 5-FEB-92   C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: Add parameters MIN_EP_S and DEPSTAR.                       C
C  Author: M. Syamlal                                 Date: 7-FEB-92   C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References: None                                C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C     SWITCH enables us to turn on/off the modification to the
C     particulate phase viscosity. If we want to simulate gas-particle
C     flow then SWITCH=1 to incorporate the effect of drag on the
C     particle viscosity. If we want to simulate granular flow
C     without the effects of an interstitial gas, SWITCH=0.
C     (Same for conductivity)
      DOUBLE PRECISION ALPHA, SWITCH
      PARAMETER(ALPHA = 1.6d0, SWITCH=1d0)
C
C                      Gravitational acceleration
      DOUBLE PRECISION GRAVITY
C
C                      Universal gas constant
      DOUBLE PRECISION GAS_CONST
C
C                      Coeficient of restitution
      DOUBLE PRECISION C_e
C
C                      (1+C_e)/2.
      DOUBLE PRECISION eta
C
C                      Coeficient of friction
      DOUBLE PRECISION C_f
C
C                      Packed bed (close packed) void fraction
      DOUBLE PRECISION EP_star
C
C                      Packed bed  solids fraction (1 - EP_star)
      DOUBLE PRECISION EP_s_cp
C
C                      Pi, the ubiquitous irrational number
      DOUBLE PRECISION Pi
C
C                      Square root of Pi
      DOUBLE PRECISION SQRT_Pi
C
C                      Maximum pressure correction allowed in one iteration
      DOUBLE PRECISION MAX_DELP
C
C                      User defined constants
      DOUBLE PRECISION C (DIMENSION_C)
C
C                      Names of user defined constants (for output file only)
      CHARACTER*20     C_NAME (DIMENSION_C)
C
C                      Angle of internal friction (degrees)
      DOUBLE PRECISION Phi
C
C                      Angle of wall-particle friction (degrees)
      DOUBLE PRECISION Phi_w
C
C                      (k=) Sin(PHI) in Plastic-flow stress formulation
      DOUBLE PRECISION Sin_Phi
C
C                      Sin^2(PHI)
      DOUBLE PRECISION Sin2_Phi
C
C                      (3-2k^2)/6k^2 in Plastic-flow stress formulation
      DOUBLE PRECISION F_Phi
C
C                      tan(PHI_w)
      DOUBLE PRECISION tan_Phi_w

C
C                      Default value for characteristic length
C                      for turbulence
      DOUBLE PRECISION L_scale0
C
C                      Maximum value of turbulent viscosity
      DOUBLE PRECISION MU_gmax
C
C                      Excluded volume (Boyle-Massoudi stress tensor)
      DOUBLE PRECISION V_ex
C
      COMMON /CONSTANTS_DP/
     &                  GRAVITY, GAS_CONST, C_e, C_f, EP_star, Pi,
     &                  SQRT_Pi, MAX_DELP, C, L_scale0, MU_gmax,
     &                  Phi, Phi_w, Sin_Phi, Sin2_Phi, F_Phi, tan_Phi_w,
     &                  V_ex, EP_s_cp, Eta
C
      COMMON /CONSTANTS_C/
     &                    C_NAME
