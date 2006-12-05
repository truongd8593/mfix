!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: constants.inc                                          C
!  Purpose: Common block containing physical constants and constants   C
!           used in the numerical technique                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 5-FEB-92   C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Add parameters MIN_EP_S and DEPSTAR.                       C
!  Author: M. Syamlal                                 Date: 7-FEB-92   C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      MODULE constant


      Use param
      Use param1

!       multiple particle sizes
        DOUBLE PRECISION ep_s_max(DIM_M)  !maximum packing volume fraction for particles, typically 
	DOUBLE PRECISION ep_s_max_ratio(DIM_M, DIM_M), d_p_ratio(DIM_M, DIM_M)
!
! commented by sof(05-04-2005) no need for this, will be user input in mfix.dat
	DOUBLE PRECISION SEGREGATION_SLOPE_COEFFICIENT !, MAX_SOLID_1_PACKING,&
			  !MAX_SOLID_2_PACKING0.6

! success-factor for aggregation and breakage  

        DOUBLE PRECISION AGGREGATION_EFF
        DOUBLE PRECISION BREAKAGE_EFF
!     ALPHA = parameter in equation for mu_s
!
!     SWITCH enables us to turn on/off the modification to the
!     particulate phase viscosity. If we want to simulate gas-particle
!     flow then SWITCH=1 to incorporate the effect of drag on the
!     particle viscosity. If we want to simulate granular flow
!     without the effects of an interstitial gas, SWITCH=0.
!     (Same for conductivity)
      DOUBLE PRECISION ALPHA, SWITCH
      PARAMETER(ALPHA = 1.6d0, SWITCH=1d0)
!
!     SWITCH_IA changes some terms in solids viscosity and conductivity
!     in order for the results using 2 or more identical solids phase to 
!     be the same as 1 soldis phase. Set to false to use original theory
!     of Iddir-Arastoopour. Sof DEC 05 2006.
      LOGICAL, PARAMETER :: SWITCH_IA = .TRUE.
 
 
!		       PHIP  = Specularity coefficient associated with
!                              particle wall collisions
      DOUBLE PRECISION PHIP
 
!	               e_w   = particle-wall coefficient of restitution	
      DOUBLE PRECISION e_w
 
!     Fr, EPS_f_min, N_Pc, D_Pc all appear in the equation for Pc,the
!     critical solids pressure. N_Pf appears as an exponent in the
!     equation of state for Pf (frictional pressure)
!	               Fr = Constant
!                      N_Pc = exponent in numerator
!	               D_Pc = exponent in denominator
 
      DOUBLE PRECISION Fr, N_Pc, D_Pc, N_Pf
      PARAMETER(Fr = 0.5d0, N_Pc=2d0, D_Pc=5d0, N_Pf=1.03d0)
 
!	               EPS_f_min = minimum solids fraction above which
!                                  friction kicks in
      DOUBLE PRECISION EPS_f_min

! not needed anymore (sof, Nov-17-2005) 
!		       EPS_max = random close-packed solids volume fraction
!      DOUBLE PRECISION EPS_max		
!
!                      Gravitational acceleration
      DOUBLE PRECISION GRAVITY, to_SI
!
!                      Universal gas constant
      DOUBLE PRECISION GAS_CONST
!
!                      Universal gas constant in cal/mol.K
      DOUBLE PRECISION, PARAMETER :: GAS_CONST_cal = 1.987207D0
!
!                      Coeficient of restitution
      DOUBLE PRECISION C_e
!
!                      (1+C_e)/2.
      DOUBLE PRECISION eta
!
!                      Coeficient of friction
      DOUBLE PRECISION C_f
!
!                      Packed bed (close packed) void fraction
      DOUBLE PRECISION EP_star
!
!                      Packed bed  solids fraction (1 - EP_star)
      DOUBLE PRECISION EP_s_cp
!
!                      Pi, the ubiquitous irrational number
      DOUBLE PRECISION Pi
!
!                      Square root of Pi
      DOUBLE PRECISION SQRT_Pi
!
!                      Maximum pressure correction allowed in one iteration
      DOUBLE PRECISION MAX_DELP
!
!                      User defined constants
      DOUBLE PRECISION C (DIMENSION_C)
!
!                      Names of user defined constants (for output file only)
      CHARACTER*20     C_NAME (DIMENSION_C)
!
!                      Angle of internal friction (degrees)
      DOUBLE PRECISION Phi
!
!                      Angle of wall-particle friction (degrees)
      DOUBLE PRECISION Phi_w
!
!                      (k=) Sin(PHI) in Plastic-flow stress formulation
      DOUBLE PRECISION Sin_Phi
!
!                      Sin^2(PHI)
      DOUBLE PRECISION Sin2_Phi
!
!                      (3-2k^2)/6k^2 in Plastic-flow stress formulation
      DOUBLE PRECISION F_Phi
!
!                      tan(PHI_w)
      DOUBLE PRECISION tan_Phi_w
 
!
!                      Default value for characteristic length
!                      for turbulence
      DOUBLE PRECISION L_scale0
!
!                      Maximum value of turbulent viscosity
      DOUBLE PRECISION MU_gmax
!
!                      Excluded volume (Boyle-Massoudi stress tensor)
      DOUBLE PRECISION V_ex
!
!                      Coefficients for calibrating Syamlal-O'Brien drag correlation with Umf data
      DOUBLE PRECISION drag_c1, drag_d1
!


      END MODULE constant                                                                        
