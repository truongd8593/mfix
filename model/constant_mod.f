!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: constant                                                    C
!  Purpose: Common block containing physical constants and constants   C
!           used in the numerical technique                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 5-FEB-92   C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!    Gera, D., Syamlal, M., and O'Brien, T. J., "Hydrodynamics of      C
!      particle segregation in fluidized beds", Int. J. of Multiphase  C
!      Flow, Vol 30, 2004, pp. 419-428.                                C
!    Johnson, P. C., and Jackson, R., "Frictional-collisional          C
!      constitutive relations for granluar materials, with application C
!      to plane shearing", JFM, Vol. 176, 1987, pp. 67-93.             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE constant


! Modules
!---------------------------------------------------------------------//
      Use param, only: dim_m, dimension_c
!---------------------------------------------------------------------//

! Packed bed (close packed) void fraction
      DOUBLE PRECISION :: EP_star

! maximum packing volume fraction for indicate particulate phase
! its value will default to 1-ep_star
      DOUBLE PRECISION :: ep_s_max(DIM_M)
! Index to rearrange particles from coarsest to finest for use in
! function CALC_ep_star(IJK,IER)
      INTEGER :: M_MAX(DIM_M)

! SWITCH enables us to turn on/off modifications to certain kinetic
! theory models for granular solids (i.e. no gas) that have been
! adjusted to account for the presence of a fluid phase. If one wants
! to simulate gas-particle flow then set SWITCH=1. As a result, the
! effects of drag on particle viscosity/conductivity will be
! incorporated. Additional gas-solids terms may also have been
! introduced into the granular energy balance depending on the KT
! model (see source_granular_energy for details). If we want to
! simulate pure granular flow without the effects of an interstitial
! gas, set SWITCH=0.
      DOUBLE PRECISION, PARAMETER :: SWITCH=1.d0

! ALPHA is a parameter introduced into the theory of Lun_1984 for
! calculating solids viscosity. It also appears when invoking the
! solids frictional model FRICTION, which uses the Lun et al.
! theory. The factor (2+alpha)/3 was eliminated in the complete
! analysis of Lun et al. but was introduced as an adjustable
! parameter. To recover the original theory alpha should be set to
! 1. For details see Johnson and Jackson, 1987.
      DOUBLE PRECISION, PARAMETER :: ALPHA = 1.6d0


! parameter used in the solids-solids drag model invoked in the
! default KT (Lun_1984). For details see Gera et al., 2004
      DOUBLE PRECISION :: SEGREGATION_SLOPE_COEFFICIENT

! SWITCH_IA enforces consistency in the solids viscosity and
! conductivity so that the results using 2 or more identical
! solids phases are the same as an equivalent single solids
! phase. Set to false to use original (published) theory of
! Iddir-Arastoopour.
      LOGICAL, PARAMETER :: SWITCH_IA = .TRUE.


! PHIP = Specularity coefficient associated with particle wall
! collisions
      DOUBLE PRECISION :: PHIP
! PHIP0 specularity coefficient for r->0
      double precision :: phip0
! k4phi k=7/2*mu*(1+e_w)
      double precision :: k4phi
! e_w = particle-wall coefficient of restitution
      DOUBLE PRECISION :: e_w

! Parameters used in the solids frictional model FRICTION:
! - Fr, N_Pc, D_Pc, and EPS_F_min are all used in the equation for
!   Pc, the critical solids pressure:
!     Fr = Constant with dyne/cm2 units of pressure. It will be
!          automatically converted to Pa in calc_mu_s.f
!     N_Pc = exponent in numerator
!     D_Pc = exponent in denominator
!     EPS_f_min = minimum solids fraction above which friction
!                 kicks in
! - N_Pf appears as an exponent in the equation of state for Pf, the
!   frictional pressure:
! - delta is a small deviation in void fraction near packing where
!   Pc and dPc/deps are calculated.
      DOUBLE PRECISION :: EPS_f_min
      DOUBLE PRECISION :: Fr, N_Pc, D_Pc, N_Pf, delta
      PARAMETER(Fr = 0.5d0, N_Pc=2d0, D_Pc=5d0, N_Pf=1.03d0, delta=1d-2)

! Coefficient of restitution
      DOUBLE PRECISION :: C_e

! (1+C_e)/2.
      DOUBLE PRECISION :: eta

! particle-type dependent rest. coef. for use in GHD theory
      DOUBLE PRECISION :: r_p(DIM_M, DIM_M)

! Coeficient of friction
      DOUBLE PRECISION :: C_f

! Angle of internal friction (degrees)
      DOUBLE PRECISION :: Phi

! Angle of wall-particle friction (degrees)
      DOUBLE PRECISION :: Phi_w

! (k=) Sin(PHI) in frictional-flow stress formulation
      DOUBLE PRECISION :: Sin_Phi

! Sin^2(PHI) in plastic-flow stress formulation
      DOUBLE PRECISION :: Sin2_Phi

! (3-2k^2)/6k^2 in Plastic-flow stress formulation
      DOUBLE PRECISION :: F_Phi

! tan(PHI_w)
      DOUBLE PRECISION :: tan_Phi_w

! Excluded volume (Boyle-Massoudi stress tensor)
      DOUBLE PRECISION :: V_ex

! Coefficients for calibrating Syamlal-O'Brien drag correlation with
! Umf data
      DOUBLE PRECISION :: drag_c1, drag_d1

! success-factor for aggregation and breakage
      DOUBLE PRECISION :: AGGREGATION_EFF
      DOUBLE PRECISION :: BREAKAGE_EFF

! UNIT conversion factor for pressure (Barye to Pa if SI)
      DOUBLE PRECISION :: to_SI

! Gravitational acceleration
      DOUBLE PRECISION :: GRAVITY, GRAVITY_X, GRAVITY_Y, GRAVITY_Z

! Universal gas constant
      DOUBLE PRECISION :: GAS_CONST

! Universal gas constant in cal/mol.K
      DOUBLE PRECISION, PARAMETER :: GAS_CONST_cal = 1.987207D0

! Pi, the ubiquitous irrational number
      DOUBLE PRECISION, PARAMETER :: Pi = 4.D0*ATAN(1.D0)

! Square root of Pi
      DOUBLE PRECISION, PARAMETER :: SQRT_Pi = 2.D0*SQRT(ATAN(1.D0))

! Maximum pressure correction allowed in one iteration
      DOUBLE PRECISION :: MAX_DELP

! User defined constants
      DOUBLE PRECISION :: C (DIMENSION_C)

! Names of user defined constants (for output file only)
      CHARACTER(LEN=20) :: C_NAME (DIMENSION_C)

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_CONSTANTS                                           C
!  Purpose: Set various constants                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_CONSTANTS

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero, one, undefined
      USE run, only: LAM_HYS, UNITS
      USE error_manager, only: err_msg, init_err_msg, finl_err_msg
      USE error_manager, only: flush_err_msg

      IMPLICIT NONE
!---------------------------------------------------------------------//

! Note that the cell flags are not set when this routine is called.

! Enter the value of all constants in various units (CGS or SI)
      IF (UNITS == 'SI') THEN
         IF (GRAVITY == UNDEFINED) GRAVITY = 9.80665D0 ! m/s2
         GAS_CONST = 8314.56D0                !Pa.m3/kmol.K, or kg m2/s2 kmol K (Perry and Green, 1984)
         to_SI = 0.1D0                        !to convert dyne/cm2 to Pa. e.g. calc_mu_g.f
         IF (LAM_HYS == UNDEFINED) LAM_HYS = 0.000001d0    ! m
      ELSEIF (UNITS == 'CGS') THEN
         IF (GRAVITY == UNDEFINED) GRAVITY = 980.665D0 !cm/s2
         GAS_CONST = 8314.56D4                !g.cm2/s2.mol.K
         to_SI = ONE                          !does not do anything in CGS. e.g.: calc_mu_g.f
         IF (LAM_HYS == UNDEFINED) LAM_HYS = 0.0001d0    ! cm
      ELSE
! Initialize the error manager.
         CALL INIT_ERR_MSG("SET_CONSTANTS")
         WRITE(ERR_MSG,1005) UNITS
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         CALL FINL_ERR_MSG
 1005 FORMAT('Error 1005: Unknown UNITS = ',A,/ &
         'Please correct the mfix.dat file.')
      ENDIF

! If all components of the gravitational acceleration vector are
! undefined (zero), then use the default value for the negative
! y-direction. This ensures backward compatibility with the old
! (legacy) GRAVITY keyword. At this point GRAVITY is defined,
! either from mfix.dat or by default above
      IF(GRAVITY_X==ZERO.AND.GRAVITY_Y==ZERO.AND.GRAVITY_Z==ZERO) THEN
         GRAVITY_X = ZERO
         GRAVITY_Y = - GRAVITY
         GRAVITY_Z = ZERO
      ENDIF

      RETURN
      END SUBROUTINE SET_CONSTANTS

      END MODULE constant
