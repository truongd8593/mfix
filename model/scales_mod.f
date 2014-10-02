!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: scales_mod.f                                           C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

MODULE scales

  Use param
  Use param1

! reference pressure
  DOUBLE PRECISION P_ref

! pressure scale
  DOUBLE PRECISION P_scale

CONTAINS

  DOUBLE PRECISION FUNCTION SCALE(XXX)
    IMPLICIT NONE
    DOUBLE PRECISION XXX
    SCALE   = (XXX - P_ref) / P_scale
  END FUNCTION SCALE

  DOUBLE PRECISION FUNCTION UNSCALE(XXX)
    IMPLICIT NONE
    DOUBLE PRECISION XXX
    UNSCALE = (XXX * P_scale + P_ref)
  END FUNCTION UNSCALE

END MODULE scales
