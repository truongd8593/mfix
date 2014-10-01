!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  BODYFORCE_MOD                                         C
!  Purpose: Include file for all body force statement functions        C
!                                                                      C
!  Author: M. Syamlal                                 Date:  6-MAR-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
MODULE bodyforce

  USE constant
  USE sendrecv

CONTAINS

  !  Body force on gas at i+1/2, j, k
  DOUBLE PRECISION FUNCTION BFX_g(IJK)
    IMPLICIT NONE
    INTEGER ijk
    BFX_g = -GRAVITY * sin(Z(K_of(IJK)))
  END FUNCTION BFX_g

  !  Body force on gas at i, j+1/2, k
  DOUBLE PRECISION FUNCTION BFY_g(IJK)
    IMPLICIT NONE
    INTEGER ijk
    BFY_g = ZERO
  END FUNCTION BFY_g

  !  Body force on gas at i, j, k+1/2
  DOUBLE PRECISION FUNCTION BFZ_g(IJK)
    IMPLICIT NONE
    INTEGER ijk
    BFZ_g = -GRAVITY * cos(Z_T(K_of(IJK)))
  END FUNCTION BFZ_g

  !  Body force on solids m at i+1/2, j, k
  DOUBLE PRECISION FUNCTION BFX_s(IJK,M)
    IMPLICIT NONE
    INTEGER ijk,m
    BFX_s = -GRAVITY * sin(Z(K_of(IJK)))
  END FUNCTION BFX_s

  !  Body force on solids m at i, j+1/2, k
  DOUBLE PRECISION FUNCTION BFY_s(IJK,M)
    IMPLICIT NONE
    INTEGER ijk,m
    BFY_s = ZERO
  END FUNCTION BFY_s

  !  Body force on solids m at i, j, k+1/2
  DOUBLE PRECISION FUNCTION BFZ_s(IJK,M)
    IMPLICIT NONE
    INTEGER ijk,m
    BFZ_s = -GRAVITY * cos(Z_T(K_of(IJK)))
  END FUNCTION BFZ_s

END MODULE bodyforce
