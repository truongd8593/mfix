!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  BODY_FORCE.INC                                        C
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
!   Body forces for a horizontal rotary drum
!

!  The body force must be specified at the momentum cell center.  Recall
!  that the momentum cells are staggered.  Therefore use Z for the x-momentum
!  cell and Z_T for the z-momentum cell.

!
!  Body force on gas at i+1/2, j, k
      BFX_g(IJK) = -GRAVITY * sin(Z(K_of(IJK)))  
!
!  Body force on gas at i, j+1/2, k
      BFY_g(IJK) = ZERO
!
!  Body force on gas at i, j, k+1/2
      BFZ_g(IJK) = -GRAVITY * cos(Z_T(K_of(IJK)))
!
!  Body force on solids m at i+1/2, j, k
      BFX_s(IJK, M) = -GRAVITY * sin(Z(K_of(IJK)))
!
!  Body force on solids m at i, j+1/2, k
      BFY_s(IJK, M) = ZERO
!
!  Body force on solids m at i, j, k+1/2
      BFZ_s(IJK, M) = -GRAVITY * cos(Z_T(K_of(IJK)))
