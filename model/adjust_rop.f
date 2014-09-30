!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ADJUST_ROP                                              C
!  Purpose: Remove small negative values of density.                   C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-AUG-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ADJUST_ROP(ROP, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
! density
      DOUBLE PRECISION, INTENT(INOUT) :: ROP(DIMENSION_3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) ROP(IJK) = DMAX1(ZERO,ROP(IJK))
      ENDDO

      RETURN
      END SUBROUTINE ADJUST_ROP


