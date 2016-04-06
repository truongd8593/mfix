! this routine is not being used - should be removed

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_L_scale                                             C
!  Purpose: Initialize length scale for turbulence model               C
!                                                                      C
!  Author: W. Sams                                    Date: 04-MAY-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_L_SCALE

! Modules
!--------------------------------------------------------------------//
      USE param
      USE param1
      USE parallel
      USE constant
      USE turb, only: l_scale0, l_scale
      USE geometry
      USE indices
      USE compar
      IMPLICIT NONE

! Local Variables
!--------------------------------------------------------------------//
      INTEGER :: IJK

!--------------------------------------------------------------------//
      IJK = 1

!     IF (IJKMAX2 > 0) THEN
         L_SCALE(IJKSTART3:IJKEND3) = L_SCALE0
!        IJK = IJKMAX3 + 1
!     ENDIF
      RETURN
      END SUBROUTINE SET_L_SCALE
