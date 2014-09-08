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
module bodyforce

contains

!  Body force on gas at i+1/2, j, k
double precision function BFX_g(IJK)
      integer ijk
      BFX_g = GRAVITY_X
end function
!
!  Body force on gas at i, j+1/2, k
double precision function BFY_g(IJK)
      integer ijk
      BFY_g = GRAVITY_Y
end function
!
!  Body force on gas at i, j, k+1/2
double precision function BFZ_g(IJK)
      integer ijk
      BFZ_g = GRAVITY_Z
end function
!
!  Body force on solids m at i+1/2, j, k
double precision function BFX_s(IJK,M)
      integer ijk,m
      BFX_s = GRAVITY_X
end function
!
!  Body force on solids m at i, j+1/2, k
double precision function BFY_s(IJK,M)
      integer ijk,m
      BFY_s = GRAVITY_Y
end function
!
!  Body force on solids m at i, j, k+1/2
double precision function BFZ_s(IJK,M)
      integer ijk,m
      BFZ_s = GRAVITY_Z
end function

end module bodyforce
