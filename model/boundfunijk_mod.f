  MODULE BOUNDFUNIJK

  CONTAINS

      INTEGER FUNCTION BOUND_FUNIJK(LI, LJ, LK)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE geometry
      USE compar
      USE fldvar
      USE indices
      IMPLICIT NONE
!
!                      Dummy indices for I, J, K
      INTEGER          LI, LJ, LK

!                      Dummy indices for MINS AND MAXS
      INTEGER          LIMIN, LJMIN, LKMIN, LIMAX, LJMAX, LKMAX
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IJK

!-----------------------------------------------
      INCLUDE 'function.inc'

      BOUND_FUNIJK  = FUNIJK ( MIN( IEND3, MAX (ISTART3, LI) ),&
                               MIN( JEND3, MAX (JSTART3, LJ) ),&
                               MIN( KEND3, MAX (KSTART3, LK) ) )

      END FUNCTION BOUND_FUNIJK

  END MODULE BOUNDFUNIJK
