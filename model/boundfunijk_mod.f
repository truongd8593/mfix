  MODULE BOUNDFUNIJK

  CONTAINS

      INTEGER FUNCTION BOUND_FUNIJK(pLI, pLJ, pLK)

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
      INTEGER          pLI, pLJ, pLK

!                      Dummy indices for MINS AND MAXS
      INTEGER          LIMIN, LJMIN, LKMIN, LIMAX, LJMAX, LKMAX
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IJK

!-----------------------------------------------
      INCLUDE 'function.inc'

      BOUND_FUNIJK  = FUNIJK ( MIN( IEND3, MAX (ISTART3, pLI) ),&
                               MIN( JEND3, MAX (JSTART3, pLJ) ),&
                               MIN( KEND3, MAX (KSTART3, pLK) ) )

      END FUNCTION BOUND_FUNIJK

  END MODULE BOUNDFUNIJK
