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
      USE constant
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

      LIMIN = IMIN2
      LJMIN = JMIN2
      LKMIN = KMIN2

      LIMAX = IMAX2
      LJMAX = JMAX2
      LKMAX = KMAX2

      IF(CYCLIC_X) THEN
        LIMIN = IMIN3
        LIMAX = IMAX3
      ENDIF

      IF(CYCLIC_Y) THEN
        LJMIN = JMIN3
        LJMAX = JMAX3
      ENDIF

      IF(CYCLIC_Z) THEN
        LKMIN = KMIN3
        LKMAX = KMAX3
      ENDIF

      BOUND_FUNIJK  = FUNIJK ( MIN( LIMAX, MAX (LIMIN, LI) ),&
                               MIN( LJMAX, MAX (LJMIN, LJ) ),&
                               MIN( LKMAX, MAX (LKMIN, LK) ) )

      END FUNCTION BOUND_FUNIJK

  END MODULE BOUNDFUNIJK
