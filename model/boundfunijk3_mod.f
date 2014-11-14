  MODULE BOUNDFUNIJK3

  CONTAINS

      INTEGER FUNCTION BOUND_FUNIJK3(pLI, pLJ, pLK)

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
      USE function3
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

      BOUND_FUNIJK3  = FUNIJK3 ( MIN( IEND4, MAX (ISTART4, pLI) ),&
                               MIN( JEND4, MAX (JSTART4, pLJ) ),&
                               MIN( KEND4, MAX (KSTART4, pLK) ) )

      END FUNCTION BOUND_FUNIJK3

  END MODULE BOUNDFUNIJK3
