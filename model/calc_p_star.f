!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_P_star(EP_g, P_star, IER)                         C
!  Purpose: Calculate P_star in cells where solids continuity is solvedC
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-AUG-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
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
      SUBROUTINE CALC_P_STAR(EP_G, P_STAR, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE pgcor
      USE pscor
      USE ur_facs 
      USE residual
      USE compar     !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Solids pressure
      DOUBLE PRECISION P_star(DIMENSION_3)
!
!                      Gas volume fraction
      DOUBLE PRECISION EP_g(DIMENSION_3)

!HPF$ align P_star(:) with TT(:)
!HPF$ align EP_g(:) with TT(:)

!
!                      error index
      INTEGER          IER
!
!                      Indices
      INTEGER          IJK

      DOUBLE PRECISION dPs
!-----------------------------------------------
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
!
!
!

!!$omp parallel do private(ijk)
!HPF$ independent
      DO IJK = 1, IJKMAX2 
         IF (FLUID_AT(IJK)) THEN 
!
            IF (EP_G(IJK) < EP_STAR) THEN 
               P_STAR(IJK) = NEG_H(EP_G(IJK)) 
            ELSE 
               P_STAR(IJK) = ZERO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE CALC_P_STAR 
