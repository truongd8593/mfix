!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFTANGENT(SVL, TANGNT, NORM)                           C
!  Purpose: DES - Calculate tangent vector between partcile pair       C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Comments: Implements Eqn 12 from the following paper                C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFTANGENT(SVL, TANGNT, NORM)
      
      USE discretelement
      USE param1
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER K, L, II
      DOUBLE PRECISION SVL(DIMN), TANGNT(DIMN), NORM(DIMN)
      DOUBLE PRECISION TANMOD 
!     
!---------------------------------------------------------------------


         TANMOD = ZERO

! Tangent computed from relative velocity and used in cfslide

         TANMOD = SQRT(DES_DOTPRDCT(SVL,SVL))     
         IF(TANMOD.NE.0) THEN
            TANGNT(:) = SVL(:)/TANMOD
         ELSE
            TANGNT(1) = NORM(2)
            TANGNT(2) = -NORM(1)
         END IF

      RETURN
      END SUBROUTINE CFTANGENT

