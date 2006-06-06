!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNORMALWALL(L, II, NORM)                              C
!  Purpose: DES - Calculate the normal vector between particles        C
!           in collision                                               C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments: Normal for Eqn 6  from the following paper                C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNORMALWALL(L, II, NORM) 

      USE param1      
      USE discretelement
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER L, II, K
      DOUBLE PRECISION NORMOD, NORM(DIMN), DIST(DIMN)

!----------------------------------------------------------------------------

      DIST(:) = DES_POS_NEW(II,:) - DES_POS_NEW(L,:)
      NORMOD = SQRT(DES_DOTPRDCT(DIST,DIST))

      IF(NORMOD.NE.ZERO) THEN
         NORM(:)= DIST(:)/NORMOD
      ELSE 
         PRINT *,'WALL NORMOD IS ZERO', II,L
         STOP
      END IF

      RETURN
      END SUBROUTINE CFNORMALWALL 


