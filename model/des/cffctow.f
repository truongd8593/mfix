!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFCTOW(L, II, NORM)                                   C
!  Purpose: DES - Calculate the total force and torque on a particle   C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments: Implements Eqns 13 & 14 from the following paper          C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFCTOW(L, II,  NORM)
     
      USE param1 
      USE discretelement
      IMPLICIT NONE
      
      INTEGER L, II, K
      DOUBLE PRECISION NORM(DIMN), CROSSP(DIMN), FT1(DIMN)
!---------------------------------------------------------------------
!     

         FC(L,:) = FC(L,:) + FN(L,:) + FT(L,:) 
         FC(II,:) = FC(II,:) - FN(L,:) - FT(L,:)

         FT1(:) = FT(L,:)
	 CALL DES_CROSSPRDCT(CROSSP, NORM, FT1)

         TOW(L,:) = TOW(L,:) + DES_RADIUS(L)*CROSSP(:)
         TOW(II,:) = TOW(L,:) - DES_RADIUS(II)*CROSSP(:)

      RETURN
      END SUBROUTINE CFFCTOW


