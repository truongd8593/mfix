!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFCTOWALL(L, NORM)                                       C
!  Purpose: DES - Calculate the total force and torque on a particle   C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFCTOWALL(L, NORM)
     
      USE param1 
      USE discretelement
      IMPLICIT NONE
      
      INTEGER L, K
      DOUBLE PRECISION NORM(DIMN), CROSSP(DIMN), FT1(DIMN)
!---------------------------------------------------------------------
!     

        FC(L,:) = FC(L,:) + FN(L,:) + FT(L,:) 

        FT1(:) = FT(L,:)
	CALL DES_CROSSPRDCT(CROSSP, NORM, FT1)         

        TOW(L,:) = TOW(L,:) + DES_RADIUS(L)*CROSSP(:)

      RETURN
      END SUBROUTINE CFFCTOWALL


