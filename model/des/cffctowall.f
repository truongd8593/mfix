!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFCTOWALL(L, NORM)                                    C
!  Purpose: DES - Calculate the total force and torque on a particle   C
!  Reviewer: Rahul Garg                               Date: 02-Aug-07  C
!  Comments: 2-D case torque calculation corrected                     C
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
      
      IF(DIMN.EQ.3) THEN 
         CALL DES_CROSSPRDCT(CROSSP, NORM, FT1)
         TOW(L,:) = TOW(L,:) + DES_RADIUS(L)*CROSSP(:)
      ELSE 
         CROSSP(1) = NORM(1)*FT1(2) - NORM(2)*FT1(1)
         TOW(L,1) = TOW(L,1) + DES_RADIUS(L)*CROSSP(1)
      endif 

      RETURN
      END SUBROUTINE CFFCTOWALL


