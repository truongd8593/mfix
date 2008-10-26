!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFINCREMENTALOVERLAPS(Vno, Vtan, OVERLP_N, OVERLP_T)   C
!  Purpose: DES - calculate incremental overlaps between particles     C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Rahul Garg                               Date: 01-AUG-07  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFINCREMENTALOVERLAPS(L, II,Vno, Vtan, OVERLP_N, OVERLP_T, CHECK_CON)

      USE discretelement        
      IMPLICIT NONE
      INTEGER :: L, II 
      DOUBLE PRECISION OVERLP_N, OVERLP_T
      DOUBLE PRECISION Vno, Vtan, R_LM, D(DIMN), DIST
      LOGICAL CHECK_CON
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      
!-----------------------------------------------------------------------
      
      R_LM = DES_RADIUS(L) + DES_RADIUS(II)
      D(:) = DES_POS_NEW(L,:) - DES_POS_NEW(II,:)
      DIST = SQRT(DES_DOTPRDCT(D,D))
      
      IF(R_LM - DIST.gt.SMALL_NUMBER) then 

         OVERLP_N =Vno*DTSOLID
         OVERLP_T = (Vtan) * DTSOLID
         CHECk_CON = .TRUE.
         
      else 
         OVERLP_N = ZERO
         OVERLP_T = ZERO
         
         CHECk_CON = .FALSE.

      endif

      RETURN
      END SUBROUTINE CFINCREMENTALOVERLAPS


