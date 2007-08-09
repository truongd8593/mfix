!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFTOTALOVERLAPSWALL(L, II, VRT, N_OVERLAP, T_OVERLAP)  C
!  Purpose:  DES - Calculate the total overlap between particles       C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFTOTALOVERLAPSWALL(L, II, VRT, N_OVERLAP, T_OVERLAP, CHECK_CON)

      USE param1      
      USE discretelement
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER K, L, II
      DOUBLE PRECISION N_OVERLAP, T_OVERLAP
      DOUBLE PRECISION D(DIMN), VRT, R_LM, DIST
      LOGICAL CHECK_CON
!-----------------------------------------------------------------------

      R_LM = DES_RADIUS(L) + DES_RADIUS(II)
      D(:) = DES_POS_NEW(L,:) - DES_POS_NEW(II,:)
      DIST = SQRT(DES_DOTPRDCT(D,D))

      IF(R_LM - DIST.gt.SMALL_NUMBER) then 
         
         N_OVERLAP = R_LM - DIST
         T_OVERLAP = (VRT)*DTSOLID
         CHECk_CON = .TRUE.
      else 
         
         N_OVERLAP = zero
         T_OVERLAP = zero
         
         CHECk_CON = .FALSE.
      endif
      
      RETURN
      END SUBROUTINE CFTOTALOVERLAPSWALL


