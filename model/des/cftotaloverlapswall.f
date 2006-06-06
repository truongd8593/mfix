!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFTOTALOVERLAPSWALL(L, II, Vtan, OVERLP_N, OVERLP_T)   C
!  Purpose:  DES - Calculate the total overlap between particles       C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFTOTALOVERLAPSWALL(L, II, Vtan, OVERLP_N, OVERLP_T)

      USE param1      
      USE discretelement
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER K, L, II
      DOUBLE PRECISION OVERLP_N, OVERLP_T
      DOUBLE PRECISION Vtan, R_LM, D(DIMN), DIST, TEMPX, TEMPY, TEMPZ

!-----------------------------------------------------------------------

      R_LM = DES_RADIUS(L) + DES_RADIUS(II)
      D(:) = DES_POS_NEW(L,:) - DES_POS_NEW(II,:)
      DIST = SQRT(DES_DOTPRDCT(D,D))

      OVERLP_N = R_LM - DIST
      OVERLP_T = Vtan*DTSOLID

      RETURN
      END SUBROUTINE CFTOTALOVERLAPSWALL


