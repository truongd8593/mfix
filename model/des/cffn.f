!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_VEL_STAR(IER)                                    C
!  Purpose: DES - Calculate the normal force on a particle due to      C
!           interparticle collision                                    C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFN(L, Vno, OVERLP_N, NORM)
      
      USE discretelement
      IMPLICIT NONE
      
      INTEGER K, L
      DOUBLE PRECISION OVERLP_N, NORM(NDIM), Vno
!     
!---------------------------------------------------------------------
!     

      DO K = 1, DIMN
         FNS1(K) = 0 -  KN*OVERLP_N*NORM(K)
         FN(K,L) = FNS1(K) -  ETA_DES_N*Vno*NORM(K)
      END DO

!     PRINT *, FN(1,L), FN(2,L)
      
      RETURN
      END SUBROUTINE CFFN


