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
      DOUBLE PRECISION OVERLP_N, NORM(DIMN), Vno
!     
!---------------------------------------------------------------------
!     

         FNS1(:) = -KN*OVERLP_N*NORM(:)
         FN(L,:) = FNS1(:) -  ETA_DES_N*Vno*NORM(:)

      RETURN
      END SUBROUTINE CFFN


