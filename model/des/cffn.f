!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFN(L, VRN, N_OVERLAP, NORM)                           C
!  Purpose: DES - Calculate the normal force on a particle due to      C
!           interparticle collision                                    C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments: Implements Eqn 6 from the following paper                 C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFN(L, VRN, N_OVERLAP, NORM)
      
      USE discretelement
      IMPLICIT NONE
      
      INTEGER K, L
      DOUBLE PRECISION N_OVERLAP, NORM(DIMN), VRN
!     
!---------------------------------------------------------------------
!     

         FNS1(:) = -KN*N_OVERLAP*NORM(:)
         FN(L,:) = FNS1(:) -  ETA_DES_N*VRN*NORM(:)

      RETURN
      END SUBROUTINE CFFN


