!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFT(L, Vtan, OVERLP_T, TANGNT                         C
!  Purpose: DES - Calculate tangential force on a particle             C
!           due to inter-particle collision                            C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments: Implements Eqn 7 from the following paper                 C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFT(L, Vtan, OVERLP_T, TANGNT)
      
      USE discretelement
      IMPLICIT NONE
      
      INTEGER L, K
      DOUBLE PRECISION OVERLP_T, TANGNT(DIMN), Vtan
!---------------------------------------------------------------------


         FTS1(:) = -KT*OVERLP_T*TANGNT(:)
         FT(L,:) = FTS1(:) - ETA_DES_T*Vtan*TANGNT(:)

      RETURN
      END SUBROUTINE CFFT


