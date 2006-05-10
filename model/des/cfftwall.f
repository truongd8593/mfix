!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFTWALL(L, Vtan, OVERLP_T, TANGNT)                    C
!  Purpose: DES - Calculate the tangential force on a particle         C
!           due to particle-wall collision                             C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFTWALL(L, Vtan, OVERLP_T, TANGNT)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER L, K
      DOUBLE PRECISION OVERLP_T, TANGNT(DIMN), Vtan
!     
!---------------------------------------------------------------------
!     

         FT(L,:) = -(KT_W*OVERLP_T + ETA_T_W*Vtan)*TANGNT(:)

      RETURN
      END SUBROUTINE CFFTWALL


