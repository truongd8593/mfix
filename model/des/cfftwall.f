!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFTWALL(L, VRT, T_OVERLAP, TANGNT)                    C
!  Purpose: DES - Calculate the tangential force on a particle         C
!           due to particle-wall collision                             C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFTWALL(L, VRT, T_OVERLAP, TANGNT)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER L, K
      DOUBLE PRECISION T_OVERLAP, TANGNT(DIMN), VRT
!     
!---------------------------------------------------------------------
!     

         FT(L,:) = -(KT_W*T_OVERLAP + ETA_T_W*VRT)*TANGNT(:)

      RETURN
      END SUBROUTINE CFFTWALL


