!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFNWALL(L, Vno, OVERLP_N, NORM)                       C
!  Purpose: DES - Calclate the normal force on a particle dure to      C
!           particle-wall collision                                    C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFNWALL(L, Vno, OVERLP_N, NORM)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L
      DOUBLE PRECISION OVERLP_N, NORM(DIMN), Vno
!     
!---------------------------------------------------------------------

         FN(L,:) = -(KN_W*OVERLP_N + ETA_N_W*Vno)*NORM(:)

      RETURN
      END SUBROUTINE CFFNWALL


