!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFNWALL(L, VRN, N_OVERLAP, NORM)                      C
!  Purpose: DES - Calclate the normal force on a particle dure to      C
!           particle-wall collision                                    C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFNWALL(L, VRN, N_OVERLAP, NORM)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L
      DOUBLE PRECISION N_OVERLAP, NORM(DIMN), VRN
!     
!---------------------------------------------------------------------

         FN(L,:) = -(KN_W*N_OVERLAP + ETA_N_W*VRN)*NORM(:)

      RETURN
      END SUBROUTINE CFFNWALL


