!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFCTOW(L, NORM)                                       C
!  Purpose: DES - Calculate the total force and torque on a particle   C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFCTOW(L, NORM)
      
      USE discretelement
      IMPLICIT NONE
      
      INTEGER L, K
      DOUBLE PRECISION NORM(NDIM), CHECK
!---------------------------------------------------------------------
!     

      DO K = 1, DIMN
         FC(K,L) = FC(K,L) + FN(K,L) + FT(K,L) 
      END DO

      IF(DIMN.EQ.2) THEN
         NORM(3) = 0.0
         FT(3,L) = 0.0
      END IF


      TOW(1,L) = TOW(1,L) + (DES_RADIUS(L)*(NORM(2)*FT(3,L) - NORM(3)*&
      FT(2,L)))
      TOW(2,L) = TOW(2,L) + (DES_RADIUS(L)*(NORM(3)*FT(1,L) - NORM(1)*&
      FT(3,L)))
      TOW(3,L) = TOW(3,L) + (DES_RADIUS(L)*(NORM(1)*FT(2,L) - NORM(2)*&
      FT(1,L)))

!     PRINT *,'FC', FC(1,L), FC(2,L), FC(3,L), 
!     PRINT *, CALLED, TOW(1,L), TOW(2,L), TOW(3,L)

      RETURN
      END SUBROUTINE CFFCTOW


