!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFVRN(Vno, VRl, NORM)                                  C
!  Purpose: DES - Calculate the normal component of relative velocity  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFVRN(Vno, VRl, NORM)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K
      DOUBLE PRECISION Vno, VRl(NDIM), NORM(NDIM)

!---------------------------------------------------------------------------

      Vno = 0.0	

      DO K = 1, DIMN
         Vno = Vno + VRl(K)*NORM(K)
      END DO 

!     PRINT *,'Vrn', Vno

      RETURN
      END SUBROUTINE CFVRN 


