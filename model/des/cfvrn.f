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

      USE param1      
      USE discretelement
      IMPLICIT NONE

      INTEGER K
      DOUBLE PRECISION Vno, VRl(DIMN), NORM(DIMN)

!---------------------------------------------------------------------------

      Vno = ZERO	

      DO K = 1, DIMN
         Vno = Vno + VRl(K)*NORM(K)
      END DO 

      RETURN
      END SUBROUTINE CFVRN 


