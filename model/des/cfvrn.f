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

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

      INTEGER K
      DOUBLE PRECISION Vno, VRl(DIMN), NORM(DIMN)

!---------------------------------------------------------------------------

      Vno = DES_DOTPRDCT(VRl,NORM)
      
      RETURN
      END SUBROUTINE CFVRN 


