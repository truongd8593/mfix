!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFVRN(VRN, VRELTRANS, NORM)                            C
!  Purpose: DES - Calculate the normal component of relative velocity  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFVRN(VRN, VRELTRANS, NORM)

      USE param1      
      USE discretelement
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

      INTEGER K
      DOUBLE PRECISION VRN, VRELTRANS(DIMN), NORM(DIMN)

!---------------------------------------------------------------------------

      VRN = DES_DOTPRDCT(VRELTRANS,NORM)
     
      RETURN
      END SUBROUTINE CFVRN 


