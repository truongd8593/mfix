!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFVRT(VRT, VRELTRANS, TANGNT)                          C
!  Purpose: DES - Calculate the tangential component of                C
!           relative velocity                                          C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFVRT(VRT, VRELTRANS, TANGNT)

      USE param1      
      USE discretelement
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

      INTEGER K
      DOUBLE PRECISION VRT, VRELTRANS(DIMN), TANGNT(DIMN), NORM(DIMN)

!-----------------------------------------------------------------------

      VRT = DES_DOTPRDCT(VRELTRANS,TANGNT)
      
      RETURN
      END SUBROUTINE CFVRT


