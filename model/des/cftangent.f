!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFTANGENT(TANGNT, NORM, VRl)                           C
!  Purpose: DES - Calculate tangent vector between partcile pair       C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFTANGENT(TANGNT, NORM, VRl)
      
      USE discretelement
      USE param1
      IMPLICIT NONE

      INTEGER K
      DOUBLE PRECISION NORM(DIMN), TANGNT(DIMN), VRl(DIMN)
      DOUBLE PRECISION VT(DIMN), VN(DIMN)
      DOUBLE PRECISION TANMOD 
!     
!---------------------------------------------------------------------


      IF(DIMN.EQ.2) THEN
         TANGNT(1) = NORM(2)
         TANGNT(2) = -NORM(1)
      ELSE 
         TANMOD = ZERO
         VN(:) = VRl(:)*NORM(:)
         VT(:) =  VRl(:) - VN(:)
         DO K= 1, DIMN
            TANMOD = TANMOD + VT(K)**2     
         END DO       
         TANMOD = SQRT(TANMOD)
         IF(TANMOD.NE.0) THEN
            TANGNT(:) = VT(:)/TANMOD
         END IF
      END IF

      RETURN
      END SUBROUTINE CFTANGENT


