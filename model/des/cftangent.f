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
      IMPLICIT NONE

      INTEGER K
      DOUBLE PRECISION NORM(NDIM), TANGNT(NDIM), VRl(NDIM)
      DOUBLE PRECISION VT(NDIM), VN(NDIM)
      DOUBLE PRECISION TANMOD 
!     
!---------------------------------------------------------------------


      IF(DIMN.EQ.2) THEN
         TANGNT(1) = NORM(2)
         TANGNT(2) = 0 - NORM(1)
      ELSE 
         TANMOD = 0D0
         DO K = 1, DIMN     
            VN(K) = VRl(K)*NORM(K)
            VT(K) =  VRl(K) - VN(K)
            TANMOD = TANMOD + VT(K)**2            
         END DO
         TANMOD = SQRT(TANMOD)
         IF(TANMOD.NE.0) THEN
           DO K = 1, DIMN
             TANGNT(K) = VT(K)/TANMOD
           END DO
         END IF
      END IF

!     PRINT *,'TANGENT',  TANGNT(1), TANGNT(2)

      RETURN
      END SUBROUTINE CFTANGENT


