!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFSLIDE(L, TANGNT, FT1)                                C
!  Purpose: DES - calculate sliding between particles                  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFSLIDE(L, TANGNT, FT1)

      USE discretelement
      IMPLICIT NONE

      INTEGER L, K
      DOUBLE PRECISION FTMD, FNMD, TANGNT(NDIM), FT1(NDIM), FN1(NDIM)
!     
!---------------------------------------------------------------------


      FTMD = 0.0
      FNMD = 0.0
      DO K = 1, DIMN
         FTMD = FTMD + (FT1(K)**2)
         FNMD = FNMD + (FN(K,L)**2)
      END DO
      FTMD = SQRT(FTMD)
      FNMD = SQRT(FNMD)

      IF (FTMD.GT.(MEW*FNMD)) THEN
         DO K = 1, DIMN 
            FT(K,L) = 0 - MEW*FNMD*TANGNT(K)
         END DO

         DO K = 1, DIMN
            IF(FT1(K).GT.0) THEN
               IF(FT(K,L).LT.0) THEN
                  FT(K,L) = 0 - FT(K,L)
               END IF
            ELSE IF(FT1(K).LT.0) THEN
               IF(FT(K,L).GT.0) THEN
                  FT(K,L) = 0 - FT(K,L)
               END IF
            END IF
         END DO
         
      ELSE
         DO K=1,NDIM 	
            FT(K,L) = FT1(K)
         END DO
      END IF

      
      RETURN
      END SUBROUTINE CFSLIDE


