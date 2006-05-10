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

      USE param1
      USE discretelement
      IMPLICIT NONE

      INTEGER L, K
      DOUBLE PRECISION FTMD, FNMD, TANGNT(DIMN), FT1(DIMN)
!     
!---------------------------------------------------------------------


      FTMD = ZERO 
      FNMD = ZERO
      DO K = 1, DIMN
         FTMD = FTMD + (FT1(K)**2)
         FNMD = FNMD + (FN(L,K)**2)
      END DO
      FTMD = SQRT(FTMD)
      FNMD = SQRT(FNMD)

      IF (FTMD.GT.(MEW*FNMD)) THEN
         FT(L,:) = - MEW*FNMD*TANGNT(:)
         DO K = 1, DIMN
            IF(FT1(K).GT.0) THEN
               IF(FT(L,K).LT.0) THEN
                  FT(L,K) = - FT(L,K)
               END IF
            ELSE IF(FT1(K).LT.0) THEN
               IF(FT(L,K).GT.0) THEN
                  FT(L,K) = - FT(L,K)
               END IF
            END IF
         END DO
      ELSE
         FT(L,:) = FT1(:)
      END IF

      
      RETURN
      END SUBROUTINE CFSLIDE


