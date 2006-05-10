!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFSLIDEWALL(L, TANGNT)                                 C
!  Purpose: DES - Calculate slide between particles and walls          C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFSLIDEWALL(L, TANGNT)
      
      USE discretelement
      USE param1
      IMPLICIT NONE

      INTEGER L, K
      DOUBLE PRECISION FTMD, FNMD, TANGNT(DIMN), FT1(DIMN)
!     
!---------------------------------------------------------------------


      FT1(:) = FT(L,:)

      FTMD = ZERO
      FNMD = ZERO
      DO K = 1, DIMN
         FTMD = FTMD + (FT(L,K)**2)
         FNMD = FNMD + (FN(L,K)**2)
      END DO
      FTMD = SQRT(FTMD)
      FNMD = SQRT(FNMD)

      IF (FTMD.GT.(MEW_W*FNMD)) THEN
         FT(L,:) = -MEW_W*FNMD*TANGNT(:)
         DO K = 1, DIMN
            IF(FT1(K).GT.0) THEN
               IF(FT(L,K).LT.ZERO) THEN
                  FT(L,K) = -FT(L,K)
               END IF
            ELSE IF(FT1(K).LT.0) THEN
               IF(FT(L,K).GT.ZERO) THEN
                  FT(L,K) = -FT(L,K)
               END IF
            END IF
         END DO
      END IF

      RETURN
      END SUBROUTINE CFSLIDEWALL


