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
      IMPLICIT NONE

      INTEGER L, K, SIGNCONV(NDIM)
      DOUBLE PRECISION FTMD, FNMD, TANGNT(NDIM)
!     
!---------------------------------------------------------------------


      DO K=1, DIMN
         IF(FT(K,L).GT.0) THEN
            SIGNCONV(K) = 1
         ELSE
            SIGNCONV(K) = -1
         END IF
      END DO

      FTMD = 0.0
      FNMD = 0.0

      DO K = 1, DIMN
         FTMD = FTMD + (FT(K,L)**2)
         FNMD = FNMD + (FN(K,L)**2)
      END DO

      FTMD = SQRT(FTMD)
      FNMD = SQRT(FNMD)

      IF (FTMD.GT.(MEW_W*FNMD)) THEN
         DO K = 1, DIMN 
            FT(K,L) = 0 - MEW_W*FNMD*TANGNT(K)
         END DO

         DO K = 1, DIMN
            IF(SIGNCONV(K).GT.0) THEN
               IF(FT(K,L).LT.0) THEN
                  FT(K,L) = 0 - FT(K,L)
               END IF
            ELSE IF(SIGNCONV(K).LT.0) THEN
               IF(FT(K,L).GT.0) THEN
                  FT(K,L) = 0 - FT(K,L)
               END IF
            END IF
         END DO

      END IF

!     PRINT *,'SLIDE WALL'

      RETURN
      END SUBROUTINE CFSLIDEWALL


