!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NSQUARE(PARTS)                                         C
!  Purpose: DES - N-Square neighbor search                             C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE NSQUARE(PARTS)

      USE discretelement
      IMPLICIT NONE

      INTEGER L, I, K, II, TC1, TC2, TCR, TCM, PARTS, J
      DOUBLE PRECISION CT
      DOUBLE PRECISION DIST, R_LM

      IF (DO_NSQUARE) THEN
         CALL SYSTEM_CLOCK(TC1,TCR,TCM)
         IF(CALLED.EQ.0) THEN
            PRINT *,'N SQUARE'
         END IF
         DO L = 1, PARTICLES 
            DO I = 1, PARTICLES
               R_LM = 0           
               IF (I.NE.L) THEN
                  DO II = 1, DIMN
                     R_LM = R_LM + (DES_POS_NEW(II,I)-DES_POS_NEW(II,L))**2 
                  END DO
                  R_LM = SQRT(R_LM)
                  DIST = DES_RADIUS(L) + DES_RADIUS(I)
                  IF (R_LM.LE.DIST) THEN
                     NEIGHBOURS(1,L) = NEIGHBOURS(1,L)+1
                     K = NEIGHBOURS(1,L) + 1
                     IF (K.LE.MN) THEN
                        NEIGHBOURS(K,L) = I
                     ELSE 
                     PRINT *,'NEIGHBOURS: K GT MN', CALLED
                        STOP
                     END IF
                  END IF
               END IF
            END DO
         END DO

         CALL SYSTEM_CLOCK(TC2,TCR,TCM)
         CT = TC2-TC1
         IF(CT.LE.0) THEN
            CT = TC2 + TCM - TC1
         END IF
         CT = CT/TCR
         N2CT = CT
!     PRINT *,'N2:- CPU TIME TAKEN:',N2CT
      ELSE
         N2CT = 0.000
      END IF

      RETURN
      END SUBROUTINE NSQUARE


