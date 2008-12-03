!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NSQUARE                                                C
!>  Purpose: DES - N-Square neighbor search  
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE NSQUARE

      USE param1
      USE discretelement
      IMPLICIT NONE

      INTEGER L, I, K, KK, II, TC1, TC2, TCR, TCM, J
      DOUBLE PRECISION CT
      DOUBLE PRECISION DIST, R_LM

      DO L = 1, PARTICLES 
         DO I = 1, PARTICLES
            DIST = ZERO           
            IF (I.GT.L) THEN
               DO II = 1, DIMN
                  DIST = DIST + (DES_POS_NEW(I,II)-DES_POS_NEW(L,II))**2 
               END DO
               DIST = SQRT(DIST)
               R_LM = DES_RADIUS(L) + DES_RADIUS(I)
               IF (DIST.LE.R_LM) THEN
                  NEIGHBOURS(L,1) = NEIGHBOURS(L,1)+1
                  NEIGHBOURS(I,1) = NEIGHBOURS(I,1)+1
                  K = NEIGHBOURS(L,1) 
                  KK = NEIGHBOURS(I,1)
                  IF (K.LE.MN) THEN
                     NEIGHBOURS(L,K+1) = I
                     PN_DIST(L,K+1) = DIST
                     PN_RLM(L,K+1) = R_LM
                  ELSE 
                     PRINT *,'NSQUARE - NEIGHBORS GT MN'
                     PRINT *, L,':',(NEIGHBOURS(L,II), II=1,MAXNEIGHBORS) 
                     STOP
                  END IF
                  IF (KK.LE.MN) THEN
                     NEIGHBOURS(I,KK+1) = L
                     PN_DIST(I,KK+1) = DIST
                     PN_RLM(I,KK+1) = R_LM
                  ELSE 
                     PRINT *,'NSQUARE - NEIGHBORS GT MN'
                     PRINT *, I,':',(NEIGHBOURS(I,II), II=1,MAXNEIGHBORS) 
                     STOP
                  END IF
               END IF
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE NSQUARE


