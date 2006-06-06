!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CF_PERIODIC_WALL_NEIGHBOR_Y(L)                         C
!  Purpose: DES - Find neighbours at periodic y-walls                  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  REVISED Aug 24 2005                                                 C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CF_PERIODIC_WALL_NEIGHBOR_Y(L)
      USE param1
      USE discretelement

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER I, II, K, J, L
      DOUBLE PRECISION R_LM, DIST(DIMN), DISP, DIA_EQ

      DIA_EQ = 2D0*DIA_EQ

      IF(DES_POS_NEW(L,2).GE.(TY2 - DIA_EQ)) THEN

         IF(BWALL(1).GT.0) THEN
            DO J = 2, BWALL(1)+1
               I = BWALL(J)
               DES_POS_NEW(I,2) = DES_POS_NEW(I,2) + (TY2 - BY1)
               DIST(:) = DES_POS_NEW(I,:)-DES_POS_NEW(L,:)
               DISP = SQRT(DES_DOTPRDCT(DIST,DIST))
               R_LM = DES_RADIUS(L) + DES_RADIUS(I)
               IF (DISP.LE.R_LM) THEN
                  NEIGHBOURS(L,1) = NEIGHBOURS(L,1)+1
                  K = NEIGHBOURS(L,1) + 1
                  IF (K.LE.MN) THEN
                     NEIGHBOURS(L,K) = I
                     PN_DIST(L,K) = DISP
                     PN_RLM(L,K) = R_LM
                  ELSE
                     PRINT *,'CFPERIODICWALLNEIGHBOURS K.GE.MN', L, K, (NEIGHBOURS(L,K), K=1, MAXNEIGHBORS)
                     STOP
                  END IF
               END IF
               DES_POS_NEW(I,2) = DES_POS_NEW(I,2) - (TY2 - BY1)
            END DO
         END IF

      ELSE IF(DES_POS_NEW(L,2).LE.(BY1 + DIA_EQ)) THEN

         IF(TWALL(1).GT.0) THEN
            DO J = 2, TWALL(1)+1
               I = TWALL(J)
               DES_POS_NEW(I,2) = DES_POS_NEW(I,2) - (TY2 - BY1)
               DIST(:) = DES_POS_NEW(I,:)-DES_POS_NEW(L,:)
               DISP = SQRT(DES_DOTPRDCT(DIST,DIST))
               R_LM = DES_RADIUS(L) + DES_RADIUS(I)
               IF (DISP.LE.R_LM) THEN
                  NEIGHBOURS(L,1) = NEIGHBOURS(L,1)+1
                  K = NEIGHBOURS(L,1) + 1
                  IF (K.LE.MN) THEN
                     NEIGHBOURS(L,K) = I
                     PN_DIST(L,K) = DISP
                     PN_RLM(L,K) = R_LM
                  ELSE
                     PRINT *,'CFPERIODICWALLNEIGHBOURS K.GE.MN', L, K, (NEIGHBOURS(L,K), K=1, MAXNEIGHBORS)
                     STOP
                  END IF
               END IF
               DES_POS_NEW(I,2) = DES_POS_NEW(I,2) + (TY2 - BY1)
            END DO
         END IF

      END IF

      RETURN
      END SUBROUTINE CF_PERIODIC_WALL_NEIGHBOR_Y


