!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CF_PERIODIC_WALL_NEIGHBOUR_Y(L, WW, EW)                C
!  Purpose: DES - Find neighbours at periodic y-walls                  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CF_PERIODIC_WALL_NEIGHBOUR_Y(L, WW, EW)

      USE discretelement

      INTEGER WW(1000), EW(1000), I, II, K, J, L
      DOUBLE PRECISION R_LM, DIST

      IF(DES_POS_NEW(2,L).GE.(TY2 - 2*RADIUS_EQ)) THEN

         IF(WW(1).GT.0) THEN
            DO J = 2, WW(1)+1
               I = WW(J)
               DES_POS_NEW(2,I) = DES_POS_NEW(2,I) + (TY2 - BY1)
               DIST = 0
               DO II = 1, DIMN
                  DIST = DIST +&
                  (DES_POS_NEW(II,I)-DES_POS_NEW(II,L))**2
               END DO
               DIST = SQRT(DIST)
               R_LM = DES_RADIUS(L) + DES_RADIUS(I)
               IF (DIST.LE.R_LM) THEN
                  NEIGHBOURS(1,L) = NEIGHBOURS(1,L)+1
                  K = NEIGHBOURS(1,L) + 1
                  IF (K.LE.MN) THEN
                     NEIGHBOURS(K,L) = I
                  ELSE
                     PRINT *,'CFPERIODICWALLNEIGHBOURS K.GE.MN', L, K, (NEIGHBOURS(K,L), K=1, MAXNEIGHBORS)
                     STOP
                  END IF
               END IF
               DES_POS_NEW(2,I) = DES_POS_NEW(2,I) - (TY2 - BY1)
            END DO
         END IF

      ELSE IF(DES_POS_NEW(2,L).LE.(WX1 + 2*RADIUS_EQ)) THEN

         IF(EW(1).GT.0) THEN
            DO J = 2, EW(1)+1
               I = EW(J)
               DES_POS_NEW(2,I) = DES_POS_NEW(2,I) - (TY2 - BY1)
               DIST = 0
               DO II = 1, DIMN
                  DIST = DIST +&
                  (DES_POS_NEW(II,I)-DES_POS_NEW(II,L))**2
               END DO
               DIST = SQRT(DIST)
               R_LM = DES_RADIUS(L) + DES_RADIUS(I)
               IF (DIST.LE.R_LM) THEN
                  NEIGHBOURS(1,L) = NEIGHBOURS(1,L)+1
                  K = NEIGHBOURS(1,L) + 1
                  IF (K.LE.MN) THEN
                     NEIGHBOURS(K,L) = I
                  ELSE
                     PRINT *,'CFPERIODICWALLNEIGHBOURS K.GE.MN', L, K, (NEIGHBOURS(K,L), K=1, MAXNEIGHBORS)
                     STOP
                  END IF
               END IF
               DES_POS_NEW(2,I) = DES_POS_NEW(2,I) + (TY2 - BY1)
            END DO
         END IF

      END IF

      RETURN
      END SUBROUTINE CF_PERIODIC_WALL_NEIGHBOUR_Y


