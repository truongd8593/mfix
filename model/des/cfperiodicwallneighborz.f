!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CF_PERIODIC_WALL_NEIGHBOR_Z(L, WW, EW)                C
!  Purpose: DES - FInd neighbors at periodic z-walls                   C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  REVISED Aug 24 2005                                                  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CF_PERIODIC_WALL_NEIGHBOR_Z(L, WW, EW)
      USE param1
      USE discretelement

      INTEGER WW(1000), EW(1000), I, II, K, J, L
      DOUBLE PRECISION R_LM, DIST

      IF(DES_POS_NEW(L,3).GE.(NZ2 - RADIUS_EQ)) THEN

         IF(WW(1).GT.0) THEN
            DO J = 2, WW(1)+1
               I = WW(J)
               DES_POS_NEW(I,3) = DES_POS_NEW(I,3) + (NZ2 - SZ1)
               DIST = ZERO
               DO II = 1, DIMN
                  DIST = DIST +&
                  (DES_POS_NEW(I,II)-DES_POS_NEW(L,II))**2
               END DO
               DIST = SQRT(DIST)
               R_LM = DES_RADIUS(L) + DES_RADIUS(I)
               IF (DIST.LE.R_LM) THEN
                  NEIGHBOURS(L,1) = NEIGHBOURS(L,1)+1
                  K = NEIGHBOURS(L,1) + 1
                  IF (K.LE.MN) THEN
                     NEIGHBOURS(L,K) = I
                  ELSE
                     PRINT *,'CFPERIODICWALLNEIGHBOURS K.GE.MN', L, K, (NEIGHBOURS(L,K), K=1, MAXNEIGHBORS)
                     STOP
                  END IF
               END IF
               DES_POS_NEW(I,3) = DES_POS_NEW(I,3) - (NZ2 - SZ1)
            END DO
         END IF

      ELSE IF(DES_POS_NEW(L,3).LE.(SZ1 + RADIUS_EQ)) THEN

         IF(EW(1).GT.0) THEN
            DO J = 2, EW(1)+1
               I = EW(J)
               DES_POS_NEW(I,3) = DES_POS_NEW(I,3) - (NZ2 - SZ1)
               DIST = ZERO
               DO II = 1, DIMN
                  DIST = DIST +&
                  (DES_POS_NEW(I,II)-DES_POS_NEW(L,II))**2
               END DO
               DIST = SQRT(DIST)
               R_LM = DES_RADIUS(L) + DES_RADIUS(I)
               IF (DIST.LE.R_LM) THEN
                  NEIGHBOURS(L,1) = NEIGHBOURS(L,1)+1
                  K = NEIGHBOURS(L,1) + 1
                  IF (K.LE.MN) THEN
                     NEIGHBOURS(L,K) = I
                  ELSE
                     PRINT *,'CFPERIODICWALLNEIGHBOURS K.GE.MN', L, K, (NEIGHBOURS(L,K), K=1, MAXNEIGHBORS)
                     STOP
                  END IF
               END IF
               DES_POS_NEW(I,3) = DES_POS_NEW(I,3) + (NZ2 - SZ1)
            END DO
         END IF

      END IF

      RETURN
      END SUBROUTINE CF_PERIODIC_WALL_NEIGHBOR_Z


