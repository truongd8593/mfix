!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CF_PERIODIC_WALL_NEIGHBOR_X(L)                         C
!  Purpose: DES - Find neighbors at periodic x-walls                   C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  REVISED Aug 24 2005                                                 C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CF_PERIODIC_WALL_NEIGHBOR_X(L)
      USE param1
      USE discretelement

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER I, II, K, J, L
      DOUBLE PRECISION R_LM, DIST(DIMN), DISP, DIA_EQ

      DIA_EQ = 2D0*DIA_EQ

      IF(DES_POS_NEW(L,1).GE.(EX2 - DIA_EQ)) THEN

         IF(WWALL(1).GT.0) THEN
            DO J = 2, WWALL(1)+1
               I = WWALL(J)
               IF((DES_POS_NEW(I,2).GE.(DES_POS_NEW(L,2)-DIA_EQ)).AND.&
		  (DES_POS_NEW(I,2).LE.(DES_POS_NEW(L,2)+DIA_EQ))) THEN
               DES_POS_NEW(I,1) = DES_POS_NEW(I,1) + (EX2 - WX1)
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
               DES_POS_NEW(I,1) = DES_POS_NEW(I,1) - (EX2 - WX1)
               END IF
            END DO
         END IF

      ELSE IF(DES_POS_NEW(L,1).LE.(WX1 + DIA_EQ)) THEN

         IF(EWALL(1).GT.0) THEN
            DO J = 2, EWALL(1)+1
               I = EWALL(J)
               IF((DES_POS_NEW(I,2).GE.(DES_POS_NEW(L,2)-DIA_EQ)).AND.&
		  (DES_POS_NEW(I,2).LE.(DES_POS_NEW(L,2)+DIA_EQ))) THEN
               DES_POS_NEW(I,1) = DES_POS_NEW(I,1) - (EX2 - WX1)
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
               DES_POS_NEW(I,1) = DES_POS_NEW(I,1) + (EX2 - WX1)
               END IF
           END DO
       END IF

      END IF

      RETURN
      END SUBROUTINE CF_PERIODIC_WALL_NEIGHBOR_X


