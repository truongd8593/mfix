!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFTOTALOVERLAPS(L, II, Vtan, OVERLP_N, OVERLP_T)       C
!  Purpose:  DES - Calculate the total overlap between particles       C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFTOTALOVERLAPS(L, II, Vtan, OVERLP_N, OVERLP_T)

      USE param1      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L, II
      DOUBLE PRECISION OVERLP_N, OVERLP_T
      DOUBLE PRECISION Vtan, R_LM, DIST, TEMPX, TEMPY, TEMPZ

!-----------------------------------------------------------------------

      DIST = ZERO 
      R_LM = DES_RADIUS(L) + DES_RADIUS(II)

      DO K = 1, DIMN
         DIST = DIST + (DES_POS_NEW(L,K)-DES_POS_NEW(II,K))**2
      END DO
      DIST = SQRT(DIST)

      IF(DES_PERIODIC_WALLS) THEN
        TEMPX = DES_POS_NEW(II,1)
        TEMPY = DES_POS_NEW(II,2)
        IF(DIMN.EQ.3) TEMPZ = DES_POS_NEW(II,3)
        IF(DES_PERIODIC_WALLS_X) THEN        
          IF(DIST.GT.R_LM) THEN
             IF(DES_POS_NEW(II,1).GE.(EX2 - 2*RADIUS_EQ)) THEN
                DES_POS_NEW(II,1) = DES_POS_NEW(II,1) - (EX2 - WX1)
             ELSE IF(DES_POS_NEW(II,1).LE.(WX1 + 2*RADIUS_EQ)) THEN
                DES_POS_NEW(II,1) = DES_POS_NEW(II,1) + (EX2 - WX1)
             END IF
          END IF
        END IF
        IF(DES_PERIODIC_WALLS_Y) THEN
           IF(DIST.GT.R_LM) THEN
              IF(DES_POS_NEW(II,2).GE.(TY2 - 2*RADIUS_EQ)) THEN
                 DES_POS_NEW(II,2) = DES_POS_NEW(II,2) - (TY2 - BY1)
              ELSE IF(DES_POS_NEW(II,2).LE.(BY1 + 2*RADIUS_EQ)) THEN
                 DES_POS_NEW(II,2) = DES_POS_NEW(II,2) + (TY2 - BY1)
              END IF
           END IF
        END IF
        IF(DES_PERIODIC_WALLS_Z) THEN
           IF(DIMN.EQ.3) THEN
           IF(DIST.GT.R_LM) THEN
              IF(DES_POS_NEW(II,3).GE.(NZ2 - 2*RADIUS_EQ)) THEN
                 DES_POS_NEW(II,3) = DES_POS_NEW(II,3) - (NZ2 - SZ1)
              ELSE IF(DES_POS_NEW(II,3).LE.(SZ1 + 2*RADIUS_EQ)) THEN
                 DES_POS_NEW(II,3) = DES_POS_NEW(II,3) + (NZ2 - SZ1)
              END IF
           END IF
           END IF
        END IF                   
      DIST = 0.0
      DO K = 1, DIMN
         DIST = DIST + (DES_POS_NEW(L,K)-DES_POS_NEW(II,K))**2
      END DO
      DIST = SQRT(DIST)
      DES_POS_NEW(II,1) = TEMPX
      DES_POS_NEW(II,2) = TEMPY
      IF (DIMN.EQ.3) DES_POS_NEW(II,3) = TEMPZ
      END IF
 
      OVERLP_N = R_LM - DIST
      OVERLP_T = Vtan*DTSOLID

      RETURN
      END SUBROUTINE CFTOTALOVERLAPS


