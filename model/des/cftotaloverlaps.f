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
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L, II
      DOUBLE PRECISION OVERLP_N, OVERLP_T
      DOUBLE PRECISION Vtan, R_LM, DIST, TEMPX, TEMPY, TEMPZ

!-----------------------------------------------------------------------

      DIST = 0.0
      R_LM = DES_RADIUS(L) + DES_RADIUS(II)

      DO K = 1, DIMN
         DIST = DIST + (DES_POS_NEW(K,L)-DES_POS_NEW(K,II))**2
      END DO
      DIST = SQRT(DIST)

      IF(DES_PERIODIC_WALLS) THEN
        TEMPX = DES_POS_NEW(1,II)
        TEMPY = DES_POS_NEW(2,II)
        TEMPZ = DES_POS_NEW(3,II)
        IF(DES_PERIODIC_WALLS_X) THEN        
          IF(DIST.GT.R_LM) THEN
             IF(DES_POS_NEW(1,II).GE.(EX2 - 2*RADIUS_EQ)) THEN
                DES_POS_NEW(1,II) = DES_POS_NEW(1,II) - (EX2 - WX1)
             ELSE IF(DES_POS_NEW(1,II).LE.(WX1 + 2*RADIUS_EQ)) THEN
                DES_POS_NEW(1,II) = DES_POS_NEW(1,II) + (EX2 - WX1)
             END IF
          END IF
        END IF
        IF(DES_PERIODIC_WALLS_Y) THEN
           IF(DIST.GT.R_LM) THEN
              IF(DES_POS_NEW(2,II).GE.(TY2 - 2*RADIUS_EQ)) THEN
                 DES_POS_NEW(2,II) = DES_POS_NEW(2,II) - (TY2 - BY1)
              ELSE IF(DES_POS_NEW(2,II).LE.(BY1 + 2*RADIUS_EQ)) THEN
                 DES_POS_NEW(2,II) = DES_POS_NEW(2,II) + (TY2 - BY1)
              END IF
           END IF
        END IF
        IF(DES_PERIODIC_WALLS_Z) THEN
           IF(DIST.GT.R_LM) THEN
              IF(DES_POS_NEW(3,II).GE.(NZ2 - 2*RADIUS_EQ)) THEN
                 DES_POS_NEW(3,II) = DES_POS_NEW(3,II) - (NZ2 - SZ1)
              ELSE IF(DES_POS_NEW(3,II).LE.(SZ1 + 2*RADIUS_EQ)) THEN
                 DES_POS_NEW(3,II) = DES_POS_NEW(3,II) + (NZ2 - SZ1)
              END IF
           END IF
        END IF                   
      DIST = 0.0
      DO K = 1, DIMN
         DIST = DIST + (DES_POS_NEW(K,L)-DES_POS_NEW(K,II))**2
      END DO
      DIST = SQRT(DIST)
      DES_POS_NEW(1,II) = TEMPX
      DES_POS_NEW(2,II) = TEMPY
      DES_POS_NEW(3,II) = TEMPZ
      END IF
 
      OVERLP_N = R_LM - DIST
      OVERLP_T = Vtan*DTSOLID
!     PRINT *,'TOTAL OVERLAP', DES_RADIUS(L), DES_RADIUS(II), OVERLP_N, OVERLP_T

      RETURN
      END SUBROUTINE CFTOTALOVERLAPS


