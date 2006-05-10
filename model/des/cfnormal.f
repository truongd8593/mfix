!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNORMAL(L, II, NORM)                                  C
!  Purpose: DES - Calculate the normal vector between particles        C
!           in collision                                               C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNORMAL(L, II, NORM) 

      USE param1      
      USE discretelement
      IMPLICIT NONE

      INTEGER L, II, K
      DOUBLE PRECISION NORMOD, NORM(DIMN), TEMPX, TEMPY, TEMPZ, TEMPD

!----------------------------------------------------------------------------

      NORMOD = ZERO 

      IF(DES_PERIODIC_WALLS) THEN
      TEMPX = DES_POS_NEW(II,1)
      TEMPY = DES_POS_NEW(II,2)
      IF(DIMN.EQ.3) TEMPZ = DES_POS_NEW(II,3)        
      IF(DES_PERIODIC_WALLS_X) THEN 
        TEMPD = DES_POS_NEW(L,1) - DES_POS_NEW(II,1)
        TEMPD = SQRT(TEMPD**2)
        IF(TEMPD.GT.(DES_RADIUS(L)+DES_RADIUS(II))) THEN
           IF(DES_POS_NEW(II,1).GE.(EX2 - 2*RADIUS_EQ)) THEN
              DES_POS_NEW(II,1) = DES_POS_NEW(II,1) - (EX2-WX1)
           ELSE IF(DES_POS_NEW(II,1).LE.(WX1 + 2*RADIUS_EQ)) THEN
                   DES_POS_NEW(II,1) = DES_POS_NEW(II,1) + (EX2-WX1)
           END IF
        END IF
      ELSE IF(DES_PERIODIC_WALLS_Y) THEN
             TEMPD = DES_POS_NEW(L,2) - DES_POS_NEW(II,2)
             TEMPD = SQRT(TEMPD**2)
             IF(TEMPD.GT.(DES_RADIUS(L)+DES_RADIUS(II))) THEN
                IF(DES_POS_NEW(II,2).GE.(TY2 - 2*RADIUS_EQ)) THEN
                   DES_POS_NEW(II,2) = DES_POS_NEW(II,2) - (TY2-BY1)
                ELSE IF(DES_POS_NEW(II,2).LE.(BY1 + 2*RADIUS_EQ)) THEN
                        DES_POS_NEW(II,2) = DES_POS_NEW(II,2) + (TY2-BY1)
                END IF
             END IF
      ELSE IF(DES_PERIODIC_WALLS_Z) THEN
           IF(DIMN.EQ.3) THEN
             TEMPD = DES_POS_NEW(L,3) - DES_POS_NEW(II,3)
             TEMPD = SQRT(TEMPD**2)
             IF(TEMPD.GT.(DES_RADIUS(L)+DES_RADIUS(II))) THEN
                IF(DES_POS_NEW(II,3).GE.(NZ2 - 2*RADIUS_EQ)) THEN
                   DES_POS_NEW(II,3) = DES_POS_NEW(II,3) - (NZ2-SZ1)
                ELSE IF(DES_POS_NEW(II,3).LE.(SZ1 + 2*RADIUS_EQ)) THEN
                        DES_POS_NEW(II,3) = DES_POS_NEW(II,3) + (NZ2-SZ1)
                END IF
             END IF
           ELSE
             PRINT *,'2D PROBLEM: PERIODIC WALLS IN Z'
             STOP
           END IF
      END IF
      END IF
      
      DO K = 1, DIMN
         NORMOD=NORMOD+(DES_POS_NEW(II,K)-DES_POS_NEW(L,K))**2
      END DO

      NORMOD = SQRT(NORMOD)

      DO K = 1, DIMN
         IF(NORMOD.NE.ZERO) THEN
            NORM(K)=(DES_POS_NEW(II,K)-DES_POS_NEW(L,K))/NORMOD
         ELSE 
            PRINT *,'NORMOD IS ZERO', II,L
            PRINT *, (DES_POS_NEW(II,L),L=1,DIMN)
            PRINT *, (DES_POS_NEW(L,II),II=1,DIMN)
            STOP
         END IF
      END DO

      IF(DES_PERIODIC_WALLS) THEN
           DES_POS_NEW(II,1) = TEMPX
           DES_POS_NEW(II,2) = TEMPY
           IF (DIMN.EQ.3) DES_POS_NEW(II,3) = TEMPZ              
      END IF  

      RETURN
      END SUBROUTINE CFNORMAL 


