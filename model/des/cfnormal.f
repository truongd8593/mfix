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
      
      USE discretelement
      IMPLICIT NONE

      INTEGER L, II, K
      DOUBLE PRECISION NORMOD, NORM(NDIM), TEMPX, TEMPY, TEMPZ, TEMPD

!----------------------------------------------------------------------------

      NORMOD = 0.0

      IF(DES_PERIODIC_WALLS) THEN
      TEMPX = DES_POS_NEW(1,II)
      TEMPY = DES_POS_NEW(2,II)
      TEMPZ = DES_POS_NEW(3,II)        
      IF(DES_PERIODIC_WALLS_X) THEN 
        TEMPD = DES_POS_NEW(1,L) - DES_POS_NEW(1,II)
        TEMPD = SQRT(TEMPD**2)
        IF(TEMPD.GT.(DES_RADIUS(L)+DES_RADIUS(II))) THEN
           IF(DES_POS_NEW(1,II).GE.(EX2 - 2*RADIUS_EQ)) THEN
              DES_POS_NEW(1,II) = DES_POS_NEW(1,II) - (EX2-WX1)
           ELSE IF(DES_POS_NEW(1,II).LE.(WX1 + 2*RADIUS_EQ)) THEN
                   DES_POS_NEW(1,II) = DES_POS_NEW(1,II) + (EX2-WX1)
           END IF
        END IF
      ELSE IF(DES_PERIODIC_WALLS_Y) THEN
             TEMPD = DES_POS_NEW(2,L) - DES_POS_NEW(2,II)
             TEMPD = SQRT(TEMPD**2)
             IF(TEMPD.GT.(DES_RADIUS(L)+DES_RADIUS(II))) THEN
                IF(DES_POS_NEW(2,II).GE.(TY2 - 2*RADIUS_EQ)) THEN
                   DES_POS_NEW(2,II) = DES_POS_NEW(2,II) - (TY2-BY1)
                ELSE IF(DES_POS_NEW(2,II).LE.(BY1 + 2*RADIUS_EQ)) THEN
                        DES_POS_NEW(2,II) = DES_POS_NEW(2,II) + (TY2-BY1)
                END IF
             END IF
      ELSE IF(DES_PERIODIC_WALLS_Z) THEN
             TEMPD = DES_POS_NEW(3,L) - DES_POS_NEW(3,II)
             TEMPD = SQRT(TEMPD**2)
             IF(TEMPD.GT.(DES_RADIUS(L)+DES_RADIUS(II))) THEN
                IF(DES_POS_NEW(3,II).GE.(NZ2 - 2*RADIUS_EQ)) THEN
                   DES_POS_NEW(3,II) = DES_POS_NEW(3,II) - (NZ2-SZ1)
                ELSE IF(DES_POS_NEW(3,II).LE.(SZ1 + 2*RADIUS_EQ)) THEN
                        DES_POS_NEW(3,II) = DES_POS_NEW(3,II) + (NZ2-SZ1)
                END IF
             END IF
      END IF
      END IF
      
      DO K = 1, DIMN
         NORMOD=NORMOD+(DES_POS_NEW(K,L)-DES_POS_NEW(K,II))**2
      END DO

      NORMOD = SQRT(NORMOD)

      DO K = 1, DIMN
         IF(NORMOD.NE.0) THEN
            NORM(K)=(DES_POS_NEW(K,II)-DES_POS_NEW(K,L))/NORMOD
         ELSE 
            PRINT *,'NORMOD IS ZERO'
            STOP
         END IF
      END DO

      IF(DES_PERIODIC_WALLS) THEN
           DES_POS_NEW(1,II) = TEMPX
           DES_POS_NEW(2,II) = TEMPY
           DES_POS_NEW(3,II) = TEMPZ              
      END IF  

!      PRINT *,'NORMAL', NORM(1), NORM(2), NORM(3), NORMOD

      RETURN
      END SUBROUTINE CFNORMAL 


