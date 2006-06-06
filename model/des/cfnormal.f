!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNORMAL(L, II, J, NORM)                               C
!  Purpose: DES - Calculate the normal vector between particles        C
!           in collision                                               C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments: Normal for Eqn 6  from the following paper                C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNORMAL(L, II, J, NORM) 

      USE param1      
      USE discretelement
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER L, II, J, K
      DOUBLE PRECISION NORMOD, NORM(DIMN), DIST(DIMN)
      DOUBLE PRECISION TEMPX, TEMPY, TEMPZ, TEMPD

!----------------------------------------------------------------------------

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
      
      DIST(:) = DES_POS_NEW(II,:) - DES_POS_NEW(L,:)
      NORMOD = SQRT(DES_DOTPRDCT(DIST,DIST))
      IF(NORMOD.NE.ZERO) THEN
!         NORM(:)= DIST(:)/PN_DIST(L,J)
         NORM(:)= DIST(:)/NORMOD
      ELSE 
         PRINT *,'NORMOD IS ZERO', II,L
         STOP
      END IF

      IF(DES_PERIODIC_WALLS) THEN
           DES_POS_NEW(II,1) = TEMPX
           DES_POS_NEW(II,2) = TEMPY
           IF (DIMN.EQ.3) DES_POS_NEW(II,3) = TEMPZ              
      END IF  

      RETURN
      END SUBROUTINE CFNORMAL 


