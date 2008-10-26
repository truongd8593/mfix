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
      USE geometry
      
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
         TEMPD = ABS(DES_POS_NEW(L,1) - DES_POS_NEW(II,1))
         IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_X) THEN 
            IF(DEBUG_DES) PRINT*, 'From cfnormal.f'
            IF(DEBUG_DES) PRINT*,'II, L = ', II,L
            IF(DEBUG_DES) PRINT*,'OLD POS, II, L = ', DES_POS_NEW(II,1), DES_POS_NEW(L,1)
            IF(TEMPX.GT.DES_POS_NEW(L,1)) THEN 
               DES_POS_NEW(II,1) = DES_POS_NEW(II,1) - XE(IMAX1)
               IF(DEBUG_DES) PRINT*,'NEW POS WEST= ', DES_POS_NEW(II,1)
            ELSE
               DES_POS_NEW(II,1) = DES_POS_NEW(II,1) + XE(IMAX1)
               
               IF(DEBUG_DES) PRINT*,'NEW POS EAST = ', DES_POS_NEW(II,1)
            ENDIF
         ENDIF
         
         TEMPD = ABS(DES_POS_NEW(L,2) - DES_POS_NEW(II,2))
         IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_Y) THEN 
            IF(TEMPY.GT.DES_POS_NEW(L,2)) THEN 
               DES_POS_NEW(II,2) = DES_POS_NEW(II,2) - YN(JMAX1)
            ELSE
               DES_POS_NEW(II,2) = DES_POS_NEW(II,2) + YN(JMAX1)
            ENDIF
         ENDIF
         
         IF(DIMN.EQ.3) THEN
            TEMPD = ABS(DES_POS_NEW(L,3) - DES_POS_NEW(II,3))
            IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_Z) THEN 
               IF(TEMPZ.GT.DES_POS_NEW(L,3)) THEN 
                  DES_POS_NEW(II,3) = DES_POS_NEW(II,3) - ZT(KMAX1)
               ELSE
                  DES_POS_NEW(II,3) = DES_POS_NEW(II,3) + ZT(KMAX1)
               ENDIF
            ENDIF
         ENDIF
         
      END IF
      
      DIST(:) = DES_POS_NEW(II,:) - DES_POS_NEW(L,:)
      NORMOD = SQRT(DES_DOTPRDCT(DIST,DIST))
      IF(NORMOD.NE.ZERO) THEN
!     NORM(:)= DIST(:)/PN_DIST(L,J)
         NORM(:)= DIST(:)/NORMOD
      ELSE 
         PRINT*, 'From cfnormal.f'
         PRINT *,'NORMOD IS ZERO', II,L
         PRINT *,'NORMOD IS ZERO', DES_POS_NEW(II,:), DES_POS_NEW(L,:) 
         STOP
      END IF

      IF(DES_PERIODIC_WALLS) THEN
         DES_POS_NEW(II,1) = TEMPX
         DES_POS_NEW(II,2) = TEMPY
         IF (DIMN.EQ.3) DES_POS_NEW(II,3) = TEMPZ              
      END IF  

      RETURN
      END SUBROUTINE CFNORMAL 


