!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_VDW_WALL_INTERACTION(I,J)                        C
!  Purpose: This module calculates the attractive force between        C
!           a wall and a particle using the Hamaker van der Waals modelC
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      SUBROUTINE CHECK_VDW_WALL_INTERACTION(I,J)

!-----MODULES USED
      USE discretelement
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER I,J,N,K, II
      INTEGER LINK
      INTEGER CHECK_LINK
      DOUBLE PRECISION RADIUS
      DOUBLE PRECISION FORCE, RELPOS(2)
      DOUBLE PRECISION DIST
      LOGICAL ALREADY_EXISTS

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START CHECK VAN DER WAALS WALL INTERACTION'
      END IF

      IF(J.eq.PARTICLES+1)THEN !West wall
         DES_POS_NEW(1,J)=WX1
         DES_POS_NEW(2,J)=DES_POS_NEW(2,I)
         DES_POS_NEW(3,J)=DES_POS_NEW(3,I)
         DES_VEL_NEW(1,J)=0.0
         DES_VEL_NEW(2,J)=0.0
         DES_VEL_NEW(3,J)=0.0
      END IF         

      IF(J.eq.PARTICLES+2)THEN !Bottom wall
         DES_POS_NEW(1,J)=DES_POS_NEW(1,I)
         DES_POS_NEW(2,J)=BY1
         DES_POS_NEW(3,J)=DES_POS_NEW(3,I)
         DES_VEL_NEW(1,J)=0.0
         DES_VEL_NEW(2,J)=0.0
         DES_VEL_NEW(3,J)=0.0
      END IF  

      IF(J.eq.PARTICLES+3)THEN !East wall
         DES_POS_NEW(1,J)=EX2
         DES_POS_NEW(2,J)=DES_POS_NEW(2,I)
         DES_POS_NEW(3,J)=DES_POS_NEW(3,I)
         DES_VEL_NEW(1,J)=0.0
         DES_VEL_NEW(2,J)=0.0
         DES_VEL_NEW(3,J)=0.0
      END IF  

      IF(J.eq.PARTICLES+4)THEN !Top wall
         DES_POS_NEW(1,J)=DES_POS_NEW(1,I)
         DES_POS_NEW(2,J)=TY2
         DES_POS_NEW(3,J)=DES_POS_NEW(3,I)
         DES_VEL_NEW(1,J)=0.0
         DES_VEL_NEW(2,J)=0.0
         DES_VEL_NEW(3,J)=0.0
      END IF  

      IF(J.eq.PARTICLES+5)THEN !North wall
         DES_POS_NEW(1,J)=DES_POS_NEW(1,I)
         DES_POS_NEW(2,J)=DES_POS_NEW(2,I)
         DES_POS_NEW(3,J)=NZ2
         DES_VEL_NEW(1,J)=0.0
         DES_VEL_NEW(2,J)=0.0
         DES_VEL_NEW(3,J)=0.0
      END IF 

      IF(J.eq.PARTICLES+6)THEN !South wall
         DES_POS_NEW(1,J)=DES_POS_NEW(1,I)
         DES_POS_NEW(2,J)=DES_POS_NEW(2,I)
         DES_POS_NEW(3,J)=SZ1
         DES_VEL_NEW(1,J)=0.0
         DES_VEL_NEW(2,J)=0.0
         DES_VEL_NEW(3,J)=0.0
      END IF 

      RADIUS=SQRT((DES_POS_NEW(1,J)-DES_POS_NEW(1,I))**2+&
        (DES_POS_NEW(2,J)-DES_POS_NEW(2,I))**2+&
        (DES_POS_NEW(3,J)-DES_POS_NEW(3,I))**2)


      DIST=RADIUS-DES_RADIUS(I)

      IF(DIST.lt.WALL_VDW_OUTER_CUTOFF)THEN
         DO II=1,2
            RELPOS(II)=DES_POS_NEW(II,J)-DES_POS_NEW(II,I)
         END DO
         IF(DIST.gt.WALL_VDW_INNER_CUTOFF)THEN
            FORCE=WALL_HAMAKER_CONSTANT*DES_RADIUS(I)/(6*DIST*DIST)
         ELSE
            FORCE=4*3.14*WALL_SURFACE_ENERGY*DES_RADIUS(I)
         END IF !Long range or surface?

         DO K=1,2
           FC(K,I)=FC(K,I)+RELPOS(K)/RADIUS*FORCE
           FC(K,J)=FC(K,J)-RELPOS(K)/RADIUS*FORCE
         END DO 
                    
      END IF !Is particle within cutoff?


      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END CHECK VAN DER WAALS WALL INTERACTION'
      END IF

      END SUBROUTINE CHECK_VDW_WALL_INTERACTION 
