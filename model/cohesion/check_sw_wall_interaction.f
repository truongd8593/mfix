!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_SW_WALL_INTERACTION(I,J)                         C
!  Purpose: This module will check for interactions between particles  C
!           and the wall when the square-well cohesion model is being  C 
!           used                                                       C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      SUBROUTINE CHECK_SW_WALL_INTERACTION(I,J)

!-----MODULES USED
      USE discretelement
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER I,J,N
      INTEGER LINK
      INTEGER CHECK_LINK
      DOUBLE PRECISION RADIUS


      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START CHECK SQUARE WELL WALL INTERACTION'
      END IF

      IF(J.gt.PARTICLES+DIMN*2)THEN
          STOP
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
        (DES_POS_NEW(2,J)-DES_POS_NEW(2,I))**2)

         LINK = CHECK_LINK(I,J)
         IF(LINK.eq.1)THEN
            CALL LINKED_INTERACTION_EVAL(I,J)
         ELSE
            CALL UNLINKED_INTERACTION_EVAL(I,J)
         END IF


      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END CHECK SQUARE WELL WALL INTERACTION'
      END IF

      END SUBROUTINE CHECK_SW_WALL_INTERACTION
