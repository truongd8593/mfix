!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADD_PART_TO_LINK_LIST(I,J)                             C
!  Purpose: Module to Add particles to their                           C
!           respective linked lists ("linked" particles have           C
!           overlapping square-wells)                                  C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C



      SUBROUTINE ADD_PART_TO_LINK_LIST(I,J)


!-----MODULES USED
      USE discretelement
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER I,J,JJ,II
      LOGICAL ALREADY_EXISTS

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START ADD LINK'
      END IF 

      IF(I.gt.PARTICLES+DIMN*2)THEN 
        PRINT *, 'STOPPED FROM inside add_part_to_link_list @',CALLED
        PRINT *, 'I=',I
        PRINT *, 'J=',J
        STOP
      END IF 

      IF(J.gt.PARTICLES+DIMN*2)THEN 
        PRINT *, 'STOPPED FROM inside add_part_to_link_list @',CALLED
        PRINT *, 'I=',I
        PRINT *, 'J=',J
        STOP
      END IF 


!-----Add particle J to the list of particle I
      IF(I.le.PARTICLES)THEN
         LINKS(LINKS(1,I)+2,I)=J
         LINKS(1,I)=LINKS(1,I)+1
         IS_LINKED(I)=1
      END IF

!-----Add particle I to the list of particle J
      IF(J.le.PARTICLES)THEN
         LINKS(LINKS(1,J)+2,J)=I
         LINKS(1,J)=LINKS(1,J)+1
         IS_LINKED(J)=1
      END IF


      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END ADD LINK'
      END IF 

      END SUBROUTINE ADD_PART_TO_LINK_LIST