!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: UPDATE_SEARCH_GRIDS                                    C
!  Purpose: This module updates the search grids used to identify      C
!           cohesion particle interactions                             C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE UPDATE_SEARCH_GRIDS

!-----MODULES USED
      USE discretelement
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER i,j
      INTEGER new(3)
      LOGICAL ALREADY_EXISTS

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START UPDATE SEARCH GRIDS'
      END IF

!-----START LOOP OVER ALL PARTICLES
      DO i=1,PARTICLES

        DO j=1,DIMN
           new(j)=INT(DES_POS_NEW(j,i)/SEARCH_GRID_SIZE(j))+1
        END DO
 
        IF(DIMN.eq.2)THEN
           new(3)=1
        END IF
        
        !!Did particle change grids?
        IF((new(1).ne.PART_GRID(1,i)).OR.(new(2).ne.PART_GRID(2,i))&
                 .OR.(new(3).ne.PART_GRID(3,i)))THEN

      IF(COHESION_DEBUG.gt.2)THEN
         PRINT *,'***Particle=',i
         PRINT *,'*****Old i j=',PART_GRID(1,i), PART_GRID(2,i)
         PRINT *,'*****New i j=',new(1), new(j)
      END IF

  !!Remove particle from old grid
        !!Check if particle is the last in the list for the current grid
        IF(PART_GRID(4,i).lt.(PART_IN_GRID(PART_GRID(1,i),PART_GRID(2,i),PART_GRID(3,i),1)+1))THEN
           j=PART_GRID(4,i)
           DO WHILE(j.le.(PART_IN_GRID(PART_GRID(1,i),PART_GRID(2,i),PART_GRID(3,i),1)+1))

             !!Shift particles to fill in vacant slot in grid
             PART_IN_GRID(PART_GRID(1,i),PART_GRID(2,i),PART_GRID(3,i),j)=&
                PART_IN_GRID(PART_GRID(1,i),PART_GRID(2,i),PART_GRID(3,i),j+1)

             !!Update grid slot recorded with each particle
             IF(PART_IN_GRID(PART_GRID(1,i),PART_GRID(2,i),PART_GRID(3,i),j+1).gt.0)THEN
                PART_GRID(4,PART_IN_GRID(PART_GRID(1,i),PART_GRID(2,i),PART_GRID(3,i),j+1))=j
             END IF

             j=j+1
           END DO  
        ELSE
           PART_IN_GRID(PART_GRID(1,i),PART_GRID(2,i),PART_GRID(3,i),PART_GRID(4,i))=-1
        END IF

        !!Update count of particles in grid
        PART_IN_GRID(PART_GRID(1,i),PART_GRID(2,i),PART_GRID(3,i),1)=&
           PART_IN_GRID(PART_GRID(1,i),PART_GRID(2,i),PART_GRID(3,i),1)-1



   !!Add particle to new grid
        DO j=1,DIMN
           PART_GRID(j,i)=new(j)
        END DO

        !!place particle in last position in new grid   
        PART_IN_GRID(new(1),new(2),new(3)&
            ,PART_IN_GRID(new(1),new(2),new(3),1)+2)=i

        !!record particle's position in new grid with it's own info array
        PART_GRID(4,i)=PART_IN_GRID(new(1),new(2),new(3),1)+2

        !!Update count of particles in new grid
        PART_IN_GRID(new(1),new(2),new(3),1)=PART_IN_GRID(new(1),new(2),new(3),1)+1

        END IF  !!Did particle move out of current grid?  

      END DO  !!Loop over all particles

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END UPDATE SEARCH GRIDS'
      END IF

      END SUBROUTINE UPDATE_SEARCH_GRIDS
