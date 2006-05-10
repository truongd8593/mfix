!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INITIALIZE_COHESION_PARAMETERS                         C
!  Purpose: This module gives initial values to all                    C
!           cohesion parameters                                        C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      SUBROUTINE INITIALIZE_COHESION_PARAMETERS

!-----MODULES USED
      USE discretelement

!-----LOCAL DUMMY VARIABLES
      INTEGER i,j

      IF(COHESION_DEBUG.gt.0)THEN
        PRINT *, '**START INITIALIZE COHESION PARAMETERS'
      END IF

    
!-----INITIALIZATIONS
      DO i=1, PARTICLES+2*DIMN
         IS_LINKED(i)=0
         LINKS(i,1)=0
         DO j=2,MAXNEIGHBORS
           LINKS(i,j)=-1
         END DO
         WELL_WIDTH(i)=RADIUS_RATIO*DES_RADIUS(i)
         WELL_DEPTH(i)=MASTER_WELL_DEPTH
         IF(I.gt.PARTICLES)THEN ! Walls
            WELL_WIDTH(i)=WALL_RADIUS_RATIO*DES_RADIUS(1)
            WELL_DEPTH(i)=MASTER_WALL_WELL_DEPTH
            PMASS(i)=99999999
         END IF
      END DO

      !!surface energy set so that force stays constant at inner cut off
      SURFACE_ENERGY=HAMAKER_CONSTANT/&
             (24*3.14*VDW_INNER_CUTOFF*VDW_INNER_CUTOFF)
      WALL_SURFACE_ENERGY=WALL_HAMAKER_CONSTANT/&
             (24*3.14*WALL_VDW_INNER_CUTOFF*WALL_VDW_INNER_CUTOFF)

      MAX_PART_IN_GRID=3

      IF(COHESION_DEBUG.gt.0)THEN
        PRINT *, '**END INITIALIZE COHESION PARAMETERS'
      END IF

      

      END SUBROUTINE INITIALIZE_COHESION_PARAMETERS
