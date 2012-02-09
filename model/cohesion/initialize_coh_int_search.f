!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INITIALIZE_PARTICLE_SEARCH                             C
!  Purpose: Module to initialize parameters used in recognizing        C
!           particle-particle interactions                             C
!                                                                      C
!      Search grids are sized so that only one particle can            C
!      fit in each grid                                                C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      SUBROUTINE INITIALIZE_COH_INT_SEARCH


!-----MODULES USED
      USE discretelement

!-----LOCAL DUMMY VARIABLES
      INTEGER i,j,k,m
      INTEGER grid(2)
      DOUBLE PRECISION largest_radius
      DOUBLE PRECISION largest_dimension

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START INITIALIZE PARTICLE SEARCH'
      END IF
     
!-----Find largest dimension for search grid size
      largest_radius=0
      IF(SQUARE_WELL)THEN
         DO i=1,PARTICLES
             IF(WELL_WIDTH(i).gt.largest_radius)THEN
               largest_radius=WELL_WIDTH(i)
             END IF
         END DO
      ELSE
         DO i=1,PARTICLES
           IF(DES_RADIUS(i).gt.largest_radius)THEN
              largest_radius=DES_RADIUS(i)
           END IF
         END DO
      END IF

      IF(VAN_DER_WAALS)THEN
         largest_radius=0.5*(VDW_OUTER_CUTOFF+2*largest_radius)
      END IF
     
      largest_dimension=2.0*largest_radius

      IF(COHESION_DEBUG.gt.2)THEN
         PRINT *,'****Largest dimension=',largest_dimension
      END IF

!-----Determine search grid size
      SEARCH_GRIDS(1)=INT((EX2-WX1)/largest_dimension)
      SEARCH_GRIDS(2)=INT((TY2-BY1)/largest_dimension)
      SEARCH_GRIDS(3)=INT((NZ2-SZ1)/largest_dimension)


      IF(SEARCH_GRIDS(1).gt.100)THEN
        SEARCH_GRIDS(1)=100
      END IF

      IF(SEARCH_GRIDS(2).gt.100)THEN
        SEARCH_GRIDS(2)=100
      END IF

      IF(SEARCH_GRIDS(3).gt.100)THEN
        SEARCH_GRIDS(3)=100
      END IF

      SEARCH_GRID_SIZE(1)=(EX2-WX1)/( REAL(SEARCH_GRIDS(1)) )
      SEARCH_GRID_SIZE(2)=(TY2-BY1)/( REAL(SEARCH_GRIDS(2)) )
      SEARCH_GRID_SIZE(3)=(NZ2-SZ1)/( REAL(SEARCH_GRIDS(3)) )


      IF(COHESION_DEBUG.gt.2)THEN
         PRINT *,'****Search Grids=',SEARCH_GRIDS(1),SEARCH_GRIDS(2)
         PRINT *,'****Search Grid Size=',SEARCH_GRID_SIZE(1),SEARCH_GRID_SIZE(2)
      END IF

!-----Initialize all searching parameters
      IF(COHESION_DEBUG.gt.2)THEN
         PRINT *,'****Initialize search grids...'
      END IF

      DO i=1,SEARCH_GRIDS(1)
         DO j=1,SEARCH_GRIDS(2)
            DO K=1,SEARCH_GRIDS(3)
               DO m=1,MAX_PART_IN_GRID
                  IF(m.eq.1)THEN
                     PART_IN_GRID(i,j,k,m)=0  !!first slot holds count of parts in grid
                  ELSE
                     PART_IN_GRID(i,j,k,m)=-1 !!other slots hold particle #
                  END IF
               END DO
            END DO
         END DO
      END DO

!-----Put particles in grids for first time
      IF(COHESION_DEBUG.gt.2)THEN
         PRINT *,'****Place particles in grids...'
      END IF

      DO i=1,PARTICLES
        !!Add particle to new grid
        DO j=1,DIMN
           PART_GRID(i,j)=INT(DES_POS_NEW(i,j)/SEARCH_GRID_SIZE(j))+1
        END DO 
        IF(DIMN.eq.2)THEN
           PART_GRID(i,3)=1
        END IF
        PART_IN_GRID(PART_GRID(i,1),PART_GRID(i,2),PART_GRID(i,3)&
            ,PART_IN_GRID(PART_GRID(i,1),PART_GRID(i,2),PART_GRID(i,3),1)+2)=i
        PART_GRID(i,4)=PART_IN_GRID(PART_GRID(i,1),PART_GRID(i,2),PART_GRID(i,3),1)+2
        PART_IN_GRID(PART_GRID(i,1),PART_GRID(i,2),PART_GRID(i,3),1)=&
           PART_IN_GRID(PART_GRID(i,1),PART_GRID(i,2),PART_GRID(i,3),1)+1
      END DO

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END INITIALIZE PARTICLE SEARCH'
      END IF


      END SUBROUTINE INITIALIZE_COH_INT_SEARCH 
