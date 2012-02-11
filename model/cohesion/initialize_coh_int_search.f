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

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: i,j,k,m
      INTEGER :: i_ind,j_ind,k_ind,np_ijk, NP
      DOUBLE PRECISION :: largest_radius
      DOUBLE PRECISION :: largest_dimension
!-----------------------------------------------

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START INITIALIZE PARTICLE SEARCH'
      ENDIF
     
!-----Find largest dimension for search grid size
      largest_radius=0
      IF(SQUARE_WELL)THEN
         DO NP=1,PARTICLES
             IF(WELL_WIDTH(NP).gt.largest_radius)THEN
               largest_radius=WELL_WIDTH(NP)
             ENDIF
         ENDDO
      ELSE
         DO NP=1,PARTICLES
           IF(DES_RADIUS(NP).gt.largest_radius)THEN
              largest_radius=DES_RADIUS(NP)
           ENDIF
         ENDDO
      ENDIF

      IF(VAN_DER_WAALS)THEN
         largest_radius=0.5*(VDW_OUTER_CUTOFF+2.d0*largest_radius)
      ENDIF
     
      largest_dimension=2.0*largest_radius

      IF(COHESION_DEBUG.gt.2)THEN
         PRINT *,'****Largest dimension=',largest_dimension
      ENDIF

!-----Determine search grid size
      SEARCH_GRIDS(1)=INT((EX2-WX1)/largest_dimension)
      SEARCH_GRIDS(2)=INT((TY2-BY1)/largest_dimension)
      SEARCH_GRIDS(3)=INT((NZ2-SZ1)/largest_dimension)

! cap size of search grids       
      IF(SEARCH_GRIDS(1).gt.100)THEN
        SEARCH_GRIDS(1)=100
      ENDIF
      IF(SEARCH_GRIDS(2).gt.100)THEN
        SEARCH_GRIDS(2)=100
      ENDIF
      IF(SEARCH_GRIDS(3).gt.100)THEN
        SEARCH_GRIDS(3)=100
      ENDIF

      SEARCH_GRID_SIZE(1)=(EX2-WX1)/( REAL(SEARCH_GRIDS(1)) )
      SEARCH_GRID_SIZE(2)=(TY2-BY1)/( REAL(SEARCH_GRIDS(2)) )
      SEARCH_GRID_SIZE(3)=(NZ2-SZ1)/( REAL(SEARCH_GRIDS(3)) )

      IF(COHESION_DEBUG.gt.2)THEN
         PRINT *,'****Search Grids=',SEARCH_GRIDS(1),SEARCH_GRIDS(2)
         PRINT *,'****Search Grid Size=',SEARCH_GRID_SIZE(1),SEARCH_GRID_SIZE(2)
      ENDIF

!-----Initialize all searching parameters
      IF(COHESION_DEBUG.gt.2)THEN
         PRINT *,'****Initialize search grids...'
      ENDIF

      DO i=1,SEARCH_GRIDS(1)
         DO j=1,SEARCH_GRIDS(2)
            DO K=1,SEARCH_GRIDS(3)
               DO m=1,MAX_PART_IN_GRID
                  IF(m.eq.1)THEN
                     PART_IN_GRID(i,j,k,m)=0  !!first slot holds count of parts in grid
                  ELSE
                     PART_IN_GRID(i,j,k,m)=-1 !!other slots hold particle #
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

!-----Put particles in grids for first time
      IF(COHESION_DEBUG.gt.2)THEN
         PRINT *,'****Place particles in grids...'
      END IF

      DO NP=1,PARTICLES
! Add particle to new grid
! part_grid stores particles i, j, k indices of grid      
         DO j=1,DIMN
            PART_GRID(NP,j)=INT(DES_POS_NEW(NP,j)/SEARCH_GRID_SIZE(j))+1
         ENDDO 
         IF(DIMN.eq.2)THEN
            PART_GRID(NP,3)=1
         ENDIF

         i_ind = PART_GRID(NP,1)
         j_ind = PART_GRID(NP,2)
         k_ind = PART_GRID(NP,3)
! number of particles at this ijk location on grid (initially 0)       
         np_ijk = PART_IN_GRID(i_ind,j_ind,k_ind,1)

! store particle index at part_in_grid(i,j,k,>1) position        
         PART_IN_GRID(i_ind,j_ind,k_ind,np_ijk+2)=NP

! overwrite part_grid(np,4) with number of particles in grid+1?
         PART_GRID(NP,4)=np_ijk+2

! increment number of particles at ijk location on grid 
         PART_IN_GRID(i_ind,j_ind,k_ind,1)=np_ijk+1
! the grid size is generally based on value barely larger than
! the particle diameter therefore it will be unlikely that the
! number of particles at a given grid location will be greater
! than 1.  one possible exception is if the number of grids called 
! for is greater than 100 and the grid is capped at 100

      ENDDO

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END INITIALIZE PARTICLE SEARCH'
      ENDIF

      END SUBROUTINE INITIALIZE_COH_INT_SEARCH 
