!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAKE_ARRAYS_DES                                        C
!  Purpose: DES - allocating DES arrays                                C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!  Comments: Added some calls that are necessary if INTERPOLATION IS ON C
!  Reviewer: Tingwen Li                               Date: 23-Jan-08  C
!  Comments: Call cell_near_wall                                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE MAKE_ARRAYS_DES
      
      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      USE run
      IMPLICIT NONE
      
      INTEGER LN, K
      INTEGER CHECK_MPI
      INTEGER L, I, II
      DOUBLE PRECISION DIST, R_LM

      IF(COORDINATES == 'CYLINDRICAL') THEN
         WRITE (UNIT_LOG, *) ' '
         WRITE (UNIT_LOG, *) 'Cylindrical coordinates are being used. STOP'
         WRITE (UNIT_LOG, *) 'DES should only be run using cartesian coordinates.'
         WRITE (*, *) ' '
         WRITE (*, *) 'Cylindrical coordinates are being used. STOP'
         WRITE (*, *) 'DES should only be run using cartesian coordinates.'
         CALL MFIX_EXIT(myPE)
      END IF

      CHECK_MPI = NODESI * NODESJ * NODESK
      IF((CHECK_MPI.NE.1).AND.(DISCRETE_ELEMENT)) THEN
         WRITE (UNIT_LOG, *) ' '
         WRITE (UNIT_LOG, *) 'DES being run on multiple processors. STOP'
         WRITE (UNIT_LOG, *) 'DES should only be run serially on one processor.'
         WRITE (*, *) ' '
         WRITE (*, *) 'DES being run on multiple processors. STOP'
         WRITE (*, *) 'DES should only be run serially on one processor.'
         CALL MFIX_EXIT(myPE)
      END IF

      IF(DES_NEIGHBOR_SEARCH.EQ.UNDEFINED_I) THEN
         DES_NEIGHBOR_SEARCH = 1
         WRITE(*,*) 'Default N-Square search will be used'
      END IF

      IF(DES_NEIGHBOR_SEARCH.EQ.1) THEN
         DO_NSQUARE = .TRUE.
         WRITE(*,*) 'DEM USING N-SQUARE Search'
      ELSE IF(DES_NEIGHBOR_SEARCH.EQ.2) THEN
         DO_QUADTREE = .TRUE.
         WRITE(*,*) 'DEM USING QUADTREE Search'
      ELSE IF(DES_NEIGHBOR_SEARCH.EQ.3) THEN
         DO_OCTREE = .TRUE.
         WRITE(*,*) 'DEM USING OCTREE Search'
      ELSE IF(DES_NEIGHBOR_SEARCH.EQ.4) THEN
         DO_GRID_BASED_SEARCH = .TRUE.
         WRITE(*,*) 'DEM USING CELL LINKED SEARCH'
      END IF
      
      IF(RUN_TYPE == 'NEW') THEN ! Fresh run
         OPEN(UNIT=10, FILE='particle_input.dat', STATUS='OLD') 
         DO LN = 1, PARTICLES
            READ (10, *) (DES_POS_OLD(LN,K),K=1,DIMN),DES_RADIUS(LN),RO_Sol(LN),(DES_VEL_OLD(LN,K),K=1,DIMN)
            OMEGA_OLD(LN,:) = ZERO
            DES_POS_NEW(LN,:) = DES_POS_OLD(LN,:)
            DES_VEL_NEW(LN,:) = DES_VEL_OLD(LN,:)
            OMEGA_NEW(LN,:) = OMEGA_OLD(LN,:)
         END DO
      ELSE IF(RUN_TYPE == 'RESTART_1') THEN !  Read Restart
         CALL READ_DES_RESTART
         DES_POS_NEW(:,:) = DES_POS_OLD(:,:)
         DES_VEL_NEW(:,:) = DES_VEL_OLD(:,:)
         OMEGA_NEW(:,:) = OMEGA_OLD(:,:)
         DESRESDT = 0.0d0
         WRITE(*,*) 'DES_RES file read at Time= ', TIME
         WRITE(UNIT_LOG,*) 'DES_RES file read at Time= ', TIME
         IF(USE_COHESION) THEN
            WRITE(UNIT_LOG,*) 'Restart 1 is not implemented with DES-COHESION'
            WRITE(*,*) 'Restart 1 is not implemented with DES-COHESION'
            CALL MFIX_EXIT(myPE)
         END IF
         
      ELSE IF (RUN_TYPE == 'RESTART_2') THEN 
         WRITE(UNIT_LOG,*) 'Restart 2 is not implemented with DES'
         WRITE(*,*) 'Restart 2 is not implemented with DES'
         CALL MFIX_EXIT(myPE)
      END IF
      CALL CFASSIGN
      CALL PARTICLES_IN_CELL     

!     call cell_near_wall and set non_rect_bc
      if(NON_RECT_BC) call cell_near_wall

      RETURN
      END SUBROUTINE MAKE_ARRAYS_DES 
