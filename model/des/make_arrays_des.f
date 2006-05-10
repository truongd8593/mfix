!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAKE_ARRAYS_DES(PARTS)                                 C
!  Purpose: DES - allocating DES arrays                                C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE MAKE_ARRAYS_DES(PARTS)
      
      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      IMPLICIT NONE
      
      INTEGER LN, K, PARTS
      INTEGER CHECK_MPI

       IF(COORDINATES == 'CYLINDRICAL') THEN
          WRITE (UNIT_LOG, *) ' '
          WRITE (UNIT_LOG, *) 'Cylindrical coordinates are being used. STOP'
          PRINT *,'DES should only be run using cartesian coordinates.'
          STOP
       END IF

       CHECK_MPI = NODESI * NODESJ * NODESK
       IF((CHECK_MPI.NE.1).AND.(DISCRETE_ELEMENT)) THEN
          WRITE (UNIT_LOG, *) ' '
          WRITE (UNIT_LOG, *) 'DES being run on multiple processors. STOP'
          PRINT *,'DES should only be run serially on one processor.'
          STOP
       END IF

       IF(DES_NEIGHBOR_SEARCH.EQ.UNDEFINED_I) THEN
          DES_NEIGHBOR_SEARCH = 1
          PRINT *,'Default N-Square search will be implemented'
       END IF
       IF(DES_NEIGHBOR_SEARCH.EQ.1) THEN
          DO_NSQUARE = .TRUE.
          PRINT *,'N-SQUARE Search'
       ELSE IF(DES_NEIGHBOR_SEARCH.EQ.2) THEN
          DO_QUADTREE = .TRUE.
          PRINT *,'QUADTREE Search'
       ELSE IF(DES_NEIGHBOR_SEARCH.EQ.3) THEN
          DO_OCTREE = .TRUE.
          PRINT *,'OCTREE Search'
       END IF
 
        OPEN(UNIT=10, FILE='particle_input.dat', STATUS='OLD')
        DO LN = 1, PARTICLES
          READ (10, *) (DES_POS_OLD(LN,K),K=1,DIMN),DES_RADIUS(LN),RO_Sol(LN),(DES_VEL_OLD(LN,K),K=1,DIMN)
          OMEGA_OLD(LN,:) = ZERO
          DES_POS_NEW(LN,:) = DES_POS_OLD(LN,:)
          DES_VEL_NEW(LN,:) = DES_VEL_OLD(LN,:)
          OMEGA_NEW(LN,:) = OMEGA_OLD(LN,:)
        END DO

      RETURN
      END SUBROUTINE MAKE_ARRAYS_DES 


