!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NEIGHBOUR(PARTS)                                       C
!  Purpose: DES - Neighbors search; N-Square, Quadtree(2D)/Octree(3D)  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE NEIGHBOUR(PARTS)

      USE param1
      USE discretelement
      IMPLICIT NONE

      INTEGER I, II, PARTS, LL, CO, NI, TEMP

      DO I = 1, PARTICLES
         PPOS(I,:) = DES_POS_NEW(I,:)
         DO II = 1, MAXNEIGHBORS
            NEIGHBOURS(I,II) = -1
         END DO
         NEIGHBOURS(I,1) = 0
      END DO

      IF (DO_NSQUARE) THEN
         CALL NSQUARE(PARTICLES)
      ELSE
         N2CT = ZERO
      END IF

      IF (DO_QUADTREE) THEN
         CALL QUADTREE(PARTICLES)
      ELSE
         QUADCT = ZERO
      END IF

      IF (DO_OCTREE) THEN
      CALL OCTREE(PARTICLES)
      ELSE
      OCTCT = ZERO
      END IF

      RETURN
      END SUBROUTINE NEIGHBOUR


