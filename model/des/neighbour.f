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

      USE discretelement
      IMPLICIT NONE

      INTEGER I, II, PARTS

      DO I = 1, PARTICLES
         DO II = 1, MN
            NEIGHBOURS(II,I) = -1
         END DO
         NEIGHBOURS(1,I) = 0
      END DO

      IF (DO_NSQUARE) THEN
         CALL NSQUARE(PARTICLES)
      ELSE
         N2CT = 0.000
      END IF

!     IF (DO_NBS) THEN
!     CALL NBS(PARTICLES)
!     ELSE
!     NBSCT = 0.000
!     END IF

      IF (DO_QUADTREE) THEN
         CALL QUADTREE(PARTICLES)
      ELSE
         QUADCT = 0.000
      END IF

      IF (DO_OCTREE) THEN
      CALL OCTREE(PARTICLES)
      ELSE
      OCTCT = 0.000
      END IF

!     DO I = 1,PARTICLES
!     IF(NEIGHBOURS(1,I).NE.0) THEN
!     PRINT *,'N','->', I,':',(NEIGHBOURS(J,I), J=1,NEIGHBOURS(1,I)+1)
!     END IF
!     END DO

      RETURN
      END SUBROUTINE NEIGHBOUR


