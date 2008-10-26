!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NEIGHBOUR                                              C
!  Purpose: DES - Neighbors search; N-Square, Quadtree(2D)/Octree(3D)  C
!           Now also Cell linked search Aug 07                         C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!  Comments: Added cell linked-list search                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE NEIGHBOUR

      USE param1
      USE discretelement
      IMPLICIT NONE

      INTEGER I, II, LL, CO, NI, TEMP
      

      DO I = 1, PARTICLES
         PPOS(I,:) = DES_POS_NEW(I,:)
         DO II = 1, MAXNEIGHBORS
            NEIGHBOURS(I,II) = -1
         END DO
         NEIGHBOURS(I,1) = 0
      END DO

      IF (DO_NSQUARE) THEN
         CALL NSQUARE
      ELSE
         N2CT = ZERO
      END IF

      IF (DO_QUADTREE) THEN
         CALL QUADTREE
      ELSE
         QUADCT = ZERO
      END IF

      IF (DO_OCTREE) THEN
         CALL OCTREE
      ELSE
         OCTCT = ZERO
      END IF
      IF (DO_GRID_BASED_SEARCH) THEN 
         CALL GRID_BASED_NEIGHBOR_SEARCH
      ENDIF 
      
      RETURN
      END SUBROUTINE NEIGHBOUR
