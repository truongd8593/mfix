!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NEIGHBOUR                                              C
!>  Purpose: DES - Neighbors search; N-Square, Quadtree(2D)/Octree(3D)  
!>           Now also Cell linked search Aug 07                         
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


! J.Musser
! Reset PPOS and NEIGHBOURS back to initialized values
      PPOS(:,:) = DES_POS_NEW(:,:)
      NEIGHBOURS(:,:) = -1
      NEIGHBOURS(:,1) = 0

      IF (DES_NEIGHBOR_SEARCH.EQ.1) THEN
         CALL NSQUARE
      ELSE
         N2CT = ZERO
      ENDIF

      IF (DES_NEIGHBOR_SEARCH.EQ.2) THEN
         CALL QUADTREE
      ELSE
         QUADCT = ZERO
      ENDIF

      IF (DES_NEIGHBOR_SEARCH.EQ.3) THEN
         CALL OCTREE
      ELSE
         OCTCT = ZERO
      ENDIF

      IF (DES_NEIGHBOR_SEARCH.EQ.4) THEN 
         CALL GRID_BASED_NEIGHBOR_SEARCH
      ENDIF 
      
      RETURN
      END SUBROUTINE NEIGHBOUR
