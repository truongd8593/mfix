!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NEIGHBOUR                                              C
!  Purpose: DES - Neighbors search;
!           N-Square,
!           Quadtree(2D)/Octree(3D)  (use at own risk)
!           Cell linked
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE NEIGHBOUR

!-----------------------------------------------
! Modules
!-----------------------------------------------      
      USE param1
      USE discretelement
      use desgrid
      Use des_thermo      
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------

! Reset PPOS and NEIGHBOURS back to initialized values
      PPOS(:,:) = DES_POS_NEW(:,:)
      NEIGHBOURS(:,:) = -1
      NEIGHBOURS(:,1) = 0
      
      IF (DES_NEIGHBOR_SEARCH.EQ.1) THEN
         CALL NSQUARE
      ELSEIF (DES_NEIGHBOR_SEARCH.EQ.4) THEN 
          CALL DESGRID_NEIGH_BUILD
      ENDIF 

! resetting do_nsearch to false here since neighbor search will have
! just been invoked
      DO_NSEARCH = .FALSE.

      RETURN
      END SUBROUTINE NEIGHBOUR


!------------------------------------------------------------------------
! subroutine       : des_dbgneigh 
! Author           : Pradeep G.
! Purpose          : For debugging the neighbour 
!                  : call des_dbgneigh at end of neighbour routine
!------------------------------------------------------------------------
      subroutine des_dbgneigh

!-----------------------------------------------
! Modules
!-----------------------------------------------       
      use discretelement 
      use compar
      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------  
! particle count (index to track accounted for particles)
      INTEGER :: lparcount
! index of l particle
      INTEGER :: lcurpar
! number of neighbors of l particle (lcurpart)
      INTEGER :: LNEIGHCNT
! index of l particle neighbor
      INTEGER :: LNEIGH
! local counter of number of neighbors of l particle
      INTEGER :: LNEIGHINDX
      character*30 filename 
!----------------------------------------------- 

      write(filename,'("dbg_neighbor",I4.4,".dat")') mype 
      open(44,file=filename)
      write(44,*) "Printing neighbour information"
      
      lparcount = 1 
! looping over the maximum allowed particles in this processor      
      do lcurpar = 1,max_pip
! exiting loop if reached the existing number of paarticles in this processor
         if(lparcount.gt.pip) exit
! cycling if this particle index does not exist
         if(.not.pea(lcurpar,1)) cycle 
! incrementing the local particle count on 'existing particles'
         lparcount = lparcount + 1
! cycling if this particle index is ghost 
         if(pea(lcurpar,4)) cycle 
         lneighcnt = neighbours(lcurpar,1) 
         if(lneighcnt .gt. 0) then 
            write(44,*)"----------------------------------------------------"
            write(44,*)"Particle ", iglobal_id(lcurpar), "Has " , lneighcnt, " neighbours"
            write(44,*)"NeighID position"
            write(44,*)"----------------------------------------------------"
         end if 
         do lneighindx=2,lneighcnt+1 
            lneigh = neighbours(lcurpar,lneighindx)
            write(44,'(I4.4,2(2x,F12.4))') iglobal_id(lneigh),des_pos_new(lneigh,1),des_pos_new(lneigh,2)  
         end do 
      end do 
      close(44)

      end subroutine DES_DBGNEIGH
