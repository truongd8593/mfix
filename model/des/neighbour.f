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
      use desgrid
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
!         CALL GRID_BASED_NEIGHBOR_SEARCH
          call desgrid_neigh_build 
      ENDIF 
! pradeep setting do_nsearch false here is most appropriate 
      do_nsearch = .false.
!REMOVE PRADEEP 
!      call des_dbgneigh 
      RETURN
      END SUBROUTINE NEIGHBOUR


!------------------------------------------------------------------------
! subroutine       : des_dbgneigh 
! Author           : Pradeep G.
! Purpose          : For debugging the neighbour 
!------------------------------------------------------------------------
      subroutine des_dbgneigh

      use discretelement 
      use compar
      implicit none
! local variables 
      integer lparcount,lcurpar,lneighindx,lneigh,lneighcnt 
      character*30 filename 

      write(filename,'("dbg_neighbor",I4.4,".dat")') mype 
      open(44,file=filename)
      write(44,*) "Printing neighbour information"
      lparcount = 1 
      do lcurpar = 1,max_pip
         if(lparcount.gt.pip) exit
         if(pea(lcurpar,1)) cycle 
         lparcount = lparcount + 1
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
      end subroutine 
