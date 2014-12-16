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

INTEGER :: cc,dd,ii,jj,iii,jjj,ddd

!-----------------------------------------------
! Reset PPOS and NEIGHBOURS back to initialized values
      PPOS(:,:) = DES_POS_NEW(:,:)

      collisions_old = collisions
      pv_coll_old = pv_coll
      pfn_coll_old = pfn_coll
      pft_coll_old = pft_coll

      old_collision_num = collision_num
      collision_num = 0

      IF (DES_NEIGHBOR_SEARCH.EQ.1) THEN
         CALL NSQUARE
      ELSEIF (DES_NEIGHBOR_SEARCH.EQ.4) THEN
          CALL DESGRID_NEIGH_BUILD
      ENDIF

      dd = 1
      iii = collisions_old(1,dd)
      jjj = collisions_old(2,dd)

      do cc = 1, collision_num
         ii = collisions(1,cc)
         jj = collisions(2,cc)

         iii = collisions_old(1,dd)
         jjj = collisions_old(2,dd)

         do while (dd .le. old_collision_num .and. (iii < ii))
            dd = dd + 1
            iii = collisions_old(1,dd)
            jjj = collisions_old(2,dd)
         enddo

         ddd = dd
         do while (ddd .le. old_collision_num .and. iii.eq.ii .and. jjj.ne.jj )
            ddd = ddd + 1
            iii = collisions_old(1,ddd)
            jjj = collisions_old(2,ddd)
         enddo

         if (ii.eq.iii .and. jj.eq.jjj) then
            pv_coll(cc) = pv_coll_old(ddd)
            pfn_coll(:,cc) = pfn_coll_old(:,ddd)
            pft_coll(:,cc) = pft_coll_old(:,ddd)
         else
            pv_coll(cc) = .false.
            pfn_coll(:,cc) = 0.0
            pft_coll(:,cc) = 0.0
         endif

      enddo

! resetting do_nsearch to false here since neighbor search will have
! just been invoked
      DO_NSEARCH = .FALSE.

      RETURN
      END SUBROUTINE NEIGHBOUR
