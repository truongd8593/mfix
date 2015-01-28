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

!$omp parallel do default(none) private(cc) &
!$omp             shared(pair_num,pairs_old,pairs,pv_pair_old,pv_pair,pfn_pair_old,pfn_pair,pft_pair_old,pft_pair)
      do cc= 1, pair_num
         pairs_old(:,cc) = pairs(:,cc)
         pv_pair_old(cc) = pv_pair(cc)
         pfn_pair_old(:,cc) = pfn_pair(:,cc)
         pft_pair_old(:,cc) = pft_pair(:,cc)
      end do
!$omp end parallel do

      old_pair_num = pair_num
      pair_num = 0

      IF (DES_NEIGHBOR_SEARCH.EQ.1) THEN
         CALL NSQUARE
      ELSEIF (DES_NEIGHBOR_SEARCH.EQ.4) THEN
          CALL DESGRID_NEIGH_BUILD
      ENDIF

      dd = 1
      iii = pairs_old(1,dd)
      jjj = pairs_old(2,dd)

!$omp parallel do default(none) private(cc,ii,jj,iii,jjj,ddd) &
!$omp          shared(pair_num,pairs,pairs_old,pv_pair,pfn_pair,pft_pair,pfn_pair_old,pft_pair_old,pv_pair_old,old_pair_num) &
!$omp          firstprivate(dd)

      do cc = 1, pair_num
         ii = pairs(1,cc)
         jj = pairs(2,cc)

         iii = pairs_old(1,dd)
         jjj = pairs_old(2,dd)

         do while (dd .le. old_pair_num .and. (iii < ii))
            dd = dd + 1
            iii = pairs_old(1,dd)
            jjj = pairs_old(2,dd)
         enddo

         ddd = dd
         do while (ddd .le. old_pair_num .and. iii.eq.ii .and. jjj.ne.jj )
            ddd = ddd + 1
            iii = pairs_old(1,ddd)
            jjj = pairs_old(2,ddd)
         enddo

         if (ii.eq.iii .and. jj.eq.jjj) then
            pv_pair(cc) = pv_pair_old(ddd)
            pfn_pair(:,cc) = pfn_pair_old(:,ddd)
            pft_pair(:,cc) = pft_pair_old(:,ddd)
         else
            pv_pair(cc) = .false.
            pfn_pair(:,cc) = 0.0
            pft_pair(:,cc) = 0.0
         endif

      enddo
!$omp end parallel do

! resetting do_nsearch to false here since neighbor search will have
! just been invoked
      DO_NSEARCH = .FALSE.

      RETURN
      END SUBROUTINE NEIGHBOUR
