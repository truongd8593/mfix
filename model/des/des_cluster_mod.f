!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_CLUSTER                                            !
!                                                                      !
!  Purpose: Common elements for MFIX-DEM cluster identification        !
!  condition.                                                          !
!                                                                      !
!  Author: J.Galvin, J.Musser                         Date:  Nov-12    !
!  Modified: S. Benyahia                              Date:  Dec-12    !
!                                                                      !
!  Comments: Info on clusters such as average void fraction, Re and    !
!  Comments: cluster size can now be printed from this file            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE DES_CLUSTER


!-----------------------------------------------
! Modules      
!-----------------------------------------------
      USE param
      USE param1
      USE compar
      USE fldvar
      USE physprop
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------




! Number of clusters.
      INTEGER ClusterCount

! Number of base particles held in current search list
      INTEGER pSearchHistoryCount

!----------------------------------------------->>>
      TYPE PARTICLE_TYPE
! Global index of particle.
         INTEGER :: ID
! Successor in particle list
         TYPE (PARTICLE_TYPE), POINTER :: next_particle
      END TYPE PARTICLE_TYPE
!-----------------------------------------------<<<

! Basic particle link-list.      
      TYPE(PARTICLE_TYPE), POINTER :: PSEARCH_HISTORY_LL

!----------------------------------------------->>>      
      TYPE CLUSTER_TYPE
! Cluster index.
         INTEGER :: ID
! Number of particles in cluster.
         INTEGER :: ParticleCount
! Successor in cluster link-list.
         TYPE (CLUSTER_TYPE), POINTER :: next_cluster
! Particle link-list         
         TYPE (PARTICLE_TYPE), POINTER :: PARTICLE_LL
      END TYPE CLUSTER_TYPE
!-----------------------------------------------<<<      

! Cluster link-list.
      TYPE(CLUSTER_TYPE), POINTER :: CLUSTER_LL

      

      contains

!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE CREATE_CLUSTER(cluster)

      TYPE(CLUSTER_TYPE), INTENT(OUT), POINTER :: cluster

      ALLOCATE(cluster)

      NULLIFY(cluster%next_cluster)
      NULLIFY(cluster%PARTICLE_LL)

      if(ClusterCount == 0) then
! no clusters have been created/identified              
         if(associated(CLUSTER_LL)) then
            print*, ' Error - cluster pointer already associated!'
            CALL MFIX_EXIT(myPE)
         else
            ClusterCount = 1
            Cluster%ParticleCount = 0 
            cluster%ID = ClusterCount
! create first cluster of linked list of clusters. with cluster_ll
! always being the 'first' in the list
            CLUSTER_LL => cluster
         endif
      else
         if(.NOT.associated(CLUSTER_LL)) then
            print*, ' Error - cluster pointer is not associated!'
            CALL MFIX_EXIT(myPE)
         else
            ClusterCount = ClusterCount + 1
            Cluster%ParticleCount = 0 
            cluster%ID = ClusterCount
! establish the link between the new cluster and existing cluster list
            cluster%next_cluster => CLUSTER_LL
! reassign/point cluster_ll to be 'first' in the linked list
            CLUSTER_LL => cluster
! as a result the linked list of clusters is created from 'bottom-up' 
! with the new cluster always being inserted before any existing
! clusters 
         endif
      ENDIF

      END SUBROUTINE CREATE_CLUSTER


!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE GetTopCluster(cluster)

      TYPE(CLUSTER_TYPE), INTENT(INOUT), POINTER :: cluster
      
      NULLIFY(cluster)
      
! first check that a list has been formed      
      if(.NOT.associated(CLUSTER_LL) .AND. ClusterCount == 0) then
         write(*,"(//,3x,A,//)") 'No clusters to delete!'
      elseif(.NOT.associated(CLUSTER_LL) .AND. ClusterCount /= 0) then
         write(*,"(//,3x,A,//)") ' Error with ClusterCount and pointer - Er:1'
      elseif(associated(CLUSTER_LL) .AND. ClusterCount == 0) then
         write(*,"(//,3x,A,//)") ' Error with ClusterCount and pointer - Er:2'
      else
! now identify top cluster              
         CALL getNextCluster(cluster) 
      endif
      END SUBROUTINE GetTopCluster


!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE DeleteTopCluster(cluster)

      TYPE(CLUSTER_TYPE), INTENT(INOUT), POINTER :: cluster

! dbg      
!      write(*,"(/A,I7)")' Deleting top cluster: ',cluster%ID
      CALL DELETE_PARTICLES_IN_CLUSTER(cluster)

      if(associated(cluster%next_cluster)) then
         CLUSTER_LL => cluster%next_cluster
      else
         NULLIFY(CLUSTER_LL)
      endif

      ClusterCount = ClusterCount - 1
      NULLIFY(cluster%next_cluster)
      NULLIFY(cluster%PARTICLE_LL)
      DEALLOCATE(cluster)
      NULLIFY(cluster)

      if(ClusterCount <0) then
         write(*,"()")' ClusterCount < 0'
         CALL MFIX_EXIT(myPE)
      endif

      END SUBROUTINE DeleteTopCluster


!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE DELETE_CLUSTERS()

      TYPE(CLUSTER_TYPE), POINTER :: cluster
      INTEGER cL, cDeleted

      NULLIFY(cluster)
      cDeleted = 0

      if(.NOT.associated(CLUSTER_LL) .AND. ClusterCount == 0) then
         write(*,"(//,3x,A,//)") 'No clusters to delete!'
      elseif(.NOT.associated(CLUSTER_LL) .AND. ClusterCount /= 0) then
         write(*,"(//,3x,A,//)") ' Error with ClusterCount and pointer - Er:1'
      elseif(associated(CLUSTER_LL) .AND. ClusterCount == 0) then
         write(*,"(//,3x,A,//)") ' Error with ClusterCount and pointer - Er:2'
      else
         do cL =1, ClusterCount
            CALL getNextCluster(cluster)
! dbg            
!            write(*,"(/A,I7)")' Deleting cluster: ',cluster%ID
            CALL DELETE_PARTICLES_IN_CLUSTER(cluster)
 
            if(associated(cluster%next_cluster)) then
               CLUSTER_LL => cluster%next_cluster
            else
               NULLIFY(CLUSTER_LL)
            endif
            cDeleted = cDeleted + 1
            NULLIFY(cluster%next_cluster)
            NULLIFY(cluster%PARTICLE_LL)
            DEALLOCATE(cluster)
            NULLIFY(cluster)
         enddo

         if(cDeleted == ClusterCount) then
! dbg                 
!            write(*,"(4X,A,I7)")'Number of clusters deleted: ',cDeleted
            ClusterCount = 0
         else
            write(*,"()")' cDeleted /= ClusterCount'
            CALL MFIX_EXIT(myPE)
         endif

      endif

      END SUBROUTINE DELETE_CLUSTERS



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE DELETE_PARTICLES_IN_CLUSTER(cluster)

      TYPE(CLUSTER_TYPE), INTENT(IN), POINTER :: cluster

      TYPE(PARTICLE_TYPE), POINTER :: particle
      INTEGER pL, pDeleted

      NULLIFY(particle)
      pDeleted = 0

      if(.NOT.associated(cluster%PARTICLE_LL) .AND. cluster%ParticleCount == 0) then
         write(*,"(//,3x,A,//)") 'No particles to delete!'
      elseif(.NOT.associated(cluster%PARTICLE_LL) .AND. cluster%ParticleCount /= 0) then
         write(*,"(//,3x,A,//)") ' Error with ParticleCount and pointer - Er:1'
      elseif(associated(cluster%PARTICLE_LL) .AND. cluster%ParticleCount == 0) then
         write(*,"(//,3x,A,//)") ' Error with ParticleCount and pointer - Er:2'
      else
         do pL =1, cluster%ParticleCount
            CALL GetNextParticle(cluster, particle)
! dbg            
!            write(*,"(6X,A,I7)")' Deleting particle: ',particle%ID

            if(associated(particle%next_particle)) then
               cluster%PARTICLE_LL => particle%next_particle
            else
               NULLIFY(cluster%PARTICLE_LL)
            endif
            pDeleted = pDeleted + 1
            NULLIFY(particle%next_particle)
            DEALLOCATE(particle)
            NULLIFY(particle)
         enddo

         if(pDeleted == cluster%ParticleCount) then
! dbg                 
!            if(pDeleted >1 ) &
!               write(*,"(4X,A,I7)")'Number of deleted particles: ',&
!               pDeleted
            cluster%ParticleCount = 0
         else
            write(*,"()")' pDeleted /= cluster%ParticleCount'
            CALL MFIX_EXIT(myPE)
         endif

      endif


      END SUBROUTINE DELETE_PARTICLES_IN_CLUSTER



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE ADD_PARTICLE_TO_CLUSTER(cluster, pID)

      TYPE(CLUSTER_TYPE), INTENT(INOUT), POINTER :: cluster
      INTEGER, INTENT(IN) :: pID

      TYPE (PARTICLE_TYPE), POINTER :: particle

! Create a particle instance.
      ALLOCATE(particle)
      NULLIFY(particle%next_particle)
! Set the particle id.
      particle%ID = pID
! Set the cluster pointer to particle.
      IF(cluster%ParticleCount == 0) THEN
         cluster%ParticleCount = 1
         cluster%PARTICLE_LL => particle
      ELSE
         cluster%ParticleCount = cluster%ParticleCount + 1
         particle%next_particle => cluster%PARTICLE_LL
         cluster%PARTICLE_LL => particle
      ENDIF


      END SUBROUTINE ADD_PARTICLE_TO_CLUSTER

!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE PRINT_CLUSTERS

      INTEGER cL, pL
      INTEGER I, J, K, IJK, L, M, xLmin, yLmin, zLmin, xLmax, yLmax, zLmax
      DOUBLE PRECISION avg_epg, avg_ug, avg_vg, avg_wg, &
                       avg_us, avg_vs, avg_ws, avg_Re, avg_Slip
      DOUBLE PRECISION xpos_min, ypos_min, zpos_min, xpos_max, ypos_max, zpos_max
      DOUBLE PRECISION clXsize, clYsize, clZsize, clDiameter
      LOGICAL countThisCluster

      TYPE(CLUSTER_TYPE), POINTER :: cluster
      TYPE(PARTICLE_TYPE), POINTER :: particle

      NULLIFY(cluster)

      PostCluster(:) = ZERO  ! initialize to zero
      open (unit=201,file='clusterInfo.dat',status='unknown',position='append')
      if(.NOT.associated(CLUSTER_LL)) THEN
         if(mype == pe_io) write(201,"(/2X,A)")'No clusters to print!'
      else
         if(mype == pe_io) write(201,"(/)")  ! skip one line to identify new entry
         
         do cL = 1, ClusterCount
            CALL getNextCluster(cluster)
            NULLIFY(particle)
! only clusters with more than 3 particles are processed
            if(cluster%ParticleCount > 3) then
               avg_epg = zero
	       avg_ug = zero
	       avg_vg = zero
	       avg_wg = zero
	       avg_us = zero
	       avg_vs = zero
	       avg_ws = zero
	       avg_Re = zero
	       xpos_min = large_number  ! initiating cluster minimum x-pos
	       ypos_min = large_number
	       zpos_min = large_number
	       xpos_max = zero          ! initiating cluster maximum x-pos
	       ypos_max = zero
	       zpos_max = zero
! need to know particle index to add particle radius to cluster	
	       xLmin = 0
	       yLmin = 0
	       zLmin = 0
	       xLmax = 0
	       yLmax = 0
	       zLmax = 0
	       countThisCluster = .true.
	       LOOP_PL : do pL = 1, cluster%ParticleCount
                  CALL GetNextParticle(cluster, particle)
		  L = particle%ID
                  I = PIJK(L,1)
                  J = PIJK(L,2)
                  K = PIJK(L,3)
                  IJK = PIJK(L,4)
		  if(IJK == 0) then
		    countThisCluster = .false.
		    EXIT LOOP_PL ! do not account for clusters across processors
		  endif
		  PostCluster(L) = float(cluster%ParticleCount)
                  M = PIJK(L,5)
		  avg_epg = avg_epg + ep_g(ijk)
		  avg_ug = avg_ug + u_g(ijk)
		  avg_vg = avg_vg + v_g(ijk)
		  avg_wg = avg_wg + w_g(ijk)
		  avg_us = avg_us + DES_VEL_NEW(L,1)
		  avg_vs = avg_vs + DES_VEL_NEW(L,2)
		  avg_ws = avg_ws + DES_VEL_NEW(L,3)

! calculating an average Re like the one below or by using average quantities
! doesn't seem to make much difference.
		  avg_Re = avg_Re + ep_g(ijk) * dsqrt((u_g(ijk)-DES_VEL_NEW(L,1))**2 + &
		  (v_g(ijk)-DES_VEL_NEW(L,2))**2 +(w_g(ijk)-DES_VEL_NEW(L,3))**2) * &
		  (2d0*DES_RADIUS(L)) / MU_g(ijk)
!
! now determining the size of a cluster	  
		  if(DES_POS_NEW(L,1) < xpos_min) then
		     xpos_min = DES_POS_NEW(L,1)
		     xLmin = L
		  endif
		  if(DES_POS_NEW(L,2) < ypos_min) then
		     ypos_min = DES_POS_NEW(L,2)
		     yLmin = L
		  endif
		  if(DES_POS_NEW(L,3) < zpos_min) then
		     zpos_min = DES_POS_NEW(L,3)
		     zLmin = L
		  endif
		  
		  if(DES_POS_NEW(L,1) > xpos_max) then
		     xpos_max = DES_POS_NEW(L,1)
		     xLmax = L
		  endif
		  if(DES_POS_NEW(L,2) > ypos_max) then
		     ypos_max = DES_POS_NEW(L,2)
		     yLmax = L
		  endif
		  if(DES_POS_NEW(L,3) > zpos_max) then
		     zpos_max = DES_POS_NEW(L,3)
		     zLmax = L
		  endif

               enddo LOOP_PL
	       if(countThisCluster) then
	         avg_epg = avg_epg / cluster%ParticleCount
	         avg_ug = avg_ug / cluster%ParticleCount
	         avg_vg = avg_vg / cluster%ParticleCount
	         avg_wg = avg_wg / cluster%ParticleCount
	         avg_us = avg_us / cluster%ParticleCount
	         avg_vs = avg_vs / cluster%ParticleCount
	         avg_ws = avg_ws / cluster%ParticleCount
	         avg_Re = avg_Re / cluster%ParticleCount
	         avg_Slip = dsqrt((avg_ug-avg_us)**2 + (avg_vg-avg_vs)**2 +(avg_wg-avg_ws)**2)
	         clXsize = xpos_max - xpos_min + DES_RADIUS(xLmin) + DES_RADIUS(xLmax)
	         clYsize = ypos_max - ypos_min + DES_RADIUS(yLmin) + DES_RADIUS(yLmax)
	         clZsize = zpos_max - zpos_min + DES_RADIUS(zLmin) + DES_RADIUS(zLmax)

! Cluster diameter is calculated by using an equivalent circle diameter
! of an ellipse in each plane. Then weighting this area by the slip
! velocity normal to it. This in effect gives a hydrodynamic diameter
! of the cluster.
	         clDiameter = dsqrt( clYsize*clZsize*(avg_ug-avg_us)**2 + & 
	                             clXsize*clZsize*(avg_vg-avg_vs)**2 + &
                                     clXsize*clYsize*(avg_wg-avg_ws)**2 ) / &
			             avg_Slip

! the scaling below was commented so to use an average diameter for poly cases
!	         clDiameter = clDiameter/(2d0*DES_RADIUS(xLmax))

! Below we are writing average properties associated with a cluster:
! cluster ID, particles count in cluster, viod fraction, Re number, 
! slip velocity, and cluster diameter. Anything else can be calculated
! and written if wanted.
                write(201,"(4X,I7,4X,I7, 4(1X,G13.6))") cluster%ID, cluster%ParticleCount, &
	           avg_epg, avg_Re, avg_Slip, clDiameter
	       endif
            endif
         enddo
! point cluster to first cluster in linked list of clusters
      endif
      close (unit=201)

      END SUBROUTINE PRINT_CLUSTERS



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE getNextCluster(cluster)

      TYPE(CLUSTER_TYPE), INTENT(INOUT), POINTER :: cluster

      if(.NOT.associated(cluster)) then
! point cluster to first cluster in linked list of clusters              
         cluster => CLUSTER_LL
      elseif(associated(cluster%next_cluster)) then
         cluster => cluster%next_cluster
      else
         print*,' You are looking for a cluster that does not exist'
         CALL MFIX_EXIT(myPE)
      endif
      END SUBROUTINE getNextCluster



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE GetNextParticle(cluster, particle)

      TYPE(CLUSTER_TYPE), INTENT(IN), POINTER :: cluster
      TYPE(PARTICLE_TYPE), INTENT(INOUT), POINTER :: particle

      if(.NOT.associated(particle)) then
         particle => cluster%PARTICLE_LL
      elseif(associated(particle%next_particle)) then
         particle => particle%next_particle
      else
         print*,' You are looking for a particle that does not exist'
         CALL MFIX_EXIT(myPE)
      endif
      END SUBROUTINE GetNextParticle



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE ADD_PARTICLE_TO_PSEARCHHISTORY(particle, pID)

      TYPE(PARTICLE_TYPE), INTENT(OUT), POINTER :: particle
      INTEGER, INTENT(IN) :: pID

      ALLOCATE(particle)

      NULLIFY(particle%next_particle)

!dbg
!      print *, '----------> addparticle_to_psearchhistory'

      if(pSearchHistoryCount == 0) then
! no particles in search list have been created/identified              
         if(associated(PSEARCH_HISTORY_LL)) then
                 print*, ' Error - particle history pointer already &
                    & associated!'
            CALL MFIX_EXIT(myPE)
         else
            pSearchHistoryCount = 1
            particle%ID = pID
! create first particle of linked list of particles. with
! psearch_history_ll always pointing to the 'first' in the list
            PSEARCH_HISTORY_LL => particle
! dbg
!            write(*,"(/A,I7)") 'Create history w/ particle ', pID
         endif
      else
         if(.NOT.associated(PSEARCH_HISTORY_LL)) then
            print*, ' Error - particle history  pointer is not &
               & associated!'
            CALL MFIX_EXIT(myPE)
         else
            pSearchHistoryCount = pSearchHistoryCount + 1
            particle%ID = pID
! establish the link between the new particle and existing particle list
            particle%next_particle => PSEARCH_HISTORY_LL
! reassign/point psearch_history_ll to be 'first' in the linked list
            PSEARCH_HISTORY_LL => particle
! as a result the linked list of particles is created from 'bottom-up' 
! with the new particle  always being inserted before any existing
! particles 
! dbg
!            print *, 'Add particle in history ', pID
            
         endif
      ENDIF

      END SUBROUTINE ADD_PARTICLE_TO_PSEARCHHISTORY



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE GetTopParticle_in_PSearchHistory(particle)

      TYPE(PARTICLE_TYPE), INTENT(INOUT), POINTER :: particle

!dbg
!     print *, '----------> gettopparticle_in_psearchhistory'

      if(.NOT.associated(particle)) then
! point particle to first particle in linked list of particles              
         particle => PSEARCH_HISTORY_LL
      elseif(associated(particle%next_particle)) then
! this will never happen if we nullify particle at start of this
! routine.
         particle => particle%next_particle
      else
         print*,' You are looking for a particle that does not exist'
         CALL MFIX_EXIT(myPE)
      endif
      END SUBROUTINE GetTopParticle_in_PSearchHistory



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE DeleteTopParticle_in_PSearchHistory()

      TYPE(PARTICLE_TYPE), POINTER :: particle

      NULLIFY(particle)
!dbg
!      print *, '----------> deletetopparticle_in_psearchhistory'

      if(.NOT.associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount == 0) then
         write(*,"(//,3x,A,//)") 'No particles in history to delete!'
      elseif(.NOT.associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount /= 0) then
         write(*,"(//,3x,A,//)") &
            ' Error with pSearchHistoryCount and pointer - Er:1'
      elseif(associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount == 0) then
         write(*,"(//,3x,A,//)") &
            ' Error with pSearchHistoryCount and pointer - Er:2'
      else
         CALL GetTopParticle_In_PSearchHistory(particle)
! dbg         
!         write(*,"(A,I7)")' Deleting particle: ',particle%ID
   
         if(associated(particle%next_particle)) then
            PSEARCH_HISTORY_LL => particle%next_particle
         else
            NULLIFY(PSEARCH_HISTORY_LL)
         endif
         pSearchHistoryCount = pSearchHistoryCount - 1
         NULLIFY(particle%next_particle)
         DEALLOCATE(particle)
         NULLIFY(particle)

         if(pSearchHistoryCount < 0) then
            write(*,"(4X,A)")'pSearchHistoryCount < 0'
            CALL MFIX_EXIT(myPE)
         endif

      endif

      END SUBROUTINE DeleteTopParticle_In_PSearchHistory



!......................................................................!
!                                                                      !
!......................................................................!
      SUBROUTINE DELETE_PSEARCHHISTORY()

      TYPE(PARTICLE_TYPE), POINTER :: particle
      INTEGER pL, pDeleted

      NULLIFY(particle)
      pDeleted = 0
!dbg
!      print *, '----------> delete_psearch_history'

      if(.NOT.associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount == 0) then
! dbg         
!         write(*,"(//,3x,A,//)") 'No particles in history to delete!'
      elseif(.NOT.associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount /= 0) then
         write(*,"(//,3x,A,//)") &
            ' Error with pSearchHistoryCount and pointer - Er:1'
      elseif(associated(PSEARCH_HISTORY_LL) .AND. &
         pSearchHistoryCount == 0) then
         write(*,"(//,3x,A,//)") &
            ' Error with pSearchHistoryCount and pointer - Er:2'
      else
         do pL =1, pSearchHistoryCount
        
            CALL GetTopParticle_In_PSearchHistory(particle)
! dbg            
!            write(*,"(/A,I7)")' Deleting particle: ',particle%ID
    
            if(associated(particle%next_particle)) then
               PSEARCH_HISTORY_LL => particle%next_particle
            else
               NULLIFY(PSEARCH_HISTORY_LL)
            endif
            pDeleted = pDeleted + 1
            NULLIFY(particle%next_particle)
            DEALLOCATE(particle)
            NULLIFY(particle)
         enddo

         if(pDeleted == pSearchHistoryCount) then
! dbg                 
!            write(*,"(4X,A,I7)")'Number of particles deleted: ',pDeleted
            pSearchHistoryCount = 0
         else
            write(*,"()")' pDeleted /= PSearchHistoryCount'
            CALL MFIX_EXIT(myPE)
         endif

      endif

      END SUBROUTINE DELETE_PSEARCHHISTORY


      
      END MODULE DES_CLUSTER
