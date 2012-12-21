!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BASE_PARTICLE_NEIGHBORS                           !
!                                                                      !
!  Purpose: This routine searches through the base particle neighbor   !
!  list to identify any neighbor that constitutes a cluster with the   !
!  base particle. If a neighbor is identified as forming a cluster     !
!  then the routine is repeated with that neighbor particle as the new !
!  base particle. A list is maintained of each particle identified and !
!  used as a base particle. This process is repeated recursively until !
!  a particle is reached that has no neighbors constituting a cluster. !
!  At this point the code sequentially backtracks through the list of  !
!  base particles and repeats the search routine, skipping over any    !
!  particle already identified as belonging to the cluster.            !
!  This traverses through all possible neighbor branches to identify   !
!  clusters with the base particle.                                    !
!                                                                      !
!                                                                      !
!  Author: J.Galvin                                   Date: Nov 2012   !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      RECURSIVE SUBROUTINE CHECK_BASE_PARTICLE_NEIGHBORS(BaseCluster,&
         BaseParticle)

!-----------------------------------------------
! Modules      
!-----------------------------------------------
      Use discretelement
      Use run
      Use des_cluster
      IMPLICIT NONE      
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! id of primary/base particle      
      INTEGER, INTENT(INOUT) :: BaseParticle
! current cluster formed from particles in search      
      TYPE(CLUSTER_TYPE), INTENT(INOUT), POINTER :: BaseCluster      
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop index
      INTEGER :: L
! particle index
      INTEGER :: SearchParticle
! sum of particle radii
      DOUBLE PRECISION :: R_LM
! distance vector between two particle centers 
      DOUBLE PRECISION :: DIST(DIMN)
! magnitude of distance between two particle centers
      DOUBLE PRECISION :: DISTMOD
! distance between particle surfaces
      DOUBLE PRECISION :: DistApart

      TYPE(PARTICLE_TYPE), POINTER :: particle  
!-----------------------------------------------
! Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------

!-----------------------------------------------

! checking if base particle has neighbors
      IF(NEIGHBOURS(BaseParticle,1)>0) THEN
! looping over base particle neighbors..
         DO L = 2, NEIGHBOURS(BaseParticle,1)+1
! defining particle index of neighboring particle
            SearchParticle = NEIGHBOURS(BaseParticle,L)
! checking only those particles that ‘exist’
            IF(PEA(SearchParticle,1)) THEN
! skipping particles that are in a cluster
               IF(InACluster(SearchParticle)) CYCLE  
! calculating distance between surface of search particle and base
! particle
               R_LM = des_radius(SearchParticle)+des_radius(BaseParticle)
               Dist(:) = des_pos_new(SearchParticle,:) - des_pos_new(BaseParticle,:)
               Distmod = sqrt(des_dotprdct(dist,dist))
               DistApart= (distmod-R_LM)
! checking if search particle and base particle form a cluster based on
! defined criterion            
               IF(DistApart<Cluster_Length_Cutoff) THEN 
! adding search particle to cluster of base particle                    
                  CALL ADD_PARTICLE_TO_CLUSTER(BaseCluster,SearchParticle)
! flag neighbor particle as belonging to a cluster               
                  InACluster(SearchParticle) = .TRUE.
! adding search particle to a list that tracks each subsequent base
! particle
                  CALL ADD_PARTICLE_TO_PSEARCHHISTORY(particle,SearchParticle)
! repeating routine with new BaseParticle=SearchParticle
                  CALL CHECK_BASE_PARTICLE_NEIGHBORS(BaseCluster,&
                     SearchParticle)
               ENDIF
            ENDIF  ! endif pea(searchparticle,1)
! continue loop with current base particle until we reach a search 
! particle that is outside cutoff distance or no more neighbors         
         ENDDO   ! end loop over neighbours for current base particle..
  

         IF(pSearchHistoryCount > 0) THEN
! removing top particle from list. at this point the current base
! particle should be the top of the search history list. so it does not
! need to be re-searched as all its neighbors have been checked. remove
! it before going to the next particle in list, if one exists.
            CALL DELETETOPPARTICLE_IN_PSEARCHHISTORY()
            IF (pSearchHistoryCount > 0) THEN
! getting top particle in psearch_history_list
               NULLIFY(particle)
               CALL GETTOPPARTICLE_IN_PSEARCHHISTORY(particle)
               SearchParticle = particle%ID
! dbg               
!               write(*,"(A,I7,/)")' Recovering search particle: ', SearchParticle
! repeating routine with new baseparticle
               CALL CHECK_BASE_PARTICLE_NEIGHBORS(BaseCluster,SearchParticle)
            ENDIF
         ENDIF
      ENDIF

! if no more particles in list, exit routine. done with base particle
! branch (i.e., searching its neighbors, and their neighbors, etc.)
     
      RETURN
      END SUBROUTINE CHECK_BASE_PARTICLE_NEIGHBORS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: IDENTIFY_SYSTEM_CLUSTERS                                !
!                                                                      !
!  Purpose: This routine will examine a snapshot of the system to      !
!           identify particle clusters based on user defined critera   !     
!                                                                      !
!  Author: J.Galvin                                   Date: Nov 2012   !      
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE IDENTIFY_SYSTEM_CLUSTERS

!-----------------------------------------------
! Modules      
!-----------------------------------------------
      Use discretelement
      Use des_cluster
      Use run
      IMPLICIT NONE      
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! particle index
      INTEGER :: L

      TYPE(PARTICLE_TYPE), POINTER :: particle

      TYPE(CLUSTER_TYPE), POINTER :: cluster

      INTERFACE
         SUBROUTINE CHECK_BASE_PARTICLE_NEIGHBORS(cluster,L)
            Use des_cluster
            TYPE(CLUSTER_TYPE), POINTER :: cluster
            INTEGER :: L
         END SUBROUTINE CHECK_BASE_PARTICLE_NEIGHBORS
      END INTERFACE

!-----------------------------------------------
! Include statement functions
!-----------------------------------------------

!-----------------------------------------------

! dbg      
      if(mype == pe_io) write(*,'(/3X,A,F8.2,/)') &
         'IDENTIFY_SYSTEM_CLUSTERS CALLED at time', time

! initializing      
      ClusterCount = 0
      PSearchHistoryCount = 0
      InACluster(:) = .FALSE.

      DO L = 1, MAX_PIP
         IF (.NOT.PEA(L,1)) CYCLE  ! check if particle exists (cycle if not)
         IF (PEA(L,4)) CYCLE       ! check if ghost particle (cycle if ghost)
         IF(.NOT.InACluster(L)) THEN    ! checking if particle already belongs to a cluster
! creating a cluster
            CALL CREATE_CLUSTER(cluster)
! adding particle L to cluster
            CALL ADD_PARTICLE_TO_CLUSTER(cluster,L)
! flagging particle L as in a cluster
            InACluster(L) = .TRUE.
! creating particle list that tracks any subsequent neighbor particle
! searches using a new base particle
            CALL ADD_PARTICLE_TO_PSEARCHHISTORY(particle,L)
! start neighbor check
            CALL CHECK_BASE_PARTICLE_NEIGHBORS(cluster,L)
         ENDIF
! If a particle already exists in a cluster then all possible cluster
! candidates of that particle should have already been checked. That is, 
! the branch of its neighbors, their neighbors, etc., should have been
! examined. So we can move to next base particle

! create routine where if cluster is only 1 particle size big delete
! it.
         CALL GetTopCluster(cluster)
         IF (cluster%ParticleCount == 1) THEN
            CALL DeleteTopCluster(cluster)
! I don't think this routine call is necessary but it is a fail-safe..            
            CALL DELETE_PSEARCHHISTORY            
         ENDIF

      ENDDO
      
      CALL PRINT_CLUSTERS()
      CALL DELETE_CLUSTERS()

! dbg      
      if(mype == pe_io) write(*,'(3X,A)') 'END IDENTIFY_SYSTEM_CLUSTERS CALLED'
      RETURN
      END SUBROUTINE IDENTIFY_SYSTEM_CLUSTERS
