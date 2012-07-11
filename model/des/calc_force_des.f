!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_FORCE_DES                                          !
!  Purpose: Calculate contact force and torque on particle from        !
!           particle-particle and particle-wall collisions. Treats     !
!           wall interaction also as a two-particle interaction but    !
!           accounting for the wall properties                         !
!                                                                      !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!  Reviewer: Sreekanth Pannala                        Date: 06-Dec-06  !
!  Reviewer: Rahul Garg                               Date: 02-Aug-07  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CALC_FORCE_DES

!-----------------------------------------------
! Modules
!----------------------------------------------- 
      USE run      
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------      
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0
! particle no. indices      
      INTEGER :: I, LL
! loop counters (neighbor particles, walls, dimension)
      INTEGER :: II, IW, IDIMN
! loop counter and indices for neighbor information      
      INTEGER :: NI, NLIM, N_NOCON, NEIGH_L
! the overlap occuring between particle-particle or particle-wall
! collision in the normal and tangential directions
      DOUBLE PRECISION :: OVERLAP_N, OVERLAP_T
! square root of the overlap
      DOUBLE PRECISION :: SQRT_OVERLAP
      DOUBLE PRECISION :: FRAC_OVERLAP1, FRAC_OVERLAP2
! distance vector between two particle centers or between a particle
! center and wall when the two surfaces are just at contact (i.e. no
! overlap)
      DOUBLE PRECISION :: R_LM
! the normal and tangential components of the translational relative
! velocity
      DOUBLE PRECISION :: V_REL_TRANS_NORM, V_REL_TRANS_TANG, &
                          V_REL_NORM_OLD, VRN_OLD(DIMN)
! time elapsed to travel the calculated normal overlap given the 
! normal relative velocity
      DOUBLE PRECISION :: DTSOLID_TMP
! distance vector between two particle centers or between a particle
! center and wall at current and previous time steps
      DOUBLE PRECISION :: DIST(DIMN), DIST_OLD(DIMN)
! magnitude of distance between two particle centers or between a
! particle center and wall at current and previous time steps
      DOUBLE PRECISION :: DISTMOD, DISTMOD_OLD
! unit normal vector along the line of contact between contacting
! particles or particle-wall at current and previous time steps
      DOUBLE PRECISION :: NORMAL(DIMN), NORM_OLD(DIMN), NORMAL_OLD(DIMN)
! tangent to the plane of contact at current time step
      DOUBLE PRECISION :: TANGENT(DIMN)
! variables for tangential displacement calculation:
! unit vector for axis of rotation and its magnitude
      DOUBLE PRECISION :: TMP_AX(DIMN), TMP_MAG
! local variables for the accumulated tangential displacement that occurs
! for the particle-particle or particle-wall collision (current time
! step and previous time step)
      DOUBLE PRECISION :: SIGMAT(DIMN),SIGMAT_OLD(DIMN)
! tangent to the plane of contact at current and previous time step
! (used for 2D calculations)      
      DOUBLE PRECISION :: TANG_OLD(DIMN),TANG_NEW(DIMN)
! normal and tangential forces      
      DOUBLE PRECISION :: FNS1(DIMN), FNS2(DIMN)
      DOUBLE PRECISION :: FTS1(DIMN), FTS2(DIMN)
! temporary storage of tangential FORCE
      DOUBLE PRECISION :: FT_TMP(DIMN)
! temporary storage of tangential DISPLACEMENT
      DOUBLE PRECISION :: PFT_TMP(DIMN)
! magnitude of normal/tangential force used for reporting only
      DOUBLE PRECISION :: FTMD, FNMD

! store solids phase index of particle (i.e. pijk(np,5))
      INTEGER :: PHASEI, PHASELL
! local values used spring constants and damping coefficients
      DOUBLE PRECISION :: ETAN_DES, ETAN_DES_W, ETAT_DES, ETAT_DES_W,&
                          KN_DES, KN_DES_W, KT_DES, KT_DES_W
! local values used for calculating cohesive forces
      DOUBLE PRECISION :: FORCE_COH, EQ_RADIUS, DistApart, &
                          Norm_Dist, magGravity 
! flag indicating when a particle is in contact with a given wall
! (1=contact, 0=nocontact)
      INTEGER :: WALLCONTACT
! local variables to store wall position and velocity 
      DOUBLE PRECISION :: WALL_POS(DIMN), WALL_VEL(DIMN)
! set to T when a sliding contact occurs
      LOGICAL :: PARTICLE_SLIDE  
! logic flag telling whether contact pair is old (previously detected) 
      LOGICAL :: ALREADY_NEIGHBOURS
    
! logic flag for local debug warnings
      LOGICAL :: DES_LOC_DEBUG
      LOGICAL :: ALREADY_EXISTS

!$      double precision omp_start, omp_end
!$      double precision omp_get_wtime     
    
!-----------------------------------------------      
! Functions
!-----------------------------------------------      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
!-----------------------------------------------      


! initialize local variables
      DES_LOC_DEBUG = .FALSE.
      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1


! initialize cohesive forces
      IF(USE_COHESION) THEN
         Fcohesive(:,:) = ZERO
         PostCohesive (:) = ZERO
      ENDIF


!$omp   parallel default(shared)                                  & 
!$omp   private(ll,fts1,fts2,fns1,fns2,ft_tmp,pft_tmp,            &
!$omp          PARTICLE_SLIDE,nlim,                               &
!$omp          n_nocon,ni,iw,wallcontact,i,             &
!$omp          already_neighbours,neigh_l,                        &
!$omp          WALL_POS,WALL_VEL,                                 &
!$omp          r_lm,dist,distmod,frac_overlap1,normal,            &
!$omp          v_rel_trans_norm,tangent,                          &
!$omp          v_rel_trans_tang,                                  &
!$omp          overlap_n,overlap_t,dtsolid_tmp,phasell,           &
!$omp          sqrt_overlap,kn_des_w,kt_des_w,etan_des_w,         &
!$omp          etat_des_w,sigmat_old,norm_old,tmp_ax,tmp_mag,     &
!$omp          tang_old,tang_new,sigmat,                          &
!$omp          ftmd,fnmd,                                         &
!$omp          ii,frac_overlap2,                                  &
!$omp          dist_old,distmod_old,normal_old,                   &
!$omp          vrn_old,v_rel_norm_old,                            &
!$omp          phasei,kn_des,kt_des,etan_des,etat_des,            &
!$omp          Idimn,force_coh,eq_radius,distapart,norm_dist,maggravity)          
!$omp do reduction(max:NEIGH_MAX,OVERLAP_MAX) schedule (dynamic,50)


! Looping over all particles in the processor      
!----------------------------------------------------------------->>>
      DO LL = 1, MAX_PIP
     
! skipping non-existent particles      
         IF(.NOT.PEA(LL,1)) CYCLE
! skipping ghost particles         
         IF(PEA(LL,4)) CYCLE
         
! debugging purposes         
         IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
            IF (.NOT.DES_LOC_DEBUG) THEN
               DES_LOC_DEBUG = .TRUE.
               WRITE(*,1000)
            ENDIF
            WRITE(*,'(5X,A,I10)') 'On Particle ', LL
            WRITE(*,'(5X,A,2(ES15.7))') &
               'X,Y POS: ', DES_POS_NEW(LL,1), DES_POS_NEW(LL,2)
            WRITE(*,'(5X,A,2(ES15.7))') &
               'X,Y VEL: ', DES_VEL_NEW(LL,1), DES_VEL_NEW(LL,2)
         ENDIF

! Initializing local variables
         TANGENT(:) = ZERO
         NORMAL(:) = ZERO
         FTS1(:) = ZERO
         FTS2(:) = ZERO
         FNS1(:) = ZERO
         FNS2(:) = ZERO
         PFT_TMP(:) = ZERO
         FT_TMP(:) = ZERO
         PARTICLE_SLIDE = .FALSE.


! Check neighbor history of particle LL and update arrays as needed
! ---------------------------------------------------------------->>>         
         IF(PN(LL,1).GE.1) THEN
            NLIM = PN(LL,1)+1
            N_NOCON = 0
! For each particle listed as in contact with particle LL in array PN,
! check the flag array PV to determine if particles remained in contact 
! after the previous call of CALC_FORCE_DES. 
            DO NI = 2, NLIM
               IF(PV(LL,NI-N_NOCON).EQ.0) THEN
! For each particle in PN(LL,2:MAXNEIGHBORS) that is no longer in 
! contact, shift the remaining particle contact information PN, PN,
! PFT left by one and reduce PN(LL,1) by one.
                  PN(LL,(NI-N_NOCON):MAXNEIGHBORS-1) = &
                     PN(LL,(NI-N_NOCON+1):MAXNEIGHBORS) 
                  PV(LL,(NI-N_NOCON):(MAXNEIGHBORS-1)) = &
                     PV(LL,(NI-N_NOCON+1):MAXNEIGHBORS)
                  PFT(LL,(NI-N_NOCON):(MAXNEIGHBORS-1),:) = &
                     PFT(LL,(NI-N_NOCON+1):MAXNEIGHBORS,:)
! Save the normal direction at previous time step
                  PFN(LL,(NI-N_NOCON):(MAXNEIGHBORS-1),:) = &
                     PFN(LL,(NI-N_NOCON+1):MAXNEIGHBORS,:)
                  N_NOCON = N_NOCON + 1
                  PN(LL,1) = PN(LL,1) - 1
               ENDIF
            ENDDO
         ENDIF

! Initializing rest of the neighbor list which is not in contact and
! clean up after the above array left shifts
         NLIM = MAX(2,PN(LL,1) + 2) 
         PN(LL,NLIM:MAXNEIGHBORS) = -1
         PFT(LL,NLIM:MAXNEIGHBORS,:) = ZERO
         PFN(LL,NLIM:MAXNEIGHBORS,:) = ZERO   
         IF (PN(LL,1) .GT. NEIGH_MAX) NEIGH_MAX = PN(LL,1)


! Initializing the neighbor list contact information when particles are
! not in contact; i.e. when particle LL has no neighbors         
         IF (PN(LL,1).EQ.0) THEN
            PFT(LL,:,:) = ZERO
            PFN(LL,:,:) = ZERO
         ENDIF
         
! Reset the flag array PV; during each call to calc_force_des this
! variable tracks whether particle LL has any current neighbors
! the array is used in the next call to calc_force_des to update
! particle LL neighbor history above
         PV(LL,2:MAXNEIGHBORS) = 0

! End check neighbor history and update of corresponding history arrays         
! ----------------------------------------------------------------<<<

 

! Check particle LL for wall contacts
!----------------------------------------------------------------->>>
! Treats wall interaction also as a two-particle interaction but accounting
! for the wall properties; make sure the particle is not classified as
! a new 'entering' particle or is already marked as a potential exiting particle
         IF(WALLDTSPLIT .AND. .NOT.PEA(LL,2) .AND. .NOT.PEA(LL,3)) THEN

            DO IW = 1, NWALLS
               WALLCONTACT = 0

! Check to see if a particle is in contact with any of the walls
               CALL CFWALLCONTACT(IW, LL, WALLCONTACT)
              
               IF(WALLCONTACT.EQ.1) THEN
                  I = MAX_PIP + IW
                  ALREADY_NEIGHBOURS=.FALSE.
                  
                  IF(PN(LL,1).GT.0) THEN                     
                     DO NEIGH_L = 2, PN(LL,1)+1
                        IF(I.EQ. PN(LL,NEIGH_L)) THEN 
                           ALREADY_NEIGHBOURS=.TRUE.
                           NI = NEIGH_L
                           EXIT
                        ENDIF                        
                     ENDDO
                  ENDIF
                  
! Assign the wall particle a position and velocity
                  CALL CFWALLPOSVEL(LL, IW, WALL_POS, WALL_VEL)

                  R_LM = DES_RADIUS(LL) + DES_RADIUS(LL) 
                  DIST(:) = WALL_POS(:) - DES_POS_NEW(LL,:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))
                  
! compute particle-wall VDW cohesive short-range forces
                  IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
                    DistApart = (DISTMOD-R_LM) ! distance between particle&wall surface
                    IF(DistApart < WALL_VDW_OUTER_CUTOFF)THEN
                      IF(DistApart > WALL_VDW_INNER_CUTOFF)THEN
                         FORCE_COH = WALL_HAMAKER_CONSTANT * DES_RADIUS(LL) / &
                             (6d0*DistApart**2) * &
                             ( Asperities/(Asperities+DES_RADIUS(LL)) + &
                             ONE/(ONE+Asperities/DistApart)**2 )
                      ELSE

                         FORCE_COH = 4d0 * PI * WALL_SURFACE_ENERGY * DES_RADIUS(LL) * &
                              ( Asperities/(Asperities+DES_RADIUS(LL)) + &
                              ONE/(ONE+Asperities/VDW_INNER_CUTOFF)**2 )
                      ENDIF
                      DO IDIMN=1,DIMN
                        Norm_Dist = DIST(IDIMN)/DISTMOD
                        Fcohesive(LL, IDIMN) = Fcohesive(LL, IDIMN) + Norm_Dist*FORCE_COH
                      ENDDO 
                    ENDIF     
                  ENDIF ! for using VDW cohesion model
                  


                  IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN 
! Particle-wall overlap detected (i.e. resolve collision)
!---------------------------------- 

! Error reporting 
                     FRAC_OVERLAP1 = (R_LM-DISTMOD)/DES_RADIUS(LL) 
                     IF (FRAC_OVERLAP1 > flag_overlap) THEN
                        WRITE(*,'(5X,A,A,ES15.7)') &
                           'WARNING: excessive overlap detected ', &
                           'at time ', S_TIME
                        WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
                           'between particle ', LL, 'and wall ',&
                           IW, 'with'
                        WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7)') &
                          'overlap = ', (R_LM-DISTMOD), &
                           ' radius = ', DES_RADIUS(LL)
                     ENDIF

                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF                             
                        WRITE(*,'(5X,A,I10,A,I5)') &
                           'DISTMOD is zero between particle ', LL, &
                           ' and wall ', IW
                        STOP
                     ENDIF

! Calculate the components of translational relative velocity for a 
! contacting particle wall and the tangent to the plane of contact
                     CALL CFRELVEL_WALL(LL, WALL_VEL, V_REL_TRANS_NORM, &
                        V_REL_TRANS_TANG, TANGENT, NORMAL, DISTMOD)
                                         
! The normal overlap calculation was changed so that it no longer
! depends on the contact history (i.e., integration of incremental
! overlap found by velocity*dtsolid).  Now overlap is based purely on
! position of neighbors.  The change was made because the former
! method was found not to conserve energy 
                     OVERLAP_N =  R_LM-DISTMOD 

                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE 
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        IF (V_REL_TRANS_NORM .GT. ZERO) THEN
                          DTSOLID_TMP = OVERLAP_N/(V_REL_TRANS_NORM)
                        ELSEIF (V_REL_TRANS_NORM .LT. ZERO) THEN
                          DTSOLID_TMP = DTSOLID
                          WRITE(*,'(5X,A,A,ES15.7)') &
                             'WARNING: normal relative velocity less ',&
                             'than zero at time ', S_TIME
                          WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
                             'for first contact between particle', LL, &
                             'and wall ', IW, 'with'
                          WRITE(*,'(7X,A,ES15.7)') &
                             'V_REL_NORM = ', V_REL_TRANS_NORM
                        ELSE
                           DTSOLID_TMP = OVERLAP_N/&
                              (V_REL_TRANS_NORM+SMALL_NUMBER)
                        ENDIF
                        OVERLAP_T = V_REL_TRANS_TANG*&
                           MIN(DTSOLID,DTSOLID_TMP)
                     ENDIF


                     phaseLL = PIJK(LL,5) 

! T.Li : Hertz vs linear spring-dashpot contact model
                     IF (TRIM(DES_COLL_MODEL) .EQ. 'HERTZIAN') THEN
                        sqrt_overlap = SQRT(OVERLAP_N)
                        KN_DES_W = hert_kwn(phaseLL)*sqrt_overlap
                        KT_DES_W = hert_kwt(phaseLL)*sqrt_overlap
                        sqrt_overlap = SQRT(sqrt_overlap)
                        ETAN_DES_W = DES_ETAN_WALL(phaseLL)*sqrt_overlap
                        ETAT_DES_W = DES_ETAT_WALL(phaseLL)*sqrt_overlap
                     ELSE
                        KN_DES_W = KN_W
                        KT_DES_W = KT_W
                        ETAN_DES_W = DES_ETAN_WALL(phaseLL)
                        ETAT_DES_W = DES_ETAT_WALL(phaseLL)
                     ENDIF

! Calculate the normal contact force
                     FNS1(:) = -KN_DES_W * OVERLAP_N*NORMAL(:)
                     FNS2(:) = -ETAN_DES_W * V_REL_TRANS_NORM*NORMAL(:)
                     FN(LL,:) = FNS1(:) + FNS2(:) 

! Calculate the tangential displacement which is integration of
! tangential relative velocity with respect to contact time. 
! Correction in the tangential direction is imposed
! Contact T.Li for details                     
! ---------------------------------------------------------------->>>
! Old procedure
!                  PFT(LL,NI,:) = PFT(LL,NI,:)+OVERLAP_T*TANGENT(:)
!                  PFT_TMP(:) = PFT(LL,NI,:)    ! update pft_tmp before it used 
!                  PFT_TMP(:) = PFT(LL,NI,:) - &
!                     DES_DOTPRDCT(PFT_TMP,NORMAL)*NORMAL(:)

! New procedure: van der Hoef et al. (2006)
                     sigmat_old(:) = pft(ll,ni,:)
                     norm_old(:) = pfn(ll,ni,:)
! calculate the unit vector for axis of rotation 
                     if(dimn.eq.3)then
                        call des_crossprdct(tmp_ax,norm_old,normal)
                        tmp_mag=des_dotprdct(tmp_ax,tmp_ax)
                        if(tmp_mag .gt. zero)then
                           tmp_ax(:)=tmp_ax(:)/sqrt(tmp_mag)
! get the old tangential direction unit vector
                           call des_crossprdct(tang_old,tmp_ax,norm_old)
! get the new tangential direction unit vector due to rotation
                           call des_crossprdct(tang_new,tmp_ax,normal)
                           sigmat(:)=des_dotprdct(sigmat_old,tmp_ax)*tmp_ax(:) &
                              + des_dotprdct(sigmat_old,tang_old)*tang_new(:)
                           sigmat(:)=sigmat(:)+overlap_t*tangent(:) 
                        else
                           sigmat(:)=sigmat_old(:)+overlap_t*tangent(:)
                        endif
                     else
                        tang_old(1) =-norm_old(2)
                        tang_old(2) = norm_old(1)
                        tang_new(1) =-normal(2)
                        tang_new(2) = normal(1)
                        sigmat(:)=des_dotprdct(sigmat_old,tang_old)*tang_new(:)
                        sigmat(:)=sigmat(:)+overlap_t*tangent(:) 
                     endif 
  
                     pft_tmp(:) =sigmat(:)
! ----------------------------------------------------------------<<<

! Calculate the tangential contact force
                     FTS1(:) = -KT_DES_W * PFT_TMP(:)
                     FTS2(:) = -ETAT_DES_W * V_REL_TRANS_TANG * TANGENT(:)
                     FT(LL,:) = FTS1(:) + FTS2(:) 

! Temporary storage of tangential contact force for reporting
                     FT_TMP(:) = FT(LL,:)
                                                     
! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall
                     CALL CFSLIDEWALL(LL, TANGENT, PARTICLE_SLIDE)
                  
! Calculate the total force FC and torque TOW on a particle in a
! particle-wall collision
                     CALL CFFCTOWALL(LL, NORMAL, DISTMOD)

! Save the old normal direction
                     pfn(ll,ni,:) = normal(:)
                  
! Save the tangential displacement history with the correction of Coulomb's law
                     IF (PARTICLE_SLIDE) THEN
! Since FT might be corrected during the call to cfslide, the tangental 
! displacement history needs to be changed accordingly   
                        PFT(LL,NI,:) = -( FT(LL,:) - FTS2(:) ) / KT_DES_W
                     ELSE
                        PFT(LL,NI,:) = PFT_TMP(:)
                     ENDIF
        
                     IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN
                     IF (.NOT.DES_LOC_DEBUG) THEN
                        DES_LOC_DEBUG = .TRUE.
                        WRITE(*,1000)
                     ENDIF                          
                     WRITE(*,*) '     STIME, DTSOLID = ', S_TIME, DTSOLID
                     WRITE(*,*) '     WALL CONTACT ON WALL =', IW
                     WRITE(*,*) '     ALREADY_NEIGHBOURS = ',&
                        ALREADY_NEIGHBOURS
                     WRITE(*,*) '     DES_VEL = ', DES_VEL_NEW(LL,1:DIMN),&
                        des_radius(LL)*OMEGA_NEW(LL,1)
                     WRITE(*,*) '     V-OMEGA R = ', &
                        DES_VEL_NEW(LL,1)+des_radius(LL)* OMEGA_NEW(LL,1),&
                        (DES_VEL_NEW(LL,1)+des_radius(LL)*OMEGA_NEW(LL,1))*DTSOLID
                     WRITE(*,*) '     M*g = ', PMASS(LL)*gravity
                     WRITE(*,*) '     KN_W, ETAN_W, KT_W, ETAT_W = ',&
                        KN_DES_W, ETAN_DES_W, KT_DES_W, ETAT_DES_W
                     WRITE(*,*) '     TANGENT= ', TANGENT
                     WRITE(*,*) '     HIST = ', PFT(LL,NI,1:2)
                     WRITE(*,*) '     PARTICLE_SLIDE ? ', PARTICLE_SLIDE
                     WRITE(*,*) '     FT and FN= ', FT( LL,:), FN(LL,:)
                     WRITE(*,*) '     KT_W*OT*TAN = ', &
                        KT_DES_W*((OVERLAP_T)) *TANGENT(:)
                     WRITE(*,*) '     OVERLAP_T = ', OVERLAP_T, TANGENT
                     FTMD = SQRT(DES_DOTPRDCT(FT_TMP,FT_TMP))
                     FNMD = SQRT(DES_DOTPRDCT(FN(LL,1:DIMN),FN(LL,1:DIMN)))
                     WRITE(*,*) '     FTMD, mu FNMD = ', FTMD, MEW_W*FNMD
                     ENDIF

                     PARTICLE_SLIDE = .FALSE.

                  ENDIF   ! IF(R_LM - DISTMOD.GT.SMALL_NUMBER) -> resolve collision                     

              ENDIF   ! end if (wallcontact.eq.1)
           ENDDO   ! end do iw = 1,nwalls
        ENDIF   ! endif(walldtsplit .and. .not.pea(LL,2) and .not. pea(ll,3))

! End check particle LL for wall contacts         
!-----------------------------------------------------------------<<<


! Check particle LL neighbour contacts         
!----------------------------------------------------------------->>>
         IF (NEIGHBOURS(LL,1).GT.0) THEN
            DO II = 2, NEIGHBOURS(LL,1)+1
               I = NEIGHBOURS(LL,II)
               IF(PEA(I,1)) THEN 

                  ALREADY_NEIGHBOURS=.FALSE.
                  
                  IF(PN(LL,1).GT.0) THEN                     
                     DO NEIGH_L = 2, PN(LL,1)+1
                        IF(I.EQ. PN(LL,NEIGH_L)) THEN 
                           ALREADY_NEIGHBOURS=.TRUE.
                           NI = NEIGH_L
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF

                  R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
                  DIST(:) = DES_POS_NEW(I,:) - DES_POS_NEW(LL,:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))

! compute particle-particle VDW cohesive short-range forces	
                  IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
                    EQ_RADIUS = 2d0 * DES_RADIUS(LL)*DES_RADIUS(I) / &
                               (DES_RADIUS(LL)+DES_RADIUS(I))  ! for use in cohesive force
                    DistApart = (DISTMOD-R_LM) ! distance between particle surface
                    IF(DistApart < VDW_OUTER_CUTOFF)THEN
                      IF(DistApart > VDW_INNER_CUTOFF)THEN
                         FORCE_COH = HAMAKER_CONSTANT * EQ_RADIUS / (12d0*DistApart**2) * &
                            ( Asperities/(Asperities+EQ_RADIUS) + &
                            ONE/(ONE+Asperities/DistApart)**2 )
                      ELSE

                         FORCE_COH = 2d0 * PI * SURFACE_ENERGY * EQ_RADIUS * &
                            ( Asperities/(Asperities+EQ_RADIUS) + &
                            ONE/(ONE+Asperities/VDW_INNER_CUTOFF)**2 )
                      ENDIF                       
                      DO IDIMN=1,DIMN
                        Norm_Dist = DIST(IDIMN)/DISTMOD
                        Fcohesive(LL, IDIMN) = Fcohesive(LL, IDIMN) + Norm_Dist*FORCE_COH
                      ENDDO 
                    ENDIF     
                  ENDIF ! for using VDW cohesion model                  
                  
                  IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN

                     IF(DEBUG_DES .AND. LL.EQ.FOCUS_PARTICLE) THEN
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000) 
                        ENDIF
                        WRITE(*,'(5X,A,(10I10))') 'NEIGHBORS: ', NEIGHBOURS(LL,:)
                     ENDIF

                     
                     FRAC_OVERLAP1 = (R_LM-DISTMOD)/DES_RADIUS(LL)
                     FRAC_OVERLAP2 = (R_LM-DISTMOD)/DES_RADIUS(I)
                     IF (FRAC_OVERLAP1 > flag_overlap .OR. &
                         FRAC_OVERLAP2 > flag_overlap) THEN
                        WRITE(*,'(5X,A,A,ES15.7)') &
                           'WARNING: excessive overlap detected ', &
                           'at time ', S_TIME
                        WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
                           'between particles ', LL, 'and ',&
                           I, 'with'
                        WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7,2X,ES15.7)') &
                           'overlap = ', (R_LM-DISTMOD), &
                           ' radii = ', DES_RADIUS(LL), DES_RADIUS(I)
                     ENDIF

                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF
                        WRITE(*,'(5X,A,I10,I10)') &
                           'DISTMOD is zero between particle-pair ',&
                           LL, I
                        STOP
                     ENDIF

! Calculate the components of translational relative velocity for a 
! contacting particle pair and the tangent to the plane of contact
                     CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, &
                        V_REL_TRANS_TANG, TANGENT, NORMAL, DISTMOD)

! Overlap calculation changed from history based to current position 
                     OVERLAP_N = R_LM-DISTMOD

                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE 
                        IF(DEBUG_DES) THEN
                           IF (.NOT.DES_LOC_DEBUG) THEN
                              DES_LOC_DEBUG = .TRUE.
                              WRITE(*,1000)
                           ENDIF
                           WRITE(*,'(5X,A,2(I10,X),A,ES15.7)') &
                              'Normal overlap for particle pair ',&
                              LL, I, ' : ', OVERLAP_N 
                        ENDIF
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1

                        IF (V_REL_TRANS_NORM .GT. ZERO) THEN
                           DTSOLID_TMP = OVERLAP_N/(V_REL_TRANS_NORM)
                        ELSEIF (V_REL_TRANS_NORM .LT. ZERO) THEN
                          DTSOLID_TMP = DTSOLID
! since they are not already neighbors PFN will not have any normal
! information stored so it needs to be calculated
                          DIST_OLD(:)=DES_POS_OLD(I,:)-DES_POS_OLD(LL,:)
                          DISTMOD_OLD=SQRT(DES_DOTPRDCT(DIST_OLD,DIST_OLD))
                          NORMAL_OLD(:)=DIST_OLD(:)/DISTMOD_OLD
                          VRN_OLD(:)=DES_VEL_OLD(LL,:)-DES_VEL_OLD(I,:)
                          V_REL_NORM_OLD=DES_DOTPRDCT(VRN_OLD,NORMAL_OLD)
                          WRITE(*,'(5X,A,A,ES15.7)') &
                             'WARNING: normal relative velocity less ',&
                             'than zero at time ', S_TIME
                          WRITE(*,'(7X,A,I10,2X,A,I10,2X,A)') &
                             'for first contact between particles', LL, &
                             'and ', I, 'with'
                          WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7)') &
                             'V_REL_NORM = ', V_REL_TRANS_NORM, &
                             'and V_REL_NORM_OLD = ', V_REL_NORM_OLD
                        ELSE
                           DTSOLID_TMP = OVERLAP_N/&
                              (V_REL_TRANS_NORM+SMALL_NUMBER)
                        ENDIF
                        OVERLAP_T = V_REL_TRANS_TANG*&
                           MIN(DTSOLID,DTSOLID_TMP)
                     ENDIF
                  ELSE
                     GOTO 300
                  ENDIF

                  phaseLL = PIJK(LL,5)                  
                  phaseI = PIJK(I,5)

! T.Li : Hertz vs linear spring-dashpot contact model
                  IF (TRIM(DES_COLL_MODEL) .EQ. 'HERTZIAN') THEN
                     sqrt_overlap = SQRT(OVERLAP_N)
                     KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap
                     KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap
                     sqrt_overlap = SQRT(sqrt_overlap)
                     ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap
                     ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap
                  ELSE
                     KN_DES = KN
                     KT_DES = KT
                     ETAN_DES = DES_ETAN(phaseLL,phaseI)
                     ETAT_DES = DES_ETAT(phaseLL,phaseI)
                  ENDIF

! Calculate the normal contact force                  
                  FNS1(:) = -KN_DES * OVERLAP_N * NORMAL(:)
                  FNS2(:) = -ETAN_DES * V_REL_TRANS_NORM*NORMAL(:)
                  FN(LL,:) = FNS1(:) + FNS2(:)       

! Contact T.Li for details
! Calculate the tangential displacement which is integration of 
! tangential relative velocity with respect to contact time. 
! Correction in the tangential direction is imposed                   
! ---------------------------------------------------------------->>>
! Old procedure
!                  PFT(LL,NI,:) = PFT(LL,NI,:) + OVERLAP_T * TANGENT(:)
!                  PFT_TMP(:) = PFT(LL,NI,:)   ! update pft_tmp before it is used
!                  PFT_TMP(:) = PFT(LL,NI,:) - &
!                     DES_DOTPRDCT(PFT_TMP,NORMAL)*NORMAL(:)

! New procedure: van der Hoef et al. (2006)
                  sigmat_old(:) = pft(ll,ni,:)
                  norm_old(:) = pfn(ll,ni,:)
! calculate the unit vector for axis of rotation
                  if(dimn.eq.3)then
                     call des_crossprdct(tmp_ax,norm_old,normal)
                     tmp_mag=des_dotprdct(tmp_ax,tmp_ax)
                     if(tmp_mag .gt. zero)then
                        tmp_ax(:)=tmp_ax(:)/sqrt(tmp_mag)
! get the old tangential direction unit vector
                        call des_crossprdct(tang_old,tmp_ax,norm_old)
! get the new tangential direction unit vector due to rotation
                        call des_crossprdct(tang_new,tmp_ax,normal)
                        sigmat(:)=des_dotprdct(sigmat_old,tmp_ax)*tmp_ax(:) &
                              + des_dotprdct(sigmat_old,tang_old)*tang_new(:)
                        sigmat(:)=sigmat(:)+overlap_t*tangent(:)
                     else
                        sigmat(:)=sigmat_old(:)+overlap_t*tangent(:)
                     endif
                  else
                     tang_old(1) =-norm_old(2)
                     tang_old(2) = norm_old(1)
                     tang_new(1) =-normal(2)
                     tang_new(2) = normal(1)
                     sigmat(:)=des_dotprdct(sigmat_old,tang_old)*tang_new(:)
                     sigmat(:)=sigmat(:)+overlap_t*tangent(:)
                  endif     
      
                  pft_tmp(:) =sigmat(:)
! ----------------------------------------------------------------<<<

! Calculate the tangential contact force                  
                  FTS1(:) = -KT_DES * PFT_TMP(:)
                  FTS2(:) = -ETAT_DES * V_REL_TRANS_TANG * TANGENT(:)
                  FT(LL,:) = FTS1(:) + FTS2(:) 

                  FT_TMP(:) = FT(LL,:)

! Check for Coulombs friction law and limit the maximum value of the 
! tangential force on a particle in contact with another particle/wall 
                  CALL CFSLIDE(LL,TANGENT,PARTICLE_SLIDE)
                  
! Calculate the total force FC and torque TOW on a particle in a 
! particle-particle collision
                  CALL CFFCTOW(LL, I, NORMAL, DISTMOD)

! Save the old normal direction
                  PFN(LL,NI,:) = NORMAL(:)

! Save the tangential displacement history with the correction of
! Coulomb's law
                  IF (PARTICLE_SLIDE) THEN
! Since FT might be corrected during the call to cfslide, the tangental 
! displacement history needs to be changed accordingly                          
                     PFT(LL,NI,:) = -( FT(LL,:) - FTS2(:) ) / KT_DES
                  ELSE
                     PFT(LL,NI,:) = PFT_TMP(:)
                  ENDIF
                  
                  IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
                     IF (.NOT.DES_LOC_DEBUG) THEN
                        DES_LOC_DEBUG = .TRUE.
                        WRITE(*,1000) 
                     ENDIF

                     PRINT*, '     EtaN, EtaT =  ', ETAN_DES, ETAT_DES
                     PRINT*, '     Percent overlap = ', (R_LM - DISTMOD)*100.d0/R_LM
                     PRINT*, '     rad ratio = ', DES_RADIUS(LL)/DES_RADIUS(I)
                     PRINT*, '     FNS1 and FNS2 = ', FNS1(:), FNS2(:)
                     PRINT*, '     PFT = ', PFT(LL,NI,:)
                     PRINT*, '     FORCEST = ', FT(LL,:)
                     PRINT*, '     FORCESN = ', FN(LL,:)
                     PRINT*, '     FORCEST = ', FT(LL,:)
                  ENDIF

                  IF(DEBUG_DES.AND.LL.eq.FOCUS_PARTICLE)THEN
                     INQUIRE(FILE='debug_file',EXIST=ALREADY_EXISTS)
                     IF(ALREADY_EXISTS)THEN
                        OPEN(UNIT=1,FILE='debug_file',STATUS='OLD',POSITION='APPEND')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FN(LL,1), 'FNy=',FN(LL,2)
                     ELSE
                        OPEN(UNIT=1,FILE='debug_file',STATUS='NEW')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FN(LL,1),'FNy=',FN(LL,2)
                     ENDIF
                     CLOSE (1)
                     PRINT*, 'PN', PN(LL,:)
                  ENDIF

                  PARTICLE_SLIDE = .FALSE.

               ENDIF         ! IF (I>LL .AND. PEA(I,1))

 300           CONTINUE
            ENDDO            ! DO II = 2, NEIGHBOURS(LL,1)+I
         ENDIF               ! IF(NEIGHBOURS(LL,1).GT.0)

! End check particle LL neighbour contacts         
! ----------------------------------------------------------------<<<

      ENDDO   
!$omp end parallel      
! end loop over paticles LL
! ----------------------------------------------------------------<<<

! Calculate gas-solids drag force on particle 
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL CALC_DES_DRAG_GS
      ENDIF

! Calculate solids-solids drag force on particle 
      IF (DES_CONTINUUM_HYBRID) THEN
         CALL CALC_DES_DRAG_SS
      ENDIF

! The square-well model is still available in the model/cohesion directory.
! However, the VDW routine is now integrated into main mfix DEM to take
! advantage of other DEM capabilities (parallel/search grid)      
      IF(USE_COHESION .AND. .NOT.VAN_DER_WAALS)THEN
         CALL CALC_COHESIVE_FORCES
      ENDIF

! just for post-processing mag. of cohesive forces on each particle
      IF(USE_COHESION)THEN
      magGravity = DSQRT(DES_DOTPRDCT(GRAV,GRAV))
      DO LL = 1, MAX_PIP
         if(.not.pea(ll,1)) cycle 
         if(pea(ll,4)) cycle 
          DO IDIMN=1,DIMN
            PostCohesive(LL) =  PostCohesive(LL) + Fcohesive(LL, IDIMN)**2
          ENDDO
          if(magGravity> ZERO) PostCohesive(LL) =  DSQRT(PostCohesive(LL)) / &
                  (PMASS(LL)*magGravity)
        ENDDO 
      ENDIF ! for cohesion model      
      
! Update the old values of particle position and velocity with the new
! values computed
      CALL CFUPDATEOLD

      IF (DES_LOC_DEBUG) WRITE(*,1001)

 1000 FORMAT(5X,'---------- START CALC_FORCE_DES ---------->')
 1001 FORMAT(5X,'<---------- END CALC_FORCE_DES ----------') 

      RETURN
      END SUBROUTINE CALC_FORCE_DES

