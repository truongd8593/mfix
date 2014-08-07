!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_FORCE_DEM                                          !
!  Author: J.Musser                                   Date: 30-Apr-14  !
!                                                                      !
!  Purpose: Driver routine for switching between particle-wall         !
!           collision algorithms.                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_COLLISION_WALL

      USE softspring_funcs_cutcell
      use discretelement, only: USE_STL_DES

      IF(USE_STL_DES) THEN
         CALL CALC_DEM_FORCE_WITH_WALL_STL
      ELSE
         CALL CALC_FORCE_DEM_WALL
      ENDIF


      RETURN
      END SUBROUTINE CALC_COLLISION_WALL




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_FORCE_DEM                                          !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Calculate contact force and torque on particle from        !
!           particle-particle and particle-wall collisions. Treats     !
!           wall interaction also as a two-particle interaction but    !
!           accounting for the wall properties                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_FORCE_DEM_WALL
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE cutcell

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
      DOUBLE PRECISION :: OVERLAP_N
      DOUBLE PRECISION :: OVERLAP_T(3)
! square root of the overlap
      DOUBLE PRECISION :: SQRT_OVERLAP
      DOUBLE PRECISION :: FRAC_OVERLAP1, FRAC_OVERLAP2
! distance vector between two particle centers or between a particle
! center and wall when the two surfaces are just at contact (i.e. no
! overlap)
      DOUBLE PRECISION :: R_LM
! the normal and tangential components of the translational relative
! velocity
      DOUBLE PRECISION :: V_REL_TRANS_NORM
! time elapsed to travel the calculated normal overlap given the
! normal relative velocity
      DOUBLE PRECISION :: DTSOLID_TMP
! distance vector between two particle centers or between a particle
! center and wall at current and previous time steps
      DOUBLE PRECISION :: DIST(3)
! magnitude of distance between two particle centers or between a
! particle center and wall at current and previous time steps
      DOUBLE PRECISION :: DISTMOD
! unit normal vector along the line of contact between contacting
! particles or particle-wall at current and previous time steps
      DOUBLE PRECISION :: NORMAL(3), NORM_OLD(3)
! tangent to the plane of contact at current time step
      DOUBLE PRECISION :: V_REL_TANG(3)
! variables for tangential displacement calculation:
! unit vector for axis of rotation and its magnitude
      DOUBLE PRECISION :: TMP_AX(3), TMP_MAG
! local variables for the accumulated tangential displacement that occurs
! for the particle-particle or particle-wall collision (current time
! step and previous time step)
      DOUBLE PRECISION :: SIGMAT(3),SIGMAT_OLD(3)
! tangent to the plane of contact at current and previous time step
! (used for 2D calculations)
      DOUBLE PRECISION :: TANG_OLD(3),TANG_NEW(3)
! normal and tangential forces
      DOUBLE PRECISION :: FNS1(3), FNS2(3)
      DOUBLE PRECISION :: FTS1(3), FTS2(3)
! temporary storage of tangential DISPLACEMENT
      DOUBLE PRECISION :: PFT_TMP(3)
! normal/tangential force
      DOUBLE PRECISION :: FN(3), FT(3)
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
      DOUBLE PRECISION :: WALL_POS(3), WALL_VEL(3)
! set to T when a sliding contact occurs
      LOGICAL :: PARTICLE_SLIDE
! logic flag telling whether contact pair is old (previously detected)
      LOGICAL :: ALREADY_NEIGHBOURS

! logic flag for local debug warnings
      LOGICAL :: DES_LOC_DEBUG
      LOGICAL :: ALREADY_EXISTS

      LOGICAL :: report_excess_overlap
!$      double precision omp_start, omp_end
!$      double precision omp_get_wtime

      report_excess_overlap = .false.
! initialize local variables
      DES_LOC_DEBUG = .FALSE.
      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1




!$omp   parallel default(shared)                                  &
!$omp   private(ll,fts1,fts2,fns1,fns2,pft_tmp,                   &
!$omp          PARTICLE_SLIDE,nlim,                               &
!$omp          n_nocon,ni,iw,wallcontact,i,                       &
!$omp          already_neighbours, report_excess_overlap,neigh_l, &
!$omp          WALL_POS,WALL_VEL,                                 &
!$omp          r_lm,dist,distmod,frac_overlap1,normal,            &
!$omp          v_rel_trans_norm,V_REL_TANG,                       &
!$omp          overlap_n,overlap_t,dtsolid_tmp,phasell,           &
!$omp          sqrt_overlap,kn_des_w,kt_des_w,etan_des_w,         &
!$omp          etat_des_w,sigmat_old,norm_old,tmp_ax,tmp_mag,     &
!$omp          tang_old,tang_new,sigmat,                          &
!$omp          ftmd,fnmd,                                         &
!$omp          ii,frac_overlap2,                                  &
!$omp          phasei,kn_des,kt_des,etan_des,etat_des,            &
!!$omp          Idimn,force_coh,eq_radius,distapart,norm_dist,maggravity)
!!$omp do reduction(max:NEIGH_MAX,OVERLAP_MAX) schedule (dynamic,50)
! Revised by Handan Liu on Jan 17 2013
!$omp          Idimn,force_coh,eq_radius,distapart,norm_dist)
!$omp do reduction(max:NEIGH_MAX) schedule (dynamic,50)


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
               'X,Y POS: ', DES_POS_NEW(1,LL), DES_POS_NEW(2,LL)
            WRITE(*,'(5X,A,2(ES15.7))') &
               'X,Y VEL: ', DES_VEL_NEW(1,LL), DES_VEL_NEW(2,LL)
         ENDIF

! Initializing local variables
         V_REL_TANG(:) = ZERO
         NORMAL(:) = ZERO
         FTS1(:) = ZERO
         FTS2(:) = ZERO
         FNS1(:) = ZERO
         FNS2(:) = ZERO
         PFT_TMP(:) = ZERO
         PARTICLE_SLIDE = .FALSE.


! Check particle LL for wall contacts
!----------------------------------------------------------------->>>
! Treats wall interaction also as a two-particle interaction but accounting
! for the wall properties; make sure the particle is not classified as
! a new 'entering' particle or is already marked as a potential exiting particle
         IF(.NOT.PEA(LL,2) .AND. .NOT.PEA(LL,3)) THEN

            DO IW = 1, NWALLS
               WALLCONTACT = 0

! Check to see if a particle is in contact with any of the walls
               CALL CFWALLCONTACT(IW, LL, WALLCONTACT)

               IF(WALLCONTACT.EQ.1) THEN
                  I = MAX_PIP + IW
                  ALREADY_NEIGHBOURS=.FALSE.

                  IF(PN_WALL(1,LL).GT.0) THEN
                     DO NEIGH_L = 2, PN_WALL(1,LL)+1
                        IF(I.EQ. PN_WALL(NEIGH_L,LL)) THEN
                           ALREADY_NEIGHBOURS=.TRUE.
                           NI = NEIGH_L
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF

! Assign the wall particle a position and velocity
                  CALL CFWALLPOSVEL(LL, IW, WALL_POS, WALL_VEL)

                  R_LM = DES_RADIUS(LL) + DES_RADIUS(LL)
                  DIST(:) = WALL_POS(:) - DES_POS_NEW(:,LL)
                  DISTMOD = SQRT(dot_product(DIST,DIST))

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
                           Fcohesive(IDIMN, LL) = Fcohesive(IDIMN, LL) + Norm_Dist*FORCE_COH
                        ENDDO
                     ENDIF
                  ENDIF ! for using VDW cohesion model


                  IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN
! Particle-wall overlap detected (i.e. resolve collision)
!----------------------------------

! Error reporting
                     FRAC_OVERLAP1 = (R_LM-DISTMOD)/DES_RADIUS(LL)
                     IF (FRAC_OVERLAP1 > flag_overlap.and.report_excess_overlap) THEN
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
                        V_REL_TANG, NORMAL, DISTMOD)

! The normal overlap calculation was changed so that it no longer
! depends on the contact history (i.e., integration of incremental
! overlap found by velocity*dtsolid).  Now overlap is based purely on
! position of neighbors.  The change was made because the former
! method was found not to conserve energy
                     OVERLAP_N =  R_LM-DISTMOD

                     IF(ALREADY_NEIGHBOURS) THEN
                        PV_WALL(NI,LL) = .TRUE.
                        OVERLAP_T = DTSOLID*V_REL_TANG
                     ELSE
                        PN_WALL(1,LL) = PN_WALL(1,LL) + 1
                        NI = PN_WALL(1,LL) + 1
                        PN_WALL(NI,LL) = I
                        PV_WALL(NI,LL) = .TRUE.
                        IF (V_REL_TRANS_NORM .GT. ZERO) THEN
                          DTSOLID_TMP = OVERLAP_N/(V_REL_TRANS_NORM)
                        ELSEIF (V_REL_TRANS_NORM .LT. ZERO) THEN
                          DTSOLID_TMP = DTSOLID
                        ELSE
                           DTSOLID_TMP = OVERLAP_N/&
                              (V_REL_TRANS_NORM+SMALL_NUMBER)
                        ENDIF
                        OVERLAP_T = MIN(DTSOLID,DTSOLID_TMP)*V_REL_TANG
                     ENDIF


                     phaseLL = PIJK(LL,5)

! T.Li : Hertz vs linear spring-dashpot contact model
                     IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
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
                     FN(:) = FNS1(:) + FNS2(:)

                     IF(USE_VDH_DEM_MODEL) then

! Calculate the tangential displacement which is integration of
! tangential relative velocity with respect to contact time.
! Correction in the tangential direction is imposed

! New procedure: van der Hoef et al. (2006)
                        sigmat_old(:) = pft_WALL(ll,ni,:)
                        norm_old(:)   = pfn_WALL(ll,ni,:)
! calculate the unit vector for axis of rotation
                        if(DO_K)then
                           call des_crossprdct(tmp_ax,norm_old,normal)
                           tmp_mag   = dot_product(tmp_ax,tmp_ax)
                           if(tmp_mag .gt. zero)then
                              tmp_ax(:) = tmp_ax(:)/sqrt(tmp_mag)
! get the old tangential direction unit vector
                              call des_crossprdct(tang_old,tmp_ax,norm_old)
! get the new tangential direction unit vector due to rotation
                              call des_crossprdct(tang_new,tmp_ax,normal)
                                 sigmat(:) = dot_product(sigmat_old,tmp_ax)*tmp_ax(:) &
                                 + dot_product(sigmat_old,tang_old)*tang_new(:)

                              sigmat(:) = sigmat(:)+overlap_t(:)
                           else
                              sigmat(:) = sigmat_old(:)+overlap_t(:)
                           endif
                        else
                           tang_old(1) = -norm_old(2)
                           tang_old(2) =  norm_old(1)
                           tang_new(1) = -normal(2)
                           tang_new(2) =  normal(1)
                           sigmat(:)= dot_product(sigmat_old,tang_old)*tang_new(:)
                           sigmat(:)= sigmat(:)+overlap_t(:)
                        endif

                        pft_tmp(:) = sigmat(:)
            ! Save the old normal direction
                        pfn_WALL(ll,ni,:)   = normal(:)
                     else ! Old procedure
                        PFT_WALL(LL,NI,:) = PFT_WALL(LL,NI,:) + OVERLAP_T(:)
                        PFT_TMP(:) = PFT_WALL(LL,NI,:) ! update pft_tmp before it used
                     !remove the normal component from the tangential force
                     !due to change of normal direction
                        PFT_TMP(:)   = PFT_WALL(LL,NI,:) - &
                           dot_product(PFT_TMP,NORMAL)*NORMAL(:)
                     endif
! ----------------------------------------------------------------<<<

! Calculate the tangential contact force
                     FTS1(:)  = -KT_DES_W * PFT_TMP(:)
                     FTS2(:)  = -ETAT_DES_W * V_REL_TANG(:)
                     FT(:) = FTS1(:) + FTS2(:)


! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall
                     CALL CFSLIDE(V_REL_TANG, PARTICLE_SLIDE, MEW_W,FT(:),FN(:))

! Calculate the total force FC and torque TOW on a particle in a
! particle-wall collision
                     CALL CFFCTOWALL(LL, NORMAL, DISTMOD, FN, FT)


! Save the tangential displacement history with the correction of Coulomb's law
                     IF (PARTICLE_SLIDE) THEN
! Since FT might be corrected during the call to cfslide, the tangental
! displacement history needs to be changed accordingly
                        PFT_WALL(LL,NI,:) = -( FT(:) - FTS2(:) ) / KT_DES_W
                     ELSE
                        PFT_WALL(LL,NI,:) = PFT_TMP(:)
                     ENDIF


                     PARTICLE_SLIDE = .FALSE.

                  ENDIF   ! IF(R_LM - DISTMOD.GT.SMALL_NUMBER) -> resolve collision

              ENDIF   ! end if (wallcontact.eq.1)
           ENDDO   ! end do iw = 1,nwalls
        ENDIF   ! endif(.not.pea(LL,2) and .not. pea(ll,3))

      ENDDO
!$omp end parallel


 1000 FORMAT(5X,'---------- START CALC_FORCE_DES ---------->')
 1001 FORMAT(5X,'<---------- END CALC_FORCE_DES ----------')

      RETURN
      END SUBROUTINE CALC_FORCE_DEM_WALL

