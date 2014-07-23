!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_FORCE_DEM                                          !
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

      SUBROUTINE CALC_FORCE_DEM

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
! normal and tangential forces
      DOUBLE PRECISION :: FNS1(3), FNS2(3)
      DOUBLE PRECISION :: FTS1(3), FTS2(3)
! temporary storage of tangential DISPLACEMENT
      DOUBLE PRECISION :: PFT_TMP(3)

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

!$      double precision omp_start, omp_end
!$      double precision omp_get_wtime

!-----------------------------------------------
! Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
!-----------------------------------------------

! initialize local variables
      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1


! initialize cohesive forces
      IF(USE_COHESION) THEN
         Fcohesive(:,:) = ZERO
         PostCohesive (:) = ZERO
      ENDIF

!$omp   parallel default(shared)                                  &
!$omp   private(ll,fts1,fts2,fns1,fns2,pft_tmp,                   &
!$omp          PARTICLE_SLIDE,nlim,                               &
!$omp          n_nocon,ni,iw,wallcontact,i,                       &
!$omp          already_neighbours, neigh_l,                       &
!$omp          WALL_POS,WALL_VEL,                                 &
!$omp          r_lm,dist,distmod,normal,ii,                       &
!$omp          v_rel_trans_norm, v_rel_tang,                      &
!$omp          overlap_n,overlap_t,dtsolid_tmp,phasell,           &
!$omp          sqrt_overlap,kn_des_w,kt_des_w,etan_des_w,         &
!$omp          etat_des_w,norm_old,tmp_ax,tmp_mag,                &
!$omp          phasei,kn_des,kt_des,etan_des,etat_des,            &
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

! Initializing local variables
         V_REL_TANG(:) = ZERO
         NORMAL(:) = ZERO
         FTS1(:) = ZERO
         FTS2(:) = ZERO
         FNS1(:) = ZERO
         FNS2(:) = ZERO
         PFT_TMP(:) = ZERO
         PARTICLE_SLIDE = .FALSE.

! Check neighbor history of particle LL and update arrays as needed
! ---------------------------------------------------------------->>>
         IF(PN(1,LL).GE.1) THEN
            NLIM = PN(1,LL)+1
            N_NOCON = 0
! For each particle listed as in contact with particle LL in array PN,
! check the flag array PV to determine if particles remained in contact
! after the previous call of CALC_FORCE_DES.
            DO NI = 2, NLIM
               IF(.NOT.PV(NI-N_NOCON,LL)) THEN
! For each particle in PN(2:MAXNEIGHBORS,LL) that is no longer in
! contact, shift the remaining particle contact information PN, PN,
! PFT left by one and reduce PN(1,LL) by one.
                  PN((NI-N_NOCON):MAXNEIGHBORS-1,LL) = &
                     PN((NI-N_NOCON+1):MAXNEIGHBORS,LL)
                  PV((NI-N_NOCON):(MAXNEIGHBORS-1),LL) = &
                     PV((NI-N_NOCON+1):MAXNEIGHBORS,LL)
                  PFT(LL,(NI-N_NOCON):(MAXNEIGHBORS-1),:) = &
                     PFT(LL,(NI-N_NOCON+1):MAXNEIGHBORS,:)
! Save the normal direction at previous time step
                  PFN(LL,(NI-N_NOCON):(MAXNEIGHBORS-1),:) = &
                     PFN(LL,(NI-N_NOCON+1):MAXNEIGHBORS,:)
                  N_NOCON = N_NOCON + 1
                  PN(1,LL) = PN(1,LL) - 1
               ENDIF
            ENDDO

            ! Initializing rest of the neighbor list which is not in contact and
            ! clean up after the above array left shifts
            IF (N_NOCON .GT. 0) THEN
               NLIM = MAX(2,PN(1,LL) + 2)
               PN(NLIM:MAXNEIGHBORS,LL) = -1
               PFT(LL,NLIM:MAXNEIGHBORS,:) = ZERO
               PFN(LL,NLIM:MAXNEIGHBORS,:) = ZERO
            ENDIF
            NEIGH_MAX = MAX(NEIGH_MAX,PN(1,LL))
         ENDIF


! Reset the flag array PV; during each call to calc_force_des this
! variable tracks whether particle LL has any current neighbors
! the array is used in the next call to calc_force_des to update
! particle LL neighbor history above
         PV(2:MAXNEIGHBORS,LL) = .FALSE.

      ENDDO
!$omp end parallel

      CALL CALC_COLLISION_WALL

! Check particle LL neighbour contacts
!----------------------------------------------------------------->>>
      DO LL = 1, MAX_PIP

! skipping non-existent particles
         IF(.NOT.PEA(LL,1)) CYCLE
! skipping ghost particles
         IF(PEA(LL,4)) CYCLE

         DO II = 2, NEIGHBOURS(LL,1)+1
            I = NEIGHBOURS(LL,II)

            !call calc_coll_force(I,LL)

            IF(.NOT.PEA(I,1)) CYCLE

            ALREADY_NEIGHBOURS=.FALSE.

            IF(PN(1,LL).GT.0) THEN
               DO NEIGH_L = 2, PN(1,LL)+1
                  IF(I.EQ. PN(NEIGH_L,LL)) THEN
                     ALREADY_NEIGHBOURS=.TRUE.
                     NI = NEIGH_L
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

            R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
            DIST(:) = DES_POS_NEW(I,:) - DES_POS_NEW(LL,:)
            DISTMOD = DES_DOTPRDCT(DIST,DIST)

            ! compute particle-particle VDW cohesive short-range forces
            IF(USE_COHESION .AND. VAN_DER_WAALS) call calc_cohesion()

            IF(DISTMOD > (R_LM + SMALL_NUMBER)**2) CYCLE
            IF(DISTMOD == 0) STOP
            DISTMOD = SQRT(DISTMOD)
            NORMAL(:)= DIST(:)/DISTMOD

            ! Overlap calculation changed from history based to current position
            OVERLAP_N = R_LM-DISTMOD

            ! Calculate the components of translational relative velocity for a
            ! contacting particle pair and the tangent to the plane of contact
            CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, &
                 V_REL_TANG, NORMAL, DISTMOD)

            IF(ALREADY_NEIGHBOURS) THEN
               PV(NI,LL) = .TRUE.
               OVERLAP_T = DTSOLID*V_REL_TANG
            ELSE
               PN(1,LL) = PN(1,LL) + 1
               NI = PN(1,LL) + 1
               PN(NI,LL) = I
               PV(NI,LL) = .TRUE.

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
            phaseI = PIJK(I,5)

            ! T.Li : Hertz vs linear spring-dashpot contact model
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
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
            FN(:,LL) = FNS1(:) + FNS2(:)

            call calc_tangential_displacement

            ! Calculate the tangential contact force
            FTS1(:) = -KT_DES * PFT_TMP(:)
            FTS2(:) = -ETAT_DES * V_REL_TANG
            FT(:,LL) = FTS1(:) + FTS2(:)

            ! Check for Coulombs friction law and limit the maximum value of the
            ! tangential force on a particle in contact with another particle/wall
            CALL CFSLIDE(LL,V_REL_TANG,PARTICLE_SLIDE,MEW,FT(:,LL),FN(:,LL))

            ! Calculate the total force FC and torque TOW on a particle in a
            ! particle-particle collision
            CALL CFFCTOW(LL, I, NORMAL, DISTMOD, FC(:,LL), FN(:,LL), FT(:,LL), TOW(:,LL))

            ! Save the tangential displacement history with the correction of
            ! Coulomb's law
            IF (PARTICLE_SLIDE) THEN
               ! Since FT might be corrected during the call to cfslide, the tangental
               ! displacement history needs to be changed accordingly
               PFT(LL,NI,:) = -( FT(:,LL) - FTS2(:) ) / KT_DES
               PARTICLE_SLIDE = .FALSE.
            ELSE
               PFT(LL,NI,:) = PFT_TMP(:)
            ENDIF

         ENDDO            ! DO II = 2, NEIGHBOURS(LL,1)+I
   ENDDO

   ! Calculate gas-solids drag force on particle
   IF(DES_CONTINUUM_COUPLED) THEN
      CALL CALC_DES_DRAG_GS
   ENDIF

   ! Calculate solids-solids drag force on particle
   IF (DES_CONTINUUM_HYBRID) THEN
      CALL CALC_DES_DRAG_SS
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

      RETURN

      contains

        subroutine calc_coll_force(I,LL)
      IMPLICIT NONE

      ! particle indices
      INTEGER :: I, LL

        end subroutine calc_coll_force


        subroutine calc_cohesion

          implicit none

          IF(DISTMOD < (R_LM+VDW_OUTER_CUTOFF)**2) THEN
             EQ_RADIUS = 2d0 * DES_RADIUS(LL)*DES_RADIUS(I) / &
                  (DES_RADIUS(LL)+DES_RADIUS(I))  ! for use in cohesive force
             IF(DISTMOD < (R_LM+VDW_INNER_CUTOFF)**2) THEN
                DistApart = (SQRT(DISTMOD)-R_LM) ! distance between particle surface
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
        end subroutine CALC_COHESION

        subroutine calc_tangential_displacement

          implicit none

! tangent to the plane of contact at current and previous time step
! (used for 2D calculations)
      DOUBLE PRECISION :: TANG_OLD(3),TANG_NEW(3)

! local variables for the accumulated tangential displacement that occurs
! for the particle-particle or particle-wall collision (current time
! step and previous time step)
      DOUBLE PRECISION :: SIGMAT(3),SIGMAT_OLD(3)

                  IF(USE_VDH_DEM_MODEL) then
! Calculate the tangential displacement which is integration of
! tangential relative velocity with respect to contact time.
! Correction in the tangential direction is imposed

! New procedure: van der Hoef et al. (2006)
                     sigmat_old(:) = pft(ll,ni,:)
                     norm_old(:) = pfn(ll,ni,:)
! calculate the unit vector for axis of rotation
                     if(DO_K)then
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
                           sigmat(:)=sigmat(:)+overlap_t
                        else
                           sigmat(:)=sigmat_old(:)+overlap_t
                        endif
                     else
                        tang_old(1) =-norm_old(2)
                        tang_old(2) = norm_old(1)
                        tang_new(1) =-normal(2)
                        tang_new(2) = normal(1)
                        sigmat(:)=des_dotprdct(sigmat_old,tang_old)*tang_new(:)
                        sigmat(:)=sigmat(:)+overlap_t
                     endif

! Save the old normal direction
                     PFN(LL,NI,:) = NORMAL(:)
                     PFT_TMP(:)   = SIGMAT(:)
                  ELSE
                     ! Old procedure
                     PFT(LL,NI,:) = PFT(LL,NI,:) + OVERLAP_T
                     PFT_TMP(:) = PFT(LL,NI,:) ! update pft_tmp before it is used
                     PFT_TMP(:) = PFT(LL,NI,:) - &
                     DES_DOTPRDCT(PFT_TMP,NORMAL)*NORMAL(:)
                  ENDIF

        end subroutine calc_tangential_displacement

      END SUBROUTINE CALC_FORCE_DEM

