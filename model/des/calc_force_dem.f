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
      USE discretelement
      USE geometry, ONLY: DO_K
      USE constant, ONLY: Pi

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0
! particle no. indices
      INTEGER :: I, LL, CC
! loop counters (neighbor particles, walls, dimension)
      INTEGER :: II, JJ, IW, IDIMN
! loop counter and indices for neighbor information
      INTEGER :: NI, NLIM, N_NOCON, NEIGH_L
! the overlap occuring between particle-particle or particle-wall
! collision in the normal direction
      DOUBLE PRECISION :: OVERLAP_N
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
! particles or particle-wall at current time step
      DOUBLE PRECISION :: NORMAL(3)
! tangent to the plane of contact at current time step
      DOUBLE PRECISION :: V_REL_TANG(3)
! normal and tangential forces
      DOUBLE PRECISION :: FN(3), FT(3)
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

      LOGICAL, PARAMETER :: report_excess_overlap = .FALSE.
!$      double precision omp_start, omp_end
!$      double precision omp_get_wtime

!-----------------------------------------------

! initialize local variables
      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1

! initialize cohesive forces
      IF(USE_COHESION) THEN
         PostCohesive (:) = ZERO
      ENDIF

!$omp   parallel default(shared)                                  &
!$omp   private(ll,fts1,fts2,fns1,fns2,pft_tmp,                   &
!$omp          PARTICLE_SLIDE,nlim,                               &
!$omp          n_nocon,ni,iw,wallcontact,i,                       &
!$omp          WALL_POS,WALL_VEL,neigh_l,                         &
!$omp          r_lm,dist,distmod,normal,ii,                       &
!$omp          v_rel_trans_norm, v_rel_tang,                      &
!$omp          overlap_n,dtsolid_tmp,phasell,                     &
!$omp          sqrt_overlap,kn_des_w,kt_des_w,etan_des_w,         &
!$omp          etat_des_w,                                        &
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

! Check neighbor history of particle LL and update arrays as needed
! ---------------------------------------------------------------->>>
         IF(PN_WALL(1,LL).GE.1) THEN
            NLIM = PN_WALL(1,LL)+1
            N_NOCON = 0
! For each particle listed as in contact with particle LL in array PN,
! check the flag array PV to determine if particles remained in contact
! after the previous call of CALC_FORCE_DES.
            DO NI = 2, NLIM
               IF(.NOT.PV_WALL(NI-N_NOCON,LL)) THEN
! For each particle in PN(2:6,LL) that is no longer in
! contact, shift the remaining particle contact information PN, PN,
! PFT left by one and reduce PN(1,LL) by one.
                  PN_WALL((NI-N_NOCON):6-1,LL) = &
                     PN_WALL((NI-N_NOCON+1):6,LL)
                  PV_WALL((NI-N_NOCON):(6-1),LL) = &
                     PV_WALL((NI-N_NOCON+1):6,LL)
                  PFT_WALL(LL,(NI-N_NOCON):(6-1),:) = &
                     PFT_WALL(LL,(NI-N_NOCON+1):6,:)
! Save the normal direction at previous time step
                  PFN_WALL(LL,(NI-N_NOCON):(6-1),:) = &
                     PFN_WALL(LL,(NI-N_NOCON+1):6,:)
                  N_NOCON = N_NOCON + 1
                  PN_WALL(1,LL) = PN_WALL(1,LL) - 1
               ENDIF
            ENDDO

            ! Initializing rest of the neighbor list which is not in contact and
            ! clean up after the above array left shifts
            IF (N_NOCON .GT. 0) THEN
               NLIM = MAX(2,PN_WALL(1,LL) + 2)
               PN_WALL(NLIM:6,LL) = -1
               PFT_WALL(LL,NLIM:6,:) = ZERO
               PFN_WALL(LL,NLIM:6,:) = ZERO
            ENDIF
            NEIGH_MAX = MAX(NEIGH_MAX,PN_WALL(1,LL))
         ENDIF


! Reset the flag array PV; during each call to calc_force_des this
! variable tracks whether particle LL has any current neighbors
! the array is used in the next call to calc_force_des to update
! particle LL neighbor history above
         PV_WALL(1:6,LL) = .FALSE.

      ENDDO
!$omp end parallel

      CALL CALC_COLLISION_WALL

! Check particle LL neighbour contacts
!----------------------------------------------------------------->>>

      FC_COLL(:,:) = 0
      TOW_COLL(:,:) = 0

!$omp parallel do default(none) private(cc,ll,i,dist,distmod,r_lm,normal,overlap_n,v_rel_tang,v_rel_trans_norm,sqrt_overlap,kn_des,kt_des,hert_kn,hert_kt,phasell,phasei,etan_des,etat_des,fns1,fns2,fts1,fts2,pft_tmp,fn,ft,particle_slide,eq_radius,distapart,force_coh) &
!$omp     shared(collisions,collision_num,des_pos_new,des_radius,des_coll_model_enum,kn,kt,pv_coll,pft_coll,pfn_coll,pijk,des_etan,des_etat,mew,fc_coll,tow_coll,use_cohesion,van_der_waals,vdw_outer_cutoff,vdw_inner_cutoff,hamaker_constant,asperities,surface_energy)
      DO CC = 1, COLLISION_NUM
         LL = COLLISIONS(1,CC)
         I  = COLLISIONS(2,CC)

         R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
         DIST(:) = DES_POS_NEW(:,I) - DES_POS_NEW(:,LL)
         DISTMOD = dot_product(DIST,DIST)

! compute particle-particle VDW cohesive short-range forces
         IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
            IF(DISTMOD < (R_LM+VDW_OUTER_CUTOFF)**2) THEN
               EQ_RADIUS = 2d0 * DES_RADIUS(LL)*DES_RADIUS(I) / &
                    (DES_RADIUS(LL)+DES_RADIUS(I))  ! for use in cohesive force
               IF(DISTMOD > (VDW_INNER_CUTOFF+R_LM)**2) THEN
                  DistApart = (SQRT(DISTMOD)-R_LM) ! distance between particle surface
                  FORCE_COH = HAMAKER_CONSTANT * EQ_RADIUS / (12d0*DistApart**2) * &
                       ( Asperities/(Asperities+EQ_RADIUS) + &
                       ONE/(ONE+Asperities/DistApart)**2 )
               ELSE
                  FORCE_COH = 2d0 * PI * SURFACE_ENERGY * EQ_RADIUS * &
                       ( Asperities/(Asperities+EQ_RADIUS) + &
                       ONE/(ONE+Asperities/VDW_INNER_CUTOFF)**2 )
               ENDIF
               FC_COLL(:, LL) = FC_COLL(:, LL) + DIST(:)*FORCE_COH/SQRT(DISTMOD)
            ENDIF
         ENDIF

         IF(DISTMOD > (R_LM + SMALL_NUMBER)**2) THEN
            PV_COLL(CC) = .false.
            PFT_COLL(:,CC) = 0.0
            PFN_COLL(:,CC) = 0.0
            CYCLE
         ENDIF

         IF(DISTMOD == 0) THEN
            WRITE(*,8550) LL, I
            STOP "division by zero"
 8550 FORMAT('DISTMOD is zero betwen particles:',2(2x,I10))
         ENDIF
         DISTMOD = SQRT(DISTMOD)
         NORMAL(:)= DIST(:)/DISTMOD

! Overlap calculation changed from history based to current position
         OVERLAP_N = R_LM-DISTMOD

         IF (report_excess_overlap) call print_excess_overlap

! Calculate the components of translational relative velocity for a
! contacting particle pair and the tangent to the plane of contact
         CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, &
              V_REL_TANG, NORMAL, DISTMOD)

         phaseLL = PIJK(LL,5)
         phaseI = PIJK(I,5)

! Hertz spring-dashpot contact model
         IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
            sqrt_overlap = SQRT(OVERLAP_N)
            KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap
            KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap
            sqrt_overlap = SQRT(sqrt_overlap)
            ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap
            ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap

! Linear spring-dashpot contact model
         ELSE
            KN_DES = KN
            KT_DES = KT
            ETAN_DES = DES_ETAN(phaseLL,phaseI)
            ETAT_DES = DES_ETAT(phaseLL,phaseI)
         ENDIF

! Calculate the normal contact force
         FNS1(:) = -KN_DES * OVERLAP_N * NORMAL(:)
         FNS2(:) = -ETAN_DES * V_REL_TRANS_NORM*NORMAL(:)
         FN(:) = FNS1(:) + FNS2(:)

         call calc_tangential_displacement(pft_tmp(:),normal(:),pfn_coll(:,cc),pft_coll(:,cc),overlap_n,v_rel_trans_norm,v_rel_tang(:),PV_COLL(CC))
         PV_COLL(CC) = .true.

! Calculate the tangential contact force
         FTS1(:) = -KT_DES * PFT_TMP(:)
         FTS2(:) = -ETAT_DES * V_REL_TANG
         FT(:) = FTS1(:) + FTS2(:)

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with another particle/wall
         PARTICLE_SLIDE = .FALSE.
         CALL CFSLIDE(V_REL_TANG(:), PARTICLE_SLIDE, MEW,           &
             FT(:), FN(:))

! Calculate the total force FC and torque TOW on a particle in a
! particle-particle collision
         CALL CFFCTOW(DES_RADIUS(LL), DES_RADIUS(I), NORMAL,        &
            DISTMOD, FC_COLL(:,CC), FN(:), FT(:), TOW_COLL(:,CC))

! Save tangential displacement history with Coulomb's law correction
         IF (PARTICLE_SLIDE) THEN
! Since FT might be corrected during the call to cfslide, the tangential
! displacement history needs to be changed accordingly
            PFT_COLL(:,CC) = -( FT(:) - FTS2(:) ) / KT_DES
         ELSE
            PFT_COLL(:,CC) = PFT_TMP(:)
         ENDIF

ENDDO
!$omp end parallel do

      magGravity = SQRT(dot_product(GRAV,GRAV))

      CC = 1
OUTER: DO LL = 1, MAX_PIP
         IF(.NOT.PEA(LL,1)) CYCLE
         IF(PEA(LL,4)) CYCLE

         II = COLLISIONS(1,CC)
         JJ = COLLISIONS(2,CC)

         if (COLLISION_NUM < CC) EXIT

         DO WHILE (II < LL)
            CC = CC+1
            if (COLLISION_NUM < CC) EXIT OUTER
            II = COLLISIONS(1,CC)
            JJ = COLLISIONS(2,CC)
         ENDDO

         DO WHILE (II .eq. LL)
            FC(:,LL) = FC(:,LL) + FC_COLL(:,CC)
            TOW(:,LL) = TOW(:,LL) + TOW_COLL(:,CC)

! just for post-processing mag. of cohesive forces on each particle
            IF(USE_COHESION)THEN
               PostCohesive(LL) =  dot_product(FC_COLL(:, CC),FC_COLL(:, CC))
               if(magGravity> ZERO) PostCohesive(LL) =  SQRT(PostCohesive(LL)) / &
                    (PMASS(LL)*magGravity)
            ENDIF ! for cohesion model

            CC = CC+1
            if (COLLISION_NUM < CC) EXIT OUTER
            II = COLLISIONS(1,CC)
            JJ = COLLISIONS(2,CC)
         ENDDO
      ENDDO OUTER

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) CALL CALC_DES_DRAG_GS

! Calculate solids-solids drag force on particle
      IF(DES_CONTINUUM_HYBRID) CALL CALC_DES_DRAG_SS

! Update the old values of particle position and velocity with the new
! values computed
      CALL CFUPDATEOLD

      RETURN

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: print_excess_overlap
!  Purpose: Print overlap warnings
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        SUBROUTINE print_excess_overlap
          IF (OVERLAP_N > flag_overlap*DES_RADIUS(LL) .OR. &
               OVERLAP_N > flag_overlap*DES_RADIUS(I)) THEN
             WRITE(*,'(5X,A,A,ES15.7)') &
                  'WARNING: excessive overlap detected ', &
                  'at time ', S_TIME
             WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
                  'between particles ', LL, 'and ',&
                  I, 'with'
             WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7,2X,ES15.7)') &
                  'overlap = ', OVERLAP_N, &
                  ' radii = ', DES_RADIUS(LL), DES_RADIUS(I)
          ENDIF
        END SUBROUTINE print_excess_overlap

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CALC_TANGENTIAL_DISPLACEMENT
!  Purpose: Calculate the tangential displacement which is integration of
!           tangential relative velocity with respect to contact time.
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        SUBROUTINE CALC_TANGENTIAL_DISPLACEMENT(PFT,norm,norm_old,sigmat_old,overlap_norm, relvel_tang_norm,relvel_tang, already_colliding)

          implicit none

! tangential displacement
          DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: pft

! local variables for the accumulated tangential displacement that occurs
! for the particle-particle or particle-wall collision (current time
! step and previous time step)
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: SIGMAT_OLD

      LOGICAL, INTENT(IN) :: already_colliding

! the overlap occuring between particle-particle or particle-wall
! collision in the tangential direction
      DOUBLE PRECISION :: OVERLAP_T(3)

      DOUBLE PRECISION, DIMENSION(3) :: SIGMAT

! unit normal vector along the line of contact between contacting
! particles or particle-wall at previous time step
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: NORM,NORM_OLD,relvel_tang
      DOUBLE PRECISION, INTENT(IN) :: overlap_norm, relvel_tang_norm

! variables for tangential displacement calculation:
! unit vector for axis of rotation and its magnitude
      DOUBLE PRECISION :: TMP_AX(3), TMP_MAG

! tangent to the plane of contact at current and previous time step
! (used for 2D calculations)
      DOUBLE PRECISION :: TANG_OLD(3),TANG_NEW(3)

      IF(already_colliding) THEN
         OVERLAP_T(:) = DTSOLID*relvel_tang(:)
      ELSE
         DTSOLID_TMP = OVERLAP_NORM/MAX(relvel_tang_norm,SMALL_NUMBER)
         OVERLAP_T(:) = MIN(DTSOLID,DTSOLID_TMP)*relvel_tang(:)
      ENDIF

      IF(USE_VDH_DEM_MODEL) then
! Calculate the tangential displacement which is integration of
! tangential relative velocity with respect to contact time.
! Correction in the tangential direction is imposed

! New procedure: van der Hoef et al. (2006)

! calculate the unit vector for axis of rotation
         if(DO_K)then
            call des_crossprdct(tmp_ax,norm_old,norm)
            tmp_mag=dot_product(tmp_ax,tmp_ax)
            if(tmp_mag .gt. zero)then
               tmp_ax(:)=tmp_ax(:)/sqrt(tmp_mag)
! get the old tangential direction unit vector
               call des_crossprdct(tang_old,tmp_ax,norm_old)
! get the new tangential direction unit vector due to rotation
               call des_crossprdct(tang_new,tmp_ax,norm)
               sigmat(:)=dot_product(sigmat_old,tmp_ax)*tmp_ax(:) &
                    + dot_product(sigmat_old,tang_old)*tang_new(:)
               sigmat(:)=sigmat(:)+overlap_t(:)
            else
               sigmat(:)=sigmat_old(:)+overlap_t(:)
            endif
         else
            tang_old(1) =-norm_old(2)
            tang_old(2) = norm_old(1)
            tang_new(1) =-norm(2)
            tang_new(2) = norm(1)
            sigmat(:)=dot_product(sigmat_old,tang_old)*tang_new(:)
            sigmat(:)=sigmat(:)+overlap_t(:)
         endif

! Save the old normal direction
         PFT(:)   = SIGMAT(:)
      ELSE
         ! Old procedure
         sigmat = sigmat_old + OVERLAP_T(:)
         pft(:) = sigmat - dot_product(sigmat,NORM)*NORM(:)
      ENDIF

    END SUBROUTINE CALC_TANGENTIAL_DISPLACEMENT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFFCTOW
!  Purpose: Calculate the total force and torque on a particle
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer: Rahul Garg                               Date: 02-Aug-07
!
!  Comments: Implement eqns 13 & 14 from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFFCTOW(RAD_L, RAD_II,  NORM, DIST_LI, FC_tmp, FN_tmp, FT_tmp, TOW_tmp)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! radii of particle-particle contact pair
      DOUBLE PRECISION, INTENT(IN) :: RAD_L, RAD_II
! distance between particle centers
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI
! unit normal vector along the line of contact pointing from
! particle L to particle II
      DOUBLE PRECISION, INTENT(IN) :: NORM(3)
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: FN_tmp, FT_tmp
      DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT) :: FC_tmp
      DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: TOW_tmp
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! distance from the contact point to the particle center
      DOUBLE PRECISION :: DIST_CL
!------------------------------------------------

! total contact force ( FC_tmp may already include cohesive force)
      FC_tmp(:) = FC_tmp(:) + FN_tmp(:) + FT_tmp(:)

! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
      DIST_CL = (DIST_LI**2 + RAD_L**2 - RAD_II**2)/&
         (2.d0*DIST_LI)

! total torque
      IF(DO_K) THEN
         CALL DES_CROSSPRDCT(TOW_tmp(:), DIST_CL*NORM, FT_TMP)
      ELSE
         TOW_tmp(1)  = DIST_CL*(NORM(1)*FT_TMP(2) - NORM(2)*FT_TMP(1))
      ENDIF

      RETURN
      END SUBROUTINE CFFCTOW

    END SUBROUTINE CALC_FORCE_DEM
