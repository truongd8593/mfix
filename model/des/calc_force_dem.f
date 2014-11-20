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
      SUBROUTINE CALC_FORCE_DEM

!---------------------------------------------------------------------//
      USE calc_collision_wall
      USE constant, ONLY: Pi
      USE des_thermo
      USE des_thermo_cond
      USE discretelement
      USE geometry, ONLY: DO_K
      USE physprop, ONLY: K_s0
      USE run

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0
! particle no. indices
      INTEGER :: I, LL, CC
! loop counter and indices for neighbor information
      INTEGER :: NI, NLIM, N_NOCON
! the overlap occuring between particle-particle or particle-wall
! collision in the normal direction
      DOUBLE PRECISION :: OVERLAP_N
! square root of the overlap
      DOUBLE PRECISION :: SQRT_OVERLAP
! heat conducted
      DOUBLE PRECISION :: QQ
! distance vector between two particle centers or between a particle
! center and wall when the two surfaces are just at contact (i.e. no
! overlap)
      DOUBLE PRECISION :: R_LM,DIST_CI,DIST_CL
! the normal and tangential components of the translational relative
! velocity
      DOUBLE PRECISION :: V_REL_TRANS_NORM
! distance vector between two particle centers or between a particle
! center and wall at current and previous time steps
      DOUBLE PRECISION :: DIST(3)
! tangent to the plane of contact at current time step
      DOUBLE PRECISION :: V_REL_TANG(3)
! normal and tangential forces
      DOUBLE PRECISION :: FN(3)
      DOUBLE PRECISION :: FNS1(3), FNS2(3)
      DOUBLE PRECISION :: FTS1(3), FTS2(3)
! temporary storage of tangential DISPLACEMENT
      DOUBLE PRECISION :: PFT_TMP(3)
! temporary storage of torque
      DOUBLE PRECISION :: TOW_TMP(3)

! store solids phase index of particle (i.e. pijk(np,5))
      INTEGER :: PHASEI, PHASELL
! local values used spring constants and damping coefficients
      DOUBLE PRECISION :: ETAN_DES, ETAT_DES
      DOUBLE PRECISION :: KN_DES, KT_DES
! local values used for calculating cohesive forces
      DOUBLE PRECISION :: FORCE_COH, EQ_RADIUS, DistApart, &
                          magGravity
! set to T when a sliding contact occurs
      LOGICAL :: PARTICLE_SLIDE

      LOGICAL, PARAMETER :: report_excess_overlap = .FALSE.
!$      double precision omp_start, omp_end
!$      double precision omp_get_wtime

!-----------------------------------------------

! Initialize cohesive forces
      IF(USE_COHESION) PostCohesive(:) = ZERO

      CALL CALC_DEM_FORCE_WITH_WALL_STL

! Check particle LL neighbour contacts
!---------------------------------------------------------------------//

      FC_COLL(:,:) = 0
      FT_COLL(:,:) = 0

!$omp parallel do default(none) private(cc,ll,i,dist,r_lm,             &
!$omp    overlap_n,v_rel_tang,v_rel_trans_norm,sqrt_overlap,           &
!$omp    kn_des,kt_des,hert_kn,hert_kt,phasell,phasei,etan_des,        &
!$omp    etat_des,fns1,fns2,fts1,fts2,pft_tmp,fn,particle_slide,       &
!$omp    eq_radius,distapart,force_coh,qq_coll,k_s0)                   &
!$omp    shared(collisions,collision_num,des_pos_new,des_radius,       &
!$omp    des_coll_model_enum,kn,kt,pv_coll,pft_coll,pfn_coll,pijk,     &
!$omp    des_etan,des_etat,mew,fc_coll,use_cohesion,dist_coll,         &
!$omp    van_der_waals,vdw_outer_cutoff,vdw_inner_cutoff,norm_coll,    &
!$omp    hamaker_constant,asperities,surface_energy,ft_coll, pea, energy_eq)
      DO CC = 1, COLLISION_NUM
         LL = COLLISIONS(1,CC)
         I  = COLLISIONS(2,CC)

         IF(.NOT.PEA(LL,1)) CYCLE
         IF(.NOT.PEA(I, 1)) CYCLE

         R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
         DIST(:) = DES_POS_NEW(:,I) - DES_POS_NEW(:,LL)
         DIST_COLL(CC) = dot_product(DIST,DIST)

! Compute particle-particle VDW cohesive short-range forces
         IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
            IF(DIST_COLL(CC) < (R_LM+VDW_OUTER_CUTOFF)**2) THEN
               EQ_RADIUS = 2d0 * DES_RADIUS(LL)*DES_RADIUS(I) /        &
                    (DES_RADIUS(LL)+DES_RADIUS(I))
               IF(DIST_COLL(CC) > (VDW_INNER_CUTOFF+R_LM)**2) THEN
                  DistApart = (SQRT(DIST_COLL(CC))-R_LM)
                  FORCE_COH = HAMAKER_CONSTANT * EQ_RADIUS /           &
                     (12d0*DistApart**2) * (Asperities/(Asperities+    &
                     EQ_RADIUS) + ONE/(ONE+Asperities/DistApart)**2 )
               ELSE
                  FORCE_COH = 2d0 * PI * SURFACE_ENERGY * EQ_RADIUS *  &
                    (Asperities/(Asperities+EQ_RADIUS) + ONE/          &
                    (ONE+Asperities/VDW_INNER_CUTOFF)**2 )
               ENDIF
               FC_COLL(:,CC) = FC_COLL(:,CC) +                         &
                  DIST(:)*FORCE_COH/SQRT(DIST_COLL(CC))
            ENDIF
         ENDIF

         IF (ENERGY_EQ) THEN
            ! Calculate conduction and radiation for thermodynamic neighbors
            QQ_COLL(CC) = ZERO
            IF(K_s0(PIJK(LL,5)) > ZERO) THEN
               QQ_COLL(CC) = DES_CONDUCTION(LL, I, sqrt(DIST_COLL(CC)), PIJK(LL,5), PIJK(LL,4))
            ENDIF
         ENDIF

         IF(DIST_COLL(CC) > (R_LM + SMALL_NUMBER)**2) THEN
            PV_COLL(CC) = .false.
            PFT_COLL(:,CC) = 0.0
            PFN_COLL(:,CC) = 0.0
            CYCLE
         ENDIF

         IF(DIST_COLL(CC) == 0) THEN
            WRITE(*,8550) LL, I
            STOP "division by zero"
 8550 FORMAT('distance between particles is zero:',2(2x,I10))
         ENDIF
         DIST_COLL(cc) = SQRT(DIST_COLL(CC))
         NORM_COLL(:,CC)= DIST(:)/DIST_COLL(cc)

! Overlap calculation changed from history based to current position
         OVERLAP_N = R_LM-DIST_COLL(CC)

         IF (report_excess_overlap) call print_excess_overlap

! Calculate the components of translational relative velocity for a
! contacting particle pair and the tangent to the plane of contact
         CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, V_REL_TANG,            &
            NORM_COLL(:,CC), DIST_COLL(CC))

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
         FNS1(:) = -KN_DES * OVERLAP_N * NORM_COLL(:,CC)
         FNS2(:) = -ETAN_DES * V_REL_TRANS_NORM*NORM_COLL(:,CC)
         FN(:) = FNS1(:) + FNS2(:)

         call calc_tangential_displacement(pft_tmp(:),NORM_COLL(:,CC), &
            pfn_coll(:,cc),pft_coll(:,cc),overlap_n,                   &
            v_rel_trans_norm,v_rel_tang(:),PV_COLL(CC))
         PV_COLL(CC) = .true.

! Calculate the tangential contact force
         FTS1(:) = -KT_DES * PFT_TMP(:)
         FTS2(:) = -ETAT_DES * V_REL_TANG
         FT_COLL(:,CC) = FTS1(:) + FTS2(:)

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with another particle/wall
         PARTICLE_SLIDE = .FALSE.
         CALL CFSLIDE(V_REL_TANG(:), PARTICLE_SLIDE, MEW,              &
             FT_COLL(:,CC), FN(:))

! Calculate the total force FC of a collision
! total contact force ( FC_COLL may already include cohesive force)
         FC_COLL(:,CC) = FC_COLL(:,CC) + FN(:) + FT_COLL(:,CC)

! Save tangential displacement history with Coulomb's law correction
         IF (PARTICLE_SLIDE) THEN
! Since FT might be corrected during the call to cfslide, the tangential
! displacement history needs to be changed accordingly
            PFT_COLL(:,CC) = -( FT_COLL(:,CC) - FTS2(:) ) / KT_DES
         ELSE
            PFT_COLL(:,CC) = PFT_TMP(:)
         ENDIF

      ENDDO
!$omp end parallel do

      magGravity = SQRT(dot_product(GRAV,GRAV))

      DO CC = 1, COLLISION_NUM
         LL = COLLISIONS(1,CC)
         I  = COLLISIONS(2,CC)

         FC(:,LL) = FC(:,LL) + FC_COLL(:,CC)
         FC(:,I) = FC(:,I) - FC_COLL(:,CC)

! calculate the distance from the particles' centers to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
      DIST_CL = DIST_COLL(CC)/2.d0 + (DES_RADIUS(LL)**2 - &
         DES_RADIUS(I)**2)/(2.d0*DIST_COLL(CC))

      DIST_CI = DIST_COLL(CC) - DIST_CL
! total torque
      IF(DO_K) THEN
! for particle i flip the signs of both norm and ft, so we get the same
         CALL DES_CROSSPRDCT(TOW_tmp(:), NORM_COLL(:,CC), FT_COLL(:,CC))
         TOW(:,LL) = TOW(:,LL) + DIST_CL*TOW_tmp(:)
         TOW(:,I) = TOW(:,I) + DIST_CI*TOW_tmp(:)
      ELSE
         TOW_tmp(1) = (NORM_COLL(1,CC)*FT_COLL(2,CC) -                 &
            NORM_COLL(2,CC)*FT_COLL(1,CC))
         TOW(1,LL) = TOW(1,LL) + DIST_CL*TOW_tmp(1)
         TOW(1,I)  = TOW(1,I)  + DIST_CI*TOW_tmp(1)
      ENDIF

! just for post-processing mag. of cohesive forces on each particle
         IF(USE_COHESION)THEN
            PostCohesive(LL) = dot_product(FC_COLL(:,CC),FC_COLL(:,CC))
            if(magGravity> ZERO .AND. PEA(LL,1)) PostCohesive(LL) =    &
               SQRT(PostCohesive(LL)) / (PMASS(LL)*magGravity)
         ENDIF ! for cohesion model

         IF(ENERGY_EQ .and. QQ_COLL(CC).ne.ZERO) THEN
            Q_Source(LL) = Q_Source(LL) + QQ_COLL(CC)
            Q_Source(I) = Q_Source(I) - QQ_COLL(CC)
         ENDIF

      ENDDO

! Calculate drag
      CALL CALC_DRAG_DES

! Update the old values of particle position and velocity with the new
! values computed
      IF (INTG_ADAMS_BASHFORTH) CALL CFUPDATEOLD

      RETURN

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: print_excess_overlap                                    !
!                                                                      !
!  Purpose: Print overlap warning messages.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PRINT_EXCESS_OVERLAP

      use error_manager

      IF(OVERLAP_N > flag_overlap*DES_RADIUS(LL) .OR.                  &
         OVERLAP_N > flag_overlap*DES_RADIUS(I)) THEN

         WRITE(ERR_MSG,1000) trim(iVAL(LL)), trim(iVAL(I)), S_TIME,    &
            DES_RADIUS(LL), DES_RADIUS(I), OVERLAP_N

         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1000 FORMAT('WARNING: Excessive overplay detected between ',          &
         'particles ',A,' and ',/A,' at time ',g11.4,'.',/             &
         'RADII:  ',g11.4,' and ',g11.4,4x,'OVERLAP: ',g11.4)

      END SUBROUTINE PRINT_EXCESS_OVERLAP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_TANGENTIAL_DISPLACEMENT                            !
!                                                                      !
!  Purpose: Calculate the tangential displacement which is the         !
!  integration of tangential relative velocity with respect to         !
!  contact time.                                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_TANGENTIAL_DISPLACEMENT(PFT, NORM, NORM_OLD,     &
         SIGMAT_OLD, OVERLAP_NORM, RELVEL_TANG_NORM, RELVEL_TANG,      &
         ALREADY_COLLIDING)

      IMPLICIT NONE

! tangential displacement
      DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: PFT

! local variables for the accumulated tangential displacement that occurs
! for the particle-particle or particle-wall collision (current time
! step and previous time step)
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: SIGMAT_OLD
      LOGICAL, INTENT(IN) :: ALREADY_COLLIDING

! the overlap occuring between particle-particle or particle-wall
! collision in the tangential direction
      DOUBLE PRECISION :: OVERLAP_T(3)

      DOUBLE PRECISION, DIMENSION(3) :: SIGMAT

! time elapsed to travel the calculated normal overlap given the
! normal relative velocity
      DOUBLE PRECISION :: DTSOLID_TMP

! unit normal vector along the line of contact between contacting
! particles or particle-wall at previous time step
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: NORM, NORM_OLD
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: RELVEL_TANG
      DOUBLE PRECISION, INTENT(IN) :: OVERLAP_NORM, RELVEL_TANG_NORM

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

    END SUBROUTINE CALC_FORCE_DEM
