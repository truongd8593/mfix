!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES                                            C
!
!  Purpose: DES - Calculate the new values of particle velocity,
!           position, angular velocity etc
!
!                                                                      C
!  Comments: Implements Eqns 1, 2, 3, 4 & 5  from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!  pradeep : changes for parallel processing
!          1. periodic boundaries might lie in different proc. so adjust
!             particle position for periodic removed
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE mfix_pic
      use geometry, only: DO_K, NO_K
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: L
      DOUBLE PRECISION :: DD(3), NEIGHBOR_SEARCH_DIST
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
      DOUBLE PRECISION :: OMEGA_MAG,OMEGA_UNIT(3),ROT_ANGLE
!-----------------------------------------------


! Adams-Bashforth defaults to Euler for the first time step.
      IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
         DO L =1, MAX_PIP
            IF(IS_NONEXISTENT(L)) CYCLE                       ! Only real particles
            IF(IS_ENTERING(L).or.IS_ENTERING_GHOST(L)) CYCLE  ! Only non-entering
            IF(IS_GHOST(L)) CYCLE                             ! Skip ghost particles
            DES_ACC_OLD(:,L) = FC(:,L)/PMASS(L) + GRAV(:)
            ROT_ACC_OLD(:,L) = TOW(:,L)
         ENDDO
      ENDIF

!$omp parallel do if(max_pip .ge. 10000) default(none)                    &
!$omp shared(MAX_PIP,INTG_EULER,INTG_ADAMS_BASHFORTH,fc,tow,              &
!$omp       omega_new,omega_old,pmass,grav,des_vel_new,des_pos_new,       &
!$omp       des_vel_old,des_pos_old,dtsolid,omoi,des_acc_old,rot_acc_old, &
!$omp       ppos,neighbor_search_rad_ratio,des_radius,DO_OLD, iGlobal_ID, &
!$omp       particle_orientation, orientation) &
!$omp private(l,dd,neighbor_search_dist,rot_angle,omega_mag,omega_unit)   &
!$omp reduction(.or.:do_nsearch) schedule (auto)

      DO L = 1, MAX_PIP
! only process particles that exist
         IF(IS_NONEXISTENT(L)) CYCLE
! skip ghost particles
         IF(IS_GHOST(L).or.IS_ENTERING_GHOST(L).or.IS_EXITING_GHOST(L)) CYCLE

! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(.NOT.IS_ENTERING(L) .AND. .NOT.IS_ENTERING_GHOST(L))THEN
            FC(:,L) = FC(:,L)/PMASS(L) + GRAV(:)
         ELSE
            FC(:,L) = ZERO
            TOW(:,L) = ZERO
         ENDIF

! Advance particle position, velocity
        IF (INTG_EULER) THEN
! first-order method
            DES_VEL_NEW(:,L) = DES_VEL_NEW(:,L) + FC(:,L)*DTSOLID
            DD(:) = DES_VEL_NEW(:,L)*DTSOLID
            DES_POS_NEW(:,L) = DES_POS_NEW(:,L) + DD(:)
            OMEGA_NEW(:,L)   = OMEGA_NEW(:,L) + TOW(:,L)*OMOI(L)*DTSOLID
         ELSEIF (INTG_ADAMS_BASHFORTH) THEN

! Second-order Adams-Bashforth/Trapezoidal scheme
            DES_VEL_NEW(:,L) = DES_VEL_OLD(:,L) + 0.5d0*&
               ( 3.d0*FC(:,L)-DES_ACC_OLD(:,L) )*DTSOLID
            OMEGA_NEW(:,L)   =  OMEGA_OLD(:,L) + 0.5d0*&
               ( 3.d0*TOW(:,L)*OMOI(L)-ROT_ACC_OLD(:,L) )*DTSOLID
            DD(:) = 0.5d0*( DES_VEL_OLD(:,L)+DES_VEL_NEW(:,L) )*DTSOLID
            DES_POS_NEW(:,L) = DES_POS_OLD(:,L) + DD(:)
            DES_ACC_OLD(:,L) = FC(:,L)
            ROT_ACC_OLD(:,L) = TOW(:,L)*OMOI(L)
         ENDIF

! Update particle orientation - Always first order
! When omega is non-zero, compute the rotation angle, and apply the
! Rodrigues' rotation formula

         IF(PARTICLE_ORIENTATION) THEN
            OMEGA_MAG = OMEGA_NEW(1,L)**2 +OMEGA_NEW(2,L)**2 + OMEGA_NEW(3,L)**2

            IF(OMEGA_MAG>ZERO) THEN
               OMEGA_MAG=DSQRT(OMEGA_MAG)
               OMEGA_UNIT(:) = OMEGA_NEW(:,L)/OMEGA_MAG
               ROT_ANGLE = OMEGA_MAG * DTSOLID

               ORIENTATION(:,L) = ORIENTATION(:,L)*DCOS(ROT_ANGLE) &
                                 + DES_CROSSPRDCT(OMEGA_UNIT,ORIENTATION(:,L))*DSIN(ROT_ANGLE) &
                                 + OMEGA_UNIT(:)*DOT_PRODUCT(OMEGA_UNIT,ORIENTATION(:,L))*(ONE-DCOS(ROT_ANGLE))
            ENDIF
         ENDIF

! Check if the particle has moved a distance greater than or equal to
! its radius during one solids time step. if so, call stop
         IF(dot_product(DD,DD).GE.DES_RADIUS(L)**2) THEN
            WRITE(*,1002) iGlobal_ID(L), sqrt(dot_product(DD,DD)), DES_RADIUS(L)
            IF (DO_OLD) WRITE(*,'(5X,A,3(ES17.9))') &
               'old particle pos = ', DES_POS_OLD(:,L)
            WRITE(*,'(5X,A,3(ES17.9))') &
               'new particle pos = ', DES_POS_NEW(:,L)
            WRITE(*,'(5X,A,3(ES17.9))')&
               'new particle vel = ', DES_VEL_NEW(:,L)
            WRITE(*,1003)
            STOP 1
         ENDIF

! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so,
! make sure that neighbor is called in des_time_march
         IF(.NOT.DO_NSEARCH) THEN
            DD(:) = DES_POS_NEW(:,L) - PPOS(:,L)
            NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO*&
               DES_RADIUS(L)
            IF(dot_product(DD,DD).GE.NEIGHBOR_SEARCH_DIST**2) DO_NSEARCH = .TRUE.
         ENDIF

! Reset total contact force and torque
         FC(:,L) = ZERO
         TOW(:,L) = ZERO

      ENDDO
!$omp end parallel do

      FIRST_PASS = .FALSE.


 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Particle ',I10, ' moved a distance ', ES17.9, &
         ' during a',/10X, 'single solids time step, which is ',&
         ' greater than',/10X,'its radius: ', ES17.9)
 1003 FORMAT(1X,70('*')/)

      RETURN
      END SUBROUTINE CFNEWVALUES
