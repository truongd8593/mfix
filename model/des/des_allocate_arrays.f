!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_ALLOCATE_ARRAYS                                     C
!  Purpose: Original allocte arrays subroutines for DES                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DES_ALLOCATE_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE constant
      USE discretelement
      Use indices
      Use geometry
      Use compar
      Use physprop
      Use des_bc
      Use pic_bc
      use funits
      use desgrid
      use desmpi
      USE mfix_pic
      Use des_thermo
      Use des_rxns
      USE cutcell

      use run, only: ENERGY_EQ
      use run, only: ANY_SPECIES_EQ

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K, IJK, M
! the number of particles in the system
      INTEGER :: NPARTICLES
!-----------------------------------------------
! Include statment functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

      CALL INIT_ERR_MSG("DES_ALLOCATE_ARRAYS")

      write(err_msg, 1000) Particles
 1000 FORMAT('Total number of particles = ', 2x, i15)
      CALL FLUSH_ERR_MSG(footer = .false.)

      NWALLS = merge(4,6,NO_K)
      MAXNEIGHBORS = MN + 1 + NWALLS

! defining the local variables nparticles
!----------------------------------------------------------------->>>

! +nwalls is included since calc_force_des temporarily uses the variables
! pos, vel, etc at elements beyond the array index given by particles.
! unclear why additional array space is needed via particles_factor
      NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS

! J.Musser : Dynamic Particle Info
! no real need to print this information out
      IF(MAX_PIS /= UNDEFINED_I .AND. MAX_PIS .GT. NPARTICLES) THEN
         !WRITE(ERR_MSG, 1001) MAX_PIS, NPARTICLES, MAX_PIS
         !CALL FLUSH_ERR_MSG(header = .false., footer = .false.)
         NPARTICLES = MAX_PIS
      ENDIF

 1001 FORMAT(/,'User supplied MAX_PIS (',I15,') > cummulative size ',&
         'of particle arrays, NPARTILCES ( ', I15 , ')', /, &
         'Therefore, setting NPARTICLES to ',I15)


! For parallel processing the array size required should be either
! specified by the user or could be determined from total particles
! with some factor
      NPARTICLES = (NPARTICLES/numPEs)

! minimum size for nparticles
      IF(NPARTICLES.LT.1000) NPARTICLES = 1000

      IF(NPARTICLES .LT. PIP) then
         write(err_msg, 1002) mype, Nparticles, pip
         CALL FLUSH_ERR_MSG(header = .false., footer = .false., Abort = .true.)
      endif
 1002 FORMAT(/,'Error 1002: For processor:', 2x, i5, /, &
      'Particle array size determined (', 2x, I15, &
      ') is less than number of particles (', 2x, i15, ')', /, &
      'Increase MAX_PIS or particles_factor in the input file')

! max_pip adjusted to accomodate temporary variables used for walls
! and DES_MPI stuff

      MAX_PIP = NPARTICLES - 2*NWALLS - 3

      WRITE(err_msg, 1003)  NPARTICLES, MAX_PIP

 1003 FORMAT('Particle array size on each proc = ',I15, /,&
      'Maximum possible physical particles (MAX_PIP) on each proc = ',I15,/, &
      'Note that this value of MAX_PIP is only ',&
      'relevant for a new run',/,'For restarts, max_pip ',&
      'will be set later on')

      CALL FLUSH_ERR_MSG(header = .false., footer = .false.)

! DES Allocatable arrays
!-----------------------------------------------
! Dynamic particle info including another index for parallel
! processing for ghost
      ALLOCATE( PEA (NPARTICLES, 4) )

      ALLOCATE (iglobal_id(nparticles))

! J.Musser: Allocate necessary arrays for DEM mass inlet/outlet BCs
      IF(DEM_BCMI /= 0 .OR. DEM_BCMO /=0) CALL ALLOCATE_DEM_MIO
! R.Garg: Allocate necessary arrays for PIC mass inlet/outlet BCs
      IF(PIC_BCMI /= 0 .OR. PIC_BCMO /=0) CALL ALLOCATE_PIC_MIO

! Particle attributes
! Radius, density, mass, moment of inertia
      Allocate(  DES_RADIUS (NPARTICLES) )
      Allocate(  RO_Sol (NPARTICLES) )
      Allocate(  PVOL (NPARTICLES) )
      Allocate(  PMASS (NPARTICLES) )
      Allocate(  OMOI (NPARTICLES) )

! Old and new particle positions, velocities (translational and
! rotational)
      Allocate(  DES_POS_OLD (DIMN,NPARTICLES) )
      Allocate(  DES_POS_NEW (DIMN,NPARTICLES) )
      Allocate(  DES_VEL_OLD (DIMN,NPARTICLES) )
      Allocate(  DES_VEL_NEW (DIMN,NPARTICLES) )
      Allocate(  DES_ACC_OLD (DIMN,NPARTICLES) )

      IF(DO_K) THEN
         Allocate(  OMEGA_OLD (DIMN,NPARTICLES) )
         Allocate(  OMEGA_NEW (DIMN,NPARTICLES) )
         ALLOCATE(  ROT_ACC_OLD (DIMN,NPARTICLES))
      ELSE
         Allocate(  OMEGA_OLD (1,NPARTICLES) )
         Allocate(  OMEGA_NEW (1,NPARTICLES) )
         ALLOCATE(  ROT_ACC_OLD (1,NPARTICLES))
      ENDIF

! Particle positions at the last call neighbor search algorithm call
      Allocate(  PPOS (DIMN,NPARTICLES) )

! Total, normal and tangetial forces
      Allocate(  FC (DIMN,NPARTICLES) )

! Torque
      IF(DO_K) THEN
         Allocate(  TOW (DIMN,NPARTICLES) )
         Allocate(  TOW_COLL (DIMN,COLLISION_MAX) )
      ELSE
         Allocate(  TOW (1,NPARTICLES) )
         Allocate(  TOW_COLL (1,COLLISION_MAX) )
      ENDIF

! Accumulated spring force
      Allocate(  PFT (NPARTICLES,0:MAXNEIGHBORS,DIMN) )
      Allocate(  PFT_WALL (NPARTICLES,6,DIMN) )

! Save the normal direction at previous time step
      Allocate(  PFN (NPARTICLES,MAXNEIGHBORS,DIMN) )
      Allocate(  PFN_WALL (NPARTICLES,6,DIMN) )

! Tracking variables for particle contact history
      Allocate(  PN (MAXNEIGHBORS, NPARTICLES) )
      Allocate(  PN_WALL (6, NPARTICLES) )
      Allocate(  PV (MAXNEIGHBORS, NPARTICLES) )
      Allocate(  PV_WALL (6, NPARTICLES) )


! Temporary variables to store wall position, velocity and normal vector
      Allocate(  WALL_NORMAL  (NWALLS,DIMN) )

! Neighbor search
      Allocate(  NEIGHBOURS (NPARTICLES, MAXNEIGHBORS) )

      OLD_COLLISION_NUM = 0
      COLLISION_NUM = 0
      COLLISION_MAX = 1024
      Allocate(  COLLISIONS (2,COLLISION_MAX) )
      Allocate(  COLLISIONS_OLD (2,COLLISION_MAX) )
      Allocate(  FC_COLL  (3,COLLISION_MAX) )
      Allocate(  PV_COLL (COLLISION_MAX) )
      Allocate(  PV_COLL_OLD (COLLISION_MAX) )
      Allocate(  PFT_COLL (3,COLLISION_MAX) )
      Allocate(  PFT_COLL_OLD (3,COLLISION_MAX) )
      Allocate(  PFN_COLL (3,COLLISION_MAX) )
      Allocate(  PFN_COLL_OLD (3,COLLISION_MAX) )

! Variable that stores the particle in cell information (ID) on the
! computational fluid grid defined by imax, jmax and kmax in mfix.dat
      ALLOCATE(PIC(DIMENSION_3))
      DO IJK=1,DIMENSION_3
        NULLIFY(pic(ijk)%p)
      ENDDO

! Particles in a computational fluid cell (for volume fraction)
      Allocate(  PINC (DIMENSION_3) )

! For each particle track its i,j,k location on computational fluid grid
! defined by imax, jmax and kmax in mfix.dat and phase no.
      Allocate(  PIJK (NPARTICLES,5) )

      IF(DES_INTERP_ON) THEN
         ALLOCATE(DRAG_AM(DIMENSION_3, DES_MMAX))
         ALLOCATE(DRAG_BM(DIMENSION_3, DIMN, DES_MMAX))
         ALLOCATE(VEL_FP(DIMN,NPARTICLES))
         ALLOCATE(F_gp(NPARTICLES ))
         F_gp(1:NPARTICLES)  = ZERO
      ENDIF

      IF(DES_INTERP_MEAN_FIELDS) THEN
         ALLOCATE(DES_ROPS_NODE(DIMENSION_3, DES_MMAX))
         ALLOCATE(DES_VEL_NODE(DIMENSION_3, DIMN, DES_MMAX))
      ENDIF

! force due to gas-pressure gradient
      ALLOCATE(P_FORCE(DIMENSION_3,DIMN))
! force due to gas-solids drag on a particle
      ALLOCATE(GD_FORCE(DIMN,NPARTICLES))

! Volume averaged solids volume in a computational fluid cell
      Allocate(  DES_U_s (DIMENSION_3, DES_MMAX) )
      Allocate(  DES_V_s (DIMENSION_3, DES_MMAX) )
      Allocate(  DES_W_s (DIMENSION_3, DES_MMAX) )

 ! Volume of nodes
       ALLOCATE(DES_VOL_NODE(DIMENSION_3))

! Variables for hybrid model
      IF (DES_CONTINUUM_HYBRID) THEN
         ALLOCATE(F_GDS(DIMENSION_3,DES_MMAX))
         ALLOCATE(F_SDS(DIMENSION_3,DIMENSION_M,DES_MMAX))
         ALLOCATE(VXF_GDS(DIMENSION_3,DES_MMAX))
         ALLOCATE(VXF_SDS(DIMENSION_3,DIMENSION_M,DES_MMAX))
         ALLOCATE(SD_FORCE(DIMN,NPARTICLES))
      ENDIF
! Bulk density in a computational fluid cell / for communication with
! MFIX continuum
      ALLOCATE( DES_ROP_S(DIMENSION_3, DES_MMAX) )
      ALLOCATE( DES_ROP_SO(DIMENSION_3, DES_MMAX) )

! MP-PIC related
      IF(MPPIC) THEN
         IF(.NOT.ALLOCATED(F_GP)) THEN
            ALLOCATE(F_GP(NPARTICLES ))
            F_GP(1:NPARTICLES)  = ZERO
         ENDIF

         IF(.NOT.ALLOCATED(VEL_FP)) ALLOCATE(VEL_FP(DIMN,NPARTICLES))

         Allocate(PS_FORCE_PIC(DIMENSION_3, DIMN))
         ALLOCATE(DES_STAT_WT(NPARTICLES))
         ALLOCATE(DES_VEL_MAX(DIMN))
         ALLOCATE(PS_GRAD(NPARTICLES, DIMN))
         ALLOCATE(AVGSOLVEL_P(DIMN, NPARTICLES))
         ALLOCATE(EPG_P(NPARTICLES))

         Allocate(PIC_U_S(DIMENSION_3, DES_MMAX))
         Allocate(PIC_V_S(DIMENSION_3, DES_MMAX))
         Allocate(PIC_W_S(DIMENSION_3, DES_MMAX))

         Allocate(PIC_P_s (DIMENSION_3, DES_MMAX) )
!         ALLOCATE(MPPIC_VPTAU(NPARTICLES, DIMN))
         PIC_U_s = zero
         PIC_V_s = zero
         PIC_W_s = zero
         PIC_P_s = zero
      ENDIF


! Granular temperature in a computational fluid cell
      Allocate(DES_THETA (DIMENSION_3, DES_MMAX) )

! Averaged velocity obtained by averaging over all the particles
      ALLOCATE(DES_VEL_AVG(DIMN) )

! Global Granular Energy
      ALLOCATE(GLOBAL_GRAN_ENERGY(DIMN) )
      ALLOCATE(GLOBAL_GRAN_TEMP(DIMN) )

! Somewhat free variable used to aid in manipulation
      ALLOCATE(MARK_PART(NPARTICLES))

! variable for bed height of solids phase M
      ALLOCATE(BED_HEIGHT(DES_MMAX))

! variable used to identify whether a particle had been put into a
! cluster
      IF (DES_CALC_CLUSTER) THEN
         ALLOCATE(InACluster(NPARTICLES))
      ENDIF

! ---------------------------------------------------------------->>>
! BEGIN COHESION
      IF(USE_COHESION) THEN
! Matrix location of particle  (should be allocated in case user wishes
! to invoke routines in /cohesion subdirectory
         Allocate(  FCohesive (DIMN,NPARTICLES) )
         Allocate(  PostCohesive (NPARTICLES) )
      ENDIF
! END COHESION
      IF(DES_CALC_CLUSTER) Allocate(  PostCluster (NPARTICLES) )
! ----------------------------------------------------------------<<<


! ---------------------------------------------------------------->>>
! BEGIN Thermodynamic Allocation
      IF(ENERGY_EQ)THEN
! Particle temperature
         Allocate( DES_T_s_OLD( NPARTICLES ) )
         Allocate( DES_T_s_NEW( NPARTICLES ) )
! Specific heat
         Allocate( DES_C_PS( NPARTICLES ) )
! Species mass fractions comprising a particle. This array may not be
! needed for all thermo problems.
         Allocate( DES_X_s( NPARTICLES, DIMENSION_N_S))
! Total rate of heat transfer to individual particles.
         Allocate( Q_Source( NPARTICLES ) )
! Average solids temperature in fluid cell
         Allocate(avgDES_T_s(DIMENSION_3) )

         Allocate(DES_ENERGY_SOURCE(DIMENSION_3) )

! Allocate the history variables for Adams-Bashforth integration
         IF (INTG_ADAMS_BASHFORTH) &
            Allocate( Q_Source0( NPARTICLES ) )
      ENDIF
! End Thermodynamic Allocation
! ----------------------------------------------------------------<<<


! ---------------------------------------------------------------->>>
! BEGIN Species Allocation
      IF(ANY_SPECIES_EQ)THEN
! Rate of solids phase production for each species
         Allocate( DES_R_sp( NPARTICLES, DIMENSION_N_s) )
! Rate of solids phase consumption for each species
         Allocate( DES_R_sc( NPARTICLES, DIMENSION_N_s) )


         Allocate( DES_R_gp( DIMENSION_3, DIMENSION_N_g ) )
         Allocate( DES_R_gc( DIMENSION_3, DIMENSION_N_g ) )
         Allocate( DES_SUM_R_g( DIMENSION_3 ) )
         Allocate( DES_R_PHASE( DIMENSION_3, DIMENSION_LM+DIMENSION_M-1 ) )
         Allocate( DES_HOR_g( DIMENSION_3 ) )


! Allocate the history variables for Adams-Bashforth integration
         IF (INTG_ADAMS_BASHFORTH) THEN
! Rate of chnage of particle mass
            Allocate( dMdt_OLD( NPARTICLES ) )
! Rate of chnage of particle mass percent species
            Allocate( dXdt_OLD( NPARTICLES, DIMENSION_N_s) )
         ENDIF

! Energy generation from reaction (cal/sec)
         Allocate( Qint( NPARTICLES ) )
      ENDIF
! End Species Allocation
! ----------------------------------------------------------------<<<

      IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(1X,A)')&
         '<---------- END DES_ALLOCATE_ARRAYS ----------'

      RETURN
      write(err_msg, '(A,/)') ''
      call flush_err_msg(header = .false.)
      CALL FINL_ERR_MSG
      END SUBROUTINE DES_ALLOCATE_ARRAYS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_DEM_MIO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE ALLOCATE_DEM_MIO

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE des_bc
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I     ! Loop counter for no. of DES_BCMI
!-----------------------------------------------

! Allocate/Initialize for inlets
      IF(DEM_BCMI /= 0)THEN

! Particle injection factor
         Allocate( PI_FACTOR (DEM_BCMI) )
! Particle injection count (injection number)
         Allocate( PI_COUNT (DEM_BCMI) )
! Particle injection time scale
         Allocate( DEM_MI_TIME (DEM_BCMI) )
! Array used for polydisperse inlets: stores the particle number
! distribution of an inlet scaled with numfrac_limit
         Allocate( DEM_BC_POLY_LAYOUT( DEM_BCMI, NUMFRAC_LIMIT ) )

         Allocate( DEM_MI(DEM_BCMI) )


! Initializiation
! Integer arrays
         PI_FACTOR(:) = -1
         PI_COUNT(:) = -1
         DEM_BC_POLY_LAYOUT(:,:) = -1
! Double precision arrays
         DEM_MI_TIME(:) = UNDEFINED

         allocate( DEM_BCMI_IJKSTART(DEM_BCMI) )
         allocate( DEM_BCMI_IJKEND(DEM_BCMI) )

         DEM_BCMI_IJKSTART = -1
         DEM_BCMI_IJKEND   = -1

      ENDIF  ! end if DEM_BCMI /= 0

! Boundary classification
!         Allocate( PARTICLE_PLCMNT (DES_BCMI) )
! Character precision arrays
!         PARTICLE_PLCMNT(:) = UNDEFINED_C


      IF(DEM_BCMO > 0)THEN
         allocate( DEM_BCMO_IJKSTART(DEM_BCMO) )
         allocate( DEM_BCMO_IJKEND(DEM_BCMO) )

         DEM_BCMO_IJKSTART = -1
         DEM_BCMO_IJKEND   = -1
      ENDIF


      RETURN
      END SUBROUTINE ALLOCATE_DEM_MIO


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_PIC_MIO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: R. Garg                                    Date: 11-Jun-14  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE ALLOCATE_PIC_MIO

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE pic_bc
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I     ! Loop counter for no. of DES_BCMI
!-----------------------------------------------

! Allocate/Initialize for inlets
      IF(PIC_BCMI /= 0)THEN


         allocate( PIC_BCMI_IJKSTART(PIC_BCMI) )
         allocate( PIC_BCMI_IJKEND  (PIC_BCMI) )
         allocate( PIC_BCMI_NORMDIR (PIC_BCMI,3) )

         ALLOCATE( PIC_BCMI_OFFSET  (PIC_BCMI,3))

         ALLOCATE( PIC_BCMI_INCL_CUTCELL(PIC_BCMI) )

         PIC_BCMI_IJKSTART = -1
         PIC_BCMI_IJKEND   = -1

      ENDIF  ! end if PIC_BCMI /= 0



      IF(PIC_BCMO > 0)THEN
         allocate( PIC_BCMO_IJKSTART(PIC_BCMO) )
         allocate( PIC_BCMO_IJKEND(PIC_BCMO) )

         PIC_BCMO_IJKSTART = -1
         PIC_BCMO_IJKEND   = -1
      ENDIF


      RETURN
      END SUBROUTINE ALLOCATE_PIC_MIO

      SUBROUTINE collision_add(ii,jj)
        USE discretelement
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ii,jj
        LOGICAL, DIMENSION(:), ALLOCATABLE :: bool_tmp
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: int_tmp
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: real_tmp

        collision_num = collision_num +1

        ! if we run out of space, reallocate to double the size of the arrays
        if (collision_num.gt.collision_max) then
           allocate(int_tmp(2,2*collision_max))
           int_tmp(:,1:collision_max) = collisions(:,1:collision_max)
           call move_alloc(int_tmp,collisions)

           allocate(int_tmp(2,2*collision_max))
           int_tmp(:,1:collision_max) = collisions_old(:,1:collision_max)
           call move_alloc(int_tmp,collisions_old)

           allocate(real_tmp(3,2*collision_max))
           real_tmp(:,1:collision_max) = fc_coll(:,1:collision_max)
           call move_alloc(real_tmp,fc_coll)

           allocate(real_tmp(3,2*collision_max))
           real_tmp(:,1:collision_max) = tow_coll(:,1:collision_max)
           call move_alloc(real_tmp,tow_coll)

           allocate(bool_tmp(2*collision_max))
           bool_tmp(1:collision_max) = pv_coll(1:collision_max)
           call move_alloc(bool_tmp,pv_coll)

           allocate(bool_tmp(2*collision_max))
           bool_tmp(1:collision_max) = pv_coll_old(1:collision_max)
           call move_alloc(bool_tmp,pv_coll_old)

           allocate(real_tmp(3,2*collision_max))
           real_tmp(:,1:collision_max) = pft_coll_old(:,1:collision_max)
           call move_alloc(real_tmp,pft_coll_old)

           allocate(real_tmp(3,2*collision_max))
           real_tmp(:,1:collision_max) = pft_coll(:,1:collision_max)
           call move_alloc(real_tmp,pft_coll)

           allocate(real_tmp(3,2*collision_max))
           real_tmp(:,1:collision_max) = pfn_coll_old(:,1:collision_max)
           call move_alloc(real_tmp,pfn_coll_old)

           allocate(real_tmp(3,2*collision_max))
           real_tmp(:,1:collision_max) = pfn_coll(:,1:collision_max)
           call move_alloc(real_tmp,pfn_coll)

           collision_max = 2*collision_max
        endif

        collisions(1,collision_num) = ii
        collisions(2,collision_num) = jj
        pv_coll(collision_num) = .false.
        pft_coll(:,collision_num) = 0.0
        pfn_coll(:,collision_num) = 0.0
      END SUBROUTINE collision_add
