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
      USE functions

      use run, only: ENERGY_EQ
      use run, only: ANY_SPECIES_EQ

      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_GARG
      use particle_filter, only: DES_INTERP_DPVM
      use particle_filter, only: DES_INTERP_GAUSS
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_WEIGHT
! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: II, IJK
! the number of particles in the system
      INTEGER :: NPARTICLES
!-----------------------------------------------

      CALL INIT_ERR_MSG("DES_ALLOCATE_ARRAYS")

      NWALLS = merge(4,6,NO_K)

! Grab the larger of PARTICLES and MAX_PIS
      IF(MAX_PIS /= UNDEFINED_I .AND. PARTICLES /= UNDEFINED_I) THEN
         NPARTICLES = max(MAX_PIS, PARTICLES)
      ELSEIF(MAX_PIS /= UNDEFINED_I)THEN
         NPARTICLES = MAX_PIS
      ELSE
         NPARTICLES = PARTICLES
      ENDIF

! For parallel processing the array size required should be either
! specified by the user or could be determined from total particles
! with some factor.
      NPARTICLES = (NPARTICLES/numPEs)*PARTICLES_FACTOR+NWALLS
      IF(NPARTICLES < 1000) NPARTICLES = 1000

! max_pip adjusted to accomodate temporary variables used for walls
! and DES_MPI stuff
      MAX_PIP = NPARTICLES - 2*NWALLS - 3

      WRITE(ERR_MSG,1000) trim(iVal(NPARTICLES)), trim(iVal(MAX_PIP))
      CALL FLUSH_ERR_MSG(HEADER = .FALSE., FOOTER = .FALSE.)

 1000 FORMAT('DES Particle array size: ',A,/&
        'DES maximum particles per process: ',A)

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
      Allocate(  DES_POS_NEW (DIMN,NPARTICLES) )
      Allocate(  DES_VEL_NEW (DIMN,NPARTICLES) )
      Allocate(  OMEGA_NEW (DIMN,NPARTICLES) )

      IF (DO_OLD) THEN
         Allocate(  DES_POS_OLD (DIMN,NPARTICLES) )
         Allocate(  DES_VEL_OLD (DIMN,NPARTICLES) )
         Allocate(  DES_ACC_OLD (DIMN,NPARTICLES) )
         Allocate(  OMEGA_OLD (DIMN,NPARTICLES) )
         Allocate(  ROT_ACC_OLD (DIMN,NPARTICLES))
      ENDIF

! Particle positions at the last call neighbor search algorithm call
      Allocate(  PPOS (DIMN,NPARTICLES) )

! Total, normal and tangetial forces
      Allocate(  FC (DIMN,NPARTICLES) )

! Torque
      Allocate(  TOW (DIMN,NPARTICLES) )

      Allocate(  PARTICLE_WALL_COLLISIONS (NPARTICLES) )

! Temporary variables to store wall position, velocity and normal vector
      Allocate(  WALL_NORMAL  (NWALLS,DIMN) )

      OLD_PAIR_NUM = 0
      PAIR_NUM = 0
      PAIR_MAX = 1024
      Allocate(  PAIRS (2,PAIR_MAX) )
      Allocate(  PAIRS_OLD (2,PAIR_MAX) )
      Allocate(  PAIR_COLLIDES (PAIR_MAX) )
      Allocate(  FC_PAIR  (3,PAIR_MAX) )
      Allocate(  QQ_PAIR (PAIR_MAX) )
      Allocate(  PV_PAIR (PAIR_MAX) )
      Allocate(  PV_PAIR_OLD (PAIR_MAX) )
      Allocate(  PFT_PAIR (3,PAIR_MAX) )
      Allocate(  PFT_PAIR_OLD (3,PAIR_MAX) )
      Allocate(  PFN_PAIR (3,PAIR_MAX) )
      Allocate(  PFN_PAIR_OLD (3,PAIR_MAX) )
      Allocate(  TOW_PAIR (3,2,PAIR_MAX) )

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

      ALLOCATE(DRAG_AM(DIMENSION_3))
      ALLOCATE(DRAG_BM(DIMENSION_3, DIMN))
      ALLOCATE(F_gp(NPARTICLES ))
      F_gp(1:NPARTICLES)  = ZERO

! Explict drag force acting on a particle.
      Allocate(DRAG_FC (DIMN,NPARTICLES) )

! force due to gas-pressure gradient
      ALLOCATE(P_FORCE(DIMN, DIMENSION_3))

! Volume averaged solids volume in a computational fluid cell
      Allocate(  DES_U_s (DIMENSION_3, DES_MMAX) )
      Allocate(  DES_V_s (DIMENSION_3, DES_MMAX) )
      Allocate(  DES_W_s (DIMENSION_3, DES_MMAX) )

! Volume of nodes
      ALLOCATE(DES_VOL_NODE(DIMENSION_3))

      ALLOCATE(F_GDS(DIMENSION_3))
      ALLOCATE(VXF_GDS(DIMENSION_3))

      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_DPVM, DES_INTERP_GAUSS)
         IF(DO_K) THEN
            ALLOCATE(FILTER_CELL(27, NPARTICLES))
            ALLOCATE(FILTER_WEIGHT(27, NPARTICLES))
         ELSE
            ALLOCATE(FILTER_CELL(9, NPARTICLES))
            ALLOCATE(FILTER_WEIGHT(9, NPARTICLES))
         ENDIF
      CASE(DES_INTERP_GARG)
         ALLOCATE(DES_ROPS_NODE(DIMENSION_3, DES_MMAX))
         ALLOCATE(DES_VEL_NODE(DIMENSION_3, DIMN, DES_MMAX))
      END SELECT

! Variables for hybrid model
      IF (DES_CONTINUUM_HYBRID) THEN
         ALLOCATE(SDRAG_AM(DIMENSION_3,DIMENSION_M))
         ALLOCATE(SDRAG_BM(DIMENSION_3, DIMN,DIMENSION_M))

         ALLOCATE(F_SDS(DIMENSION_3,DIMENSION_M))
         ALLOCATE(VXF_SDS(DIMENSION_3,DIMENSION_M))
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



! ---------------------------------------------------------------->>>
! BEGIN COHESION
      IF(USE_COHESION) THEN
! Matrix location of particle  (should be allocated in case user wishes
! to invoke routines in /cohesion subdirectory
         Allocate(  PostCohesive (NPARTICLES) )
      ENDIF
! END COHESION
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

      CALL FINL_ERR_MSG

      RETURN
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


!``````````````````````````````````````````````````````````````````````!
! Subroutine: ADD_PAIR                                                 !
!                                                                      !
! Purpose: Adds a neighbor pair to the pairs array.                    !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE add_pair(ii,jj)
      USE discretelement
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ii,jj

      pair_num = pair_num +1

! Reallocate to double the size of the arrays if needed.
      IF(PAIR_NUM > PAIR_MAX) THEN
         PAIR_MAX = PAIR_MAX*2
         CALL PAIR_GROW
      ENDIF

      pairs(1,pair_num) = ii
      pairs(2,pair_num) = jj

      RETURN
      END SUBROUTINE add_pair


!``````````````````````````````````````````````````````````````````````!
! Subroutine: PAIR_GROW                                                !
!                                                                      !
! Purpose: Grow pair arrays to pair_max. Note that pair                !
! max should be increased before calling this routine. Also, no        !
! assumption to the previous array size is made as needed for restarts.!
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE PAIR_GROW

      USE discretelement

      IMPLICIT NONE

      LOGICAL, DIMENSION(:), ALLOCATABLE :: bool_tmp
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: int_tmp
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: real_tmp
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: real_tmp3
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: real_scalar_tmp

      INTEGER :: lSIZE1, lSIZE2, lSIZE3

      lSIZE2 = size(pairs,2)
      allocate(int_tmp(2,PAIR_MAX))
      int_tmp(:,1:lSIZE2) = pairs(:,1:lSIZE2)
      call move_alloc(int_tmp,pairs)

      lSIZE2 = size(pairs_old,2)
      allocate(int_tmp(2,PAIR_MAX))
      int_tmp(:,1:lSIZE2) = pairs_old(:,1:lSIZE2)
      call move_alloc(int_tmp,pairs_old)

      lSIZE3 = size(tow_pair,3)
      allocate(real_tmp3(3,2,PAIR_MAX))
      real_tmp3(:,:,1:lSIZE3) = tow_pair(:,:,1:lSIZE3)
      call move_alloc(real_tmp3,tow_pair)

      lSIZE2 = size(fc_pair,2)
      allocate(real_tmp(3,PAIR_MAX))
      real_tmp(:,1:lSIZE2) = fc_pair(:,1:lSIZE2)
      call move_alloc(real_tmp,fc_pair)

      lSIZE1 = size(pair_collides,1)
      allocate(bool_tmp(PAIR_MAX))
      bool_tmp(1:lSIZE1) = pair_collides(1:lSIZE1)
      call move_alloc(bool_tmp,pair_collides)

      lSIZE1 = size(qq_pair,1)
      allocate(real_scalar_tmp(PAIR_MAX))
      real_scalar_tmp(1:lSIZE1) = qq_pair(1:lSIZE1)
      call move_alloc(real_scalar_tmp,qq_pair)

      lSIZE1 = size(pv_pair,1)
      allocate(bool_tmp(PAIR_MAX))
      bool_tmp(1:lSIZE1) = pv_pair(1:lSIZE1)
      call move_alloc(bool_tmp,pv_pair)

      lSIZE1 = size(pv_pair_old,1)
      allocate(bool_tmp(PAIR_MAX))
      bool_tmp(1:lSIZE1) = pv_pair_old(1:lSIZE1)
      call move_alloc(bool_tmp,pv_pair_old)

      lSIZE2 = size(pft_pair_old,2)
      allocate(real_tmp(3,PAIR_MAX))
      real_tmp(:,1:lSIZE2) = pft_pair_old(:,1:lSIZE2)
      call move_alloc(real_tmp,pft_pair_old)

      lSIZE2 = size(pft_pair,2)
      allocate(real_tmp(3,PAIR_MAX))
      real_tmp(:,1:lSIZE2) = pft_pair(:,1:lSIZE2)
      call move_alloc(real_tmp,pft_pair)

      lSIZE2 = size(pfn_pair_old,2)
      allocate(real_tmp(3,PAIR_MAX))
      real_tmp(:,1:lSIZE2) = pfn_pair_old(:,1:lSIZE2)
      call move_alloc(real_tmp,pfn_pair_old)

      lSIZE2 = size(pfn_pair,2)
      allocate(real_tmp(3,PAIR_MAX))
      real_tmp(:,1:lSIZE2) = pfn_pair(:,1:lSIZE2)
      call move_alloc(real_tmp,pfn_pair)

      RETURN
      END SUBROUTINE PAIR_GROW
