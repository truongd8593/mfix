!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                       C
!   Module name: DISCRETELEMENT                                         C
!   Purpose: DES mod file                                               C
!                                                                       C
!                                                                       C
!   Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!   Reviewer:                                          Date:            C
!                                                                       C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!   
!   Common Block containing DEM conditions 
!   

      MODULE DISCRETELEMENT


      Use param
      Use param1

!
! DES Variables      
!
      INTEGER DIMN, MAXNEIGHBORS, MAXQUADS, NMQD, NWALLS, PBP
      DOUBLE PRECISION S_TIME, DES_SPX_TIME
      DOUBLE PRECISION DTSOLID, DTSOLID_FACTOR 
      DOUBLE PRECISION P_TIME, PTC
!
!   Particle properties 
      INTEGER PARTICLES
      DOUBLE PRECISION PARTICLES_FACTOR 
!
!   Particle-particle and Particle-wall contact parameters
!             Spring contants      
      DOUBLE PRECISION KN, KN_W  ! Normal
      DOUBLE PRECISION KT, KT_W  ! Tangential
!             Damping coeffients      
      DOUBLE PRECISION ETA_DES_N, ETA_N_W  ! Normal
      DOUBLE PRECISION ETA_DES_T, ETA_T_W  ! Tangential
!             Friction coefiicients and coeff of restitution
      DOUBLE PRECISION MEW, MEW_W, E_RESTITUTION 
!
!   Wall treatment      
      INTEGER WALLCONTACT
      DOUBLE PRECISION WX1, EX2, BY1, TY2, SZ1, NZ2
!             Wall vibration parameters
      DOUBLE PRECISION  DES_GAMMA, DES_F
!   
!   Neighbor search      
      INTEGER DES_NEIGHBOR_SEARCH, MN, NQUAD
      INTEGER QLM, QLN, INIT_QUAD_COUNT, INQC
      DOUBLE PRECISION RADIUS_EQ, NEIGHBOR_SEARCH_N
      DOUBLE PRECISION NEIGHBOR_SEARCH_RAD_RATIO, NEIGHBOR_SEARCH_DIST
      DOUBLE PRECISION N2CT, NBSCT, QUADCT, OCTCT, MQUAD_FACTOR
!
!   Kinetic and Potential energy of the system
      DOUBLE PRECISION DES_KE, DES_PE
!
!   Output file count
      INTEGER IFI
!
!    
! DES Logicals
!
!   DES - Continuum       
      LOGICAL DISCRETE_ELEMENT 
      LOGICAL DES_CONTINUUM_COUPLED
!
!   Slide check
      LOGICAL PARTICLE_SLIDE
!
!   Neighbor search     
      LOGICAL DO_QUADTREE
      LOGICAL DO_OCTREE
      LOGICAL DO_NSQUARE
      LOGICAL DO_NSEARCH
!
!   Particle treatment at the walls  
      LOGICAL WALLFIXEDOVERLAP
      LOGICAL WALLDTSPLIT
      LOGICAL WALLREFLECT
!
!   Periodic Wall BC
      LOGICAL DES_PERIODIC_WALLS
      LOGICAL DES_PERIODIC_WALLS_X
      LOGICAL DES_PERIODIC_WALLS_Y
      LOGICAL DES_PERIODIC_WALLS_Z
!
!   Inlet Outlet BC 
      LOGICAL INLET_OUTLET
      LOGICAL INLET_OUTLET_X
      LOGICAL INLET_OUTLET_Y
      LOGICAL INLET_OUTLET_Z
!
!   Drag      
      LOGICAL TSUJI_DRAG
!
!
!   Allocatable arrays
!
!   Particle attributes
!     Radius, density, mass, moment of inertia      
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_RADIUS ! (PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RO_Sol ! (PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PVOL !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PMASS ! (PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OMOI ! (PARTICLES)
!
!   Old and new particle positions, velocities (translational and
!                                                             rotational)      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_OLD ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_NEW ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_OLD ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_NEW ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_OLD ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_NEW ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PPOS ! (PARTICLES,DIMN)
!
!   Total, normal and tangetial forces      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FC ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FN ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FT ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FNS1 ! (DIMN)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FTS1 ! (DIMN)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GRAV ! (DIMN)
!
!   Torque      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TOW ! (PARTICLES,DIMN)
!
!   Accumulated spring forces      
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: PFN ! (PARTICLES,DIMN,MAXNEIGHBORS)
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: PFT ! (PARTICLES,DIMN,MAXNEIGHBORS)
!
!   Wall position, velocity and normal vector
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_WALL_POS ! (NWALLS,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_WALL_VEL ! (NWALLS,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WALL_NORMAL ! (NWALLS,DIMN)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PN ! (PARTICLES, MAXNEIGHBORS)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PV ! (PARTICLES, MAXNEIGHBORS)
!
!  Periodic walls
      INTEGER, DIMENSION(:), ALLOCATABLE :: WWALL ! (PBP)
      INTEGER, DIMENSION(:), ALLOCATABLE :: EWALL ! (PBP)
      INTEGER, DIMENSION(:), ALLOCATABLE :: BWALL ! (PBP)
      INTEGER, DIMENSION(:), ALLOCATABLE :: TWALL ! (PBP)
      INTEGER, DIMENSION(:), ALLOCATABLE :: SWALL ! (PBP)
      INTEGER, DIMENSION(:), ALLOCATABLE :: NWALL ! (PBP)
!
!   Particles in a computational cell (for volume fraction)
      INTEGER, DIMENSION(:), ALLOCATABLE :: PINC 
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PIJK ! (PARTCILES,5)=>I,J,K,IJK,M 

!
!   Volume averaged solids volume in a cell      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_U_s
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_V_s
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_W_s
!
!   Drag exerted by the gas o solids
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: SOLID_DRAG
!
!   Neighbor search
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NEIGHBOURS ! (PARTICLES, MAXNEIGHBORS)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: LQUAD ! (MAXQUADS, NMQD)
      INTEGER, DIMENSION(:), ALLOCATABLE :: PQUAD ! (PARTICLES) 
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CQUAD ! (NWALLS,MAXQUADS)
!
!   Neighbor distances
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PN_DIST ! (PARTICLES,MAXNEIGHBORS)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PN_RLM ! (PARTICLES,MAXNEIGHBORS)
!
!   Granular temperature
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_THETA ! (PARTICLES,MMAX)
!
!   Cell faces
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XE ! (IMAX3)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: YN ! (JMAX3)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZT ! (KMAX3)

!
!
!********************************************************************************
!
!  COHESION      
!
!  Square-well potential parameters
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WELL_WIDTH ! (PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WELL_DEPTH ! (PARTICLES)

!  Well depth read in from mfix.dat
      DOUBLE PRECISION MASTER_WELL_DEPTH
      DOUBLE PRECISION MASTER_WALL_WELL_DEPTH

!  Ratio of well width to inner radius
      DOUBLE PRECISION RADIUS_RATIO
      DOUBLE PRECISION WALL_RADIUS_RATIO

!  Array of linked partners
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: LINKS ! (PARTICLES, MAXNEIGHBORS)

!  Array of agglomerated partners
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: AGGS ! (PARTICLES, MAXNEIGHBORS)

!  Does particle have at least one linked partner
      INTEGER, DIMENSION(:), ALLOCATABLE :: IS_LINKED ! (PARTICLES)

!  Switch to turn cohesion on and off
      LOGICAL USE_COHESION      

!  Switch to turn square well on and off
      LOGICAL SQUARE_WELL

!  Switch to run debuggin on and off
      INTEGER COHESION_DEBUG
      INTEGER COHESION_DEBUG_START

!  Specific particle to target in debugging
      INTEGER COH_DEBUG_PARTICLE

!  Specific time step to target in debugging
      INTEGER COH_DEBUG_STEP

!  1x2 array to hold the last two interacing particles
      INTEGER LAST_COLLISION(2)

!  Parameters to control writing of animation file
      INTEGER ANIMATION_WRITE_INTERVAL
      INTEGER ANIMATION_WRITE_COUNTER
      DOUBLE PRECISION TIME_FACTOR
      INTEGER MAX_ANIMATOR_STEP
      INTEGER MIN_ANIMATOR_STEP

!  Parameters to control writing of single particle log and collision log
      INTEGER MAX_LOG_STEP
      INTEGER MIN_LOG_STEP

!  Does particle have at least one aggloerated partner
      INTEGER, DIMENSION(:), ALLOCATABLE :: IS_AGGLOMERATED ! (PARTICLES)

!  Number of particle that is followed by single particle log
      INTEGER FOCUS_PARTICLE

!  Number of search grids in the x-direction
      INTEGER SEARCH_GRIDS(3)

!  Matrix of particles in each search grid
      INTEGER PART_IN_GRID(0:110, 0:110, 0:110, 0:55)

!  Matrix location of particle 
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PART_GRID ! (PARTICLES,4)

!  Maximum number of particles in a search grid
      INTEGER MAX_PART_IN_GRID

!  Size of search grids
      DouBLE PRECISION SEARCH_GRID_SIZE(3)

!  Use new algorithm for inelastic collisions
      LOGICAL USE_COL_MW

!  Van der Waals constants
      LOGICAL VAN_DER_WAALS
      DOUBLE PRECISION HAMAKER_CONSTANT
      DOUBLE PRECISION VDW_INNER_CUTOFF       ! (in cm)
      DOUBLE PRECISION VDW_OUTER_CUTOFF
      DOUBLE PRECISION WALL_HAMAKER_CONSTANT
      DOUBLE PRECISION WALL_VDW_INNER_CUTOFF
      DOUBLE PRECISION WALL_VDW_OUTER_CUTOFF
      DOUBLE PRECISION SURFACE_ENERGY
      DOUBLE PRECISION WALL_SURFACE_ENERGY


!  Parameters to control Rhodes (2001) cohesion model (10/16/03)
      LOGICAL RHODES_COHESION
      DOUBLE PRECISION RHODES_COHESION_FACTOR
      DOUBLE PRECISION RHODES_COHESION_FACTOR_WALL
      DOuBLE PRECISION RHODES_COHESION_LENGTH_SCALE
      DOUBLE PRECISION RHODES_COHESION_LENGTH_SCALE_WALL


!  Variables used in net force measurements
      LOGICAL RECORD_NET_FORCES
      INTEGER VERTICAL_NET_FORCE_BINS
      DOUBLE PRECISION NET_FORCE_COUNTER
      DOUBLE PRECISION VERTICAL_NET_FORCE_INCREMENT
      DOUBLE PRECISION NET_FORCE_TIME_INCREMENT
      DOUBLE PRECISION TIME_AT_PREVIOUS_NET_FORCE
      DOUBLE PRECISION NET_PART_TAN_FORCE(20,2)       !Contributions to net
      DOUBLE PRECISION NET_WALL_TAN_FORCE(20,2)       !  force in x- and y-
      DOUBLE PRECISION NET_PART_NORM_FORCE(20,2)      !  directions from up to
      DOUBLE PRECISION NET_WALL_NORM_FORCE(20,2)      !  20 different heights
      DOUBLE PRECISION NET_GRAVITY_FORCE(20,2)
      DOUBLE PRECISION NET_PART_COH_FORCE(20,2)
      DOUBLE PRECISION NET_WALL_COH_FORCE(20,2)
      DOUBLE PRECISION NET_FLUID_DRAG_FORCE(20,2)
      DOUBLE PRECISION NET_PRESSURE_FORCE(20,2)
      INTEGER APP_COH_DIST_INT
      INTEGER CAP_COH_DIST_INT
      INTEGER ESC_COH_DIST_INT


      END MODULE DISCRETELEMENT
