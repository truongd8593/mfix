!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: DISCRETELEMENT                                         C
!     Purpose: DES mod file                                               C
!                                                                         C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!     Reviewer:                                          Date:            C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
!     Common Block containing DEM conditions 
!     

      MODULE DISCRETELEMENT


      Use param
      Use param1
      
      parameter (NPARTICLES = 10000)
      parameter (MAXNEIGHBORS = 16)
      parameter (MAXQUADS = 50000)
      parameter (NMQD = 11)
      parameter (NDIM = 3)
      parameter (NWALLS = 6)

      INTEGER NEIGHBOURS(MAXNEIGHBORS,NPARTICLES)
      INTEGER PN(MAXNEIGHBORS,NPARTICLES), PV(MAXNEIGHBORS,NPARTICLES)
      INTEGER CALLED, PCOUNT, FACTOR, DES_NEIGHBOR_SEARCH
      INTEGER PARTICLES, MN, NQUAD, DIMN, ZONES, ZN1, ZN2, ZN3, ZN4
      INTEGER NEIGHBOR, WALLCONTACT, OUTOFBOX
      INTEGER LQUAD(NMQD,MAXQUADS), PINC(NPARTICLES)
      
      DOUBLE PRECISION DES_POS_OLD(NDIM,NPARTICLES)
      DOUBLE PRECISION DES_POS_NEW(NDIM,NPARTICLES)
      DOUBLE PRECISION DES_RADIUS(NPARTICLES), PR(NPARTICLES)
      DOUBLE PRECISION DES_VEL_OLD(NDIM,NPARTICLES)
      DOUBLE PRECISION DES_VEL_NEW(NDIM,NPARTICLES)
      DOUBLE PRECISION FC(NDIM,NPARTICLES)
      DOUBLE PRECISION FN(NDIM,NPARTICLES)
      DOUBLE PRECISION FT(NDIM,NPARTICLES) 
      DOUBLE PRECISION TOW(NDIM,NPARTICLES)
      DOUBLE PRECISION KN, KN_W
      DOUBLE PRECISION KT, KT_W 
      DOUBLE PRECISION ETA_DES_N, ETA_N_W
      DOUBLE PRECISION ETA_DES_T, ETA_T_W
      DOUBLE PRECISION MEW, MEW_W, E_RESTITUTION, DES_GAMMA, DES_F
      DOUBLE PRECISION DES_KE, DES_PE
      DOUBLE PRECISION OMEGA_OLD(NDIM,NPARTICLES)
      DOUBLE PRECISION OMEGA_NEW(NDIM,NPARTICLES)
      DOUBLE PRECISION WX1, EX2, BY1, TY2, SZ1, NZ2
      DOUBLE PRECISION DES_WALL_POS(NDIM,NWALLS), DES_WALL_VEL(NDIM,NWALLS)
      DOUBLE PRECISION RADIUS_EQ, WALL_NORMAL(NDIM,NWALLS)
      DOUBLE PRECISION PFN(NDIM,MAXNEIGHBORS,NPARTICLES)
      DOUBLE PRECISION PFT(NDIM,MAXNEIGHBORS,NPARTICLES)
      DOUBLE PRECISION CQUAD(NWALLS,MAXQUADS)
      DOUBLE PRECISION DTSOLID, RO_Sol(NPARTICLES)
      DOUBLE PRECISION N2CT, NBSCT, QUADCT, OCTCT
      DOUBLE PRECISION FNS1(NDIM), FTS1(NDIM)
      DOUBLE PRECISION PMASS(NPARTICLES), MOI(NPARTICLES), GRAV(NDIM), ROs
      DOUBLE PRECISION DES_THETA(NPARTICLES,NDIM)

      LOGICAL DISCRETE_ELEMENT 
      LOGICAL DES_CONTINUUM_COUPLED
      LOGICAL COUPLED_FLOW
      LOGICAL DO_NBS
      LOGICAL DO_QUADTREE
      LOGICAL DO_OCTREE
      LOGICAL DO_NSQUARE
      LOGICAL WALLFIXEDOVERLAP
      LOGICAL WALLDTSPLIT
      LOGICAL WALLREFLECT
      LOGICAL DES_PERIODIC_WALLS
      LOGICAL DES_PERIODIC_WALLS_X
      LOGICAL DES_PERIODIC_WALLS_Y
      LOGICAL DES_PERIODIC_WALLS_Z
      LOGICAL INLET_OUTLET
      LOGICAL TIME_ADJUST
      LOGICAL EQUIVALENT_RADIUS
      LOGICAL TSUJI_DRAG
      LOGICAL KUIPERS_DRAG


!  Square-well potential parameters
      DOUBLE PRECISION WELL_WIDTH(NPARTICLES)
      DOUBLE PRECISION WELL_DEPTH(NPARTICLES)

!  Well depth read in from mfix.dat
      DOUBLE PRECISION MASTER_WELL_DEPTH
      DOUBLE PRECISION MASTER_WALL_WELL_DEPTH

!  Ratio of well width to inner radius
      DOUBLE PRECISION RADIUS_RATIO
      DOUBLE PRECISION WALL_RADIUS_RATIO

!  Array of linked partners
      INTEGER LINKS(MAXNEIGHBORS, NPARTICLES)

!  Array of agglomerated partners
      INTEGER AGGS(MAXNEIGHBORS, NPARTICLES)

!  Does particle have at least one linked partner
      INTEGER IS_LINKED(NPARTICLES)

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
      INTEGER IS_AGGLOMERATED(NPARTICLES)

!  Number of particle that is followed by single particle log
      INTEGER FOCUS_PARTICLE

!  Number of search grids in the x-direction
      INTEGER SEARCH_GRIDS(3)

!  Matrix of particles in each search grid
      INTEGER PART_IN_GRID(100, 100, 100, 50)

!  Matrix location of particle 
      INTEGER PART_GRID(4,NPARTICLES)

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
      DOUBLE PRECISION RHODES_COHESION_LENGTH_SCALE
      DOUBLE PRECISION RHODES_COHESION_LENGTH_SCALE_WALL


!  Variables used in net force measurements
      LOGICAL RECORD_NET_FORCES
      INTEGER VERTICAL_NET_FORCE_BINS
      DOUBLE PRECISION NET_FORCE_COUNTER
      DOUBLE PRECISION VERTICAL_NET_FORCE_INCREMENT
      DOUBLE PRECISION NET_FORCE_TIME_INCREMENT
      DOUBLE PRECISION TIME_AT_PREVIOUS_NET_FORCE
      DOUBLE PRECISION NET_PART_TAN_FORCE(2,20)       !Contributions to net
      DOUBLE PRECISION NET_WALL_TAN_FORCE(2,20)       !  force in x- and y-
      DOUBLE PRECISION NET_PART_NORM_FORCE(2,20)      !  directions from up to
      DOUBLE PRECISION NET_WALL_NORM_FORCE(2,20)      !  20 different heights
      DOUBLE PRECISION NET_GRAVITY_FORCE(2,20)
      DOUBLE PRECISION NET_PART_COH_FORCE(2,20)
      DOUBLE PRECISION NET_WALL_COH_FORCE(2,20)
      DOUBLE PRECISION NET_FLUID_DRAG_FORCE(2,20)
      DOUBLE PRECISION NET_PRESSURE_FORCE(2,20)
      INTEGER APP_COH_DIST_INT
      INTEGER CAP_COH_DIST_INT
      INTEGER ESC_COH_DIST_INT


      END MODULE DISCRETELEMENT
