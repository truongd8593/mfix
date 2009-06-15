!                                                                       C
!   Module name: DISCRETELEMENT                                         C
!>   Purpose: DES mod file 
!                                                                       C
!                                                                       C
!   Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!   Reviewer: Jin Sun and Rahul Garg                   Date: 01-Aug-07  C
!   Comments: Added declaration of interpolation related data           C
!                                                                       C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!   
!     Common Block containing DEM conditions 
!     

      MODULE DISCRETELEMENT


      USE param
      USE param1

!===========START of Interpolation related data======================
!     the coefficient add to gas momentum A matrix  at cell corners
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE ::drag_am 

!     the coefficient add to gas momentum B matrix  at cell corners
      DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE ::drag_bm 

!     fluid velocity at particle position
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::vel_fp 
      
      DOUBLE PRECISION, DIMENSION(:,:,:),POINTER :: weightp    
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: f_gp 
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: wtderivp, wtbar
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: sstencil
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: gstencil, vstencil,pgradstencil
      
      LOGICAL:: intx_per, inty_per, intz_per
      CHARACTER(LEN=7):: scheme, interp_scheme
      INTEGER:: order, ob2l, ob2r
      
      TYPE iap1
      INTEGER, DIMENSION(:), POINTER:: p
      END TYPE iap1
!     id's of particles in a cell 
      TYPE(iap1), DIMENSION(:,:,:), ALLOCATABLE:: pic
!===========END of Interpolation related data======================

     
!     DES Variables      
     
      INTEGER, PARAMETER :: DES_EXTRA_UNIT = 2000, DES_VOLFRAC_UNIT = 2001
      DOUBLE PRECISION DES_SPX_TIME, DES_RES_TIME
      DOUBLE PRECISION DTSOLID_FACTOR 

      DOUBLE PRECISION DTSOLID, S_TIME 

!     Print DES Data
      LOGICAL PRINT_DES_DATA 

!     usr specified time interal that controls frequency of writing DEM 
!     output when doing pure granular flow simulation
      DOUBLE PRECISION P_TIME

!     If true, then DEM output data is written in tecplot format
      LOGICAL :: DEM_OUTPUT_DATA_TECPLOT 

      LOGICAL :: DEBUG_DES

!     Output file count
      INTEGER IFI

!     if gener_part_config is true, then the particle_input.dat file
!     does not need to be supplied nor does the total number of
!     particles as these are determined based on the specified volume
!     fraction (vol_frac) in the specified domain (des_eps_xyzstart)     
      LOGICAL :: GENER_PART_CONFIG
      DOUBLE PRECISION ::  VOL_FRAC(DIM_M), DES_EPS_XSTART, &
                           DES_EPS_YSTART, DES_EPS_ZSTART
!     the number of particles that belong to solid phase M according
!     to the D_p0 specified in the mfix.dat 
      INTEGER PART_MPHASE(DIM_M)

!     Total number of particles in simulation 
      INTEGER PARTICLES

!     Particle-particle and Particle-wall contact parameters
!     Spring contants      
      DOUBLE PRECISION KN, KN_W ! Normal
      DOUBLE PRECISION KT, KT_W, KT_FAC, KT_W_FAC ! Tangential factors = KT/KN and KT_w/KN_w, resp.
!     Damping coeffients      
      DOUBLE PRECISION ETA_DES_N, ETA_N_W ! Normal
      DOUBLE PRECISION ETA_DES_T, ETA_T_W ! Tangential
!     Tangential damping factors, eta_t = eta_t_factor * eta_N
      DOUBLE PRECISION DES_ETAT_FAC, DES_ETAT_W_FAC
!     Damping coeffients in array form 
      DOUBLE PRECISION , DIMENSION(:,:), ALLOCATABLE :: DES_ETAN, DES_ETAT         !(MMAX, MMAX)
      DOUBLE PRECISION , DIMENSION(:), ALLOCATABLE :: DES_ETAN_WALL, DES_ETAT_WALL !(MMAX)
!     Friction coeficients
      DOUBLE PRECISION MEW, MEW_W
!     coeff of restituion input in one D array, solid solid
!     Tangential rest. coef. are not used in the code and thus are removed (sof) DEC-04-2008
      DOUBLE PRECISION DES_EN_INPUT(DIM_M+DIM_M*(DIM_M-1)/2) !DES_ET_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)
!     coeff of restituion input in one D array, solid wall 
      DOUBLE PRECISION  DES_EN_WALL_INPUT(DIM_M) !  DES_ET_WALL_INPUT(DIM_M)
!     actual coeff of rest.'s rearranged 
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  REAL_EN
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  REAL_EN_WALL
     
!     Neighbor search      
      INTEGER DES_NEIGHBOR_SEARCH, NQUAD, NEIGH_MAX
      INTEGER QLM, QLN, INIT_QUAD_COUNT, INQC
      INTEGER MAXQUADS, NMQD 
      DOUBLE PRECISION RADIUS_EQ, NEIGHBOR_SEARCH_N
      DOUBLE PRECISION NEIGHBOR_SEARCH_RAD_RATIO
      DOUBLE PRECISION N2CT, QUADCT, OCTCT, MQUAD_FACTOR
      LOGICAL DO_NSEARCH

!     factor for sum of radii in des_grid_based_neighbor_search
      DOUBLE PRECISION FACTOR_RLM

      INTEGER DIMN, NWALLS
      INTEGER MN, MAXNEIGHBORS
      DOUBLE PRECISION PARTICLES_FACTOR 
      INTEGER NFACTOR
      DOUBLE PRECISION lid_vel

!     DES - Continuum       
      LOGICAL DISCRETE_ELEMENT 
      LOGICAL DES_CONTINUUM_COUPLED

!     Switch to decide whether to call drag_gs or to call des_drag_gs via
!     drag_fgs to calculate the gas-solids drag coefficient.  if false
!     then drag_gs is used, otherwise des_drag_gs via drag_fgs is used
      LOGICAL DES_INTERP_ON

!     Drag      
      LOGICAL TSUJI_DRAG

!     run time logic.  if calc_fc is true, then the contact forces (FC) are 
!     updated to include gas-solids drag and gas pressure in the call to 
!     drag_fgs.  calc_fc does not play a role in pure granular flow simulations
      LOGICAL CALC_FC

!     run time logic.  if callfromdes is true, then the pertinent mean fields
!     (in this case ROP_S and F_GS) are not computed/updated in the call to
!     drag_fgs.  it is done to speed up the simulation. callfromdes does
!     not play a role in pure granular flow simulations and is only relevant
!     when des_interp_on is set to T
      LOGICAL CALLFROMDES
   
!     Wall vibration parameters
      DOUBLE PRECISION  DES_GAMMA, DES_F

!     Particle treatment at the walls  
      LOGICAL WALLDTSPLIT
      LOGICAL WALLREFLECT
 
!     variables for cell_near_wall
      LOGICAL NON_RECT_BC   
      INTEGER NPC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_near_w

!     Periodic Wall BC
      LOGICAL DES_PERIODIC_WALLS
      LOGICAL DES_PERIODIC_WALLS_X
      LOGICAL DES_PERIODIC_WALLS_Y
      LOGICAL DES_PERIODIC_WALLS_Z
     
!     Inlet Outlet BC 
      LOGICAL INLET_OUTLET
      LOGICAL INLET_OUTLET_X
      LOGICAL INLET_OUTLET_Y
      LOGICAL INLET_OUTLET_Z

!     Position of domain boundaries
      DOUBLE PRECISION WX1, EX2, BY1, TY2, SZ1, NZ2

!     run time logic. is set to T when a sliding contact occurs
      LOGICAL PARTICLE_SLIDE

!     Kinetic and Potential energy of the system
      DOUBLE PRECISION DES_KE, DES_PE

!     Global Granular Energy
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  GLOBAL_GRAN_ENERGY,GLOBAL_GRAN_TEMP
     
!     Constant input pressure gradient 
      DOUBLE PRECISION  pgrad(3)

!     Intial particle velocity distribution's mean and Standard Deviation
      DOUBLE PRECISION pvel_mean, PVEL_StDev

      DOUBLE PRECISION AVG_RAD, RMS_RAD
      DOUBLE PRECISION :: MIN_RADIUS, MAX_RADIUS
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE::   AVE_VEL_X, AVE_VEL_Y,  AVE_VEL_Z

      INTEGER, ALLOCATABLE, DIMENSION(:) :: MARK_PART
      DOUBLE PRECISION OVERLAP_MAX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bed_height   


!     Allocatable arrays
     
!     Particle attributes
!     Radius, density, mass, moment of inertia      
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_RADIUS ! (PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RO_Sol ! (PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PVOL !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PMASS ! (PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OMOI ! (PARTICLES)
     
!     Old and new particle positions, velocities (translational and
!     rotational)      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_OLD ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_NEW ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_OLD ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_NEW ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_OLD ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_NEW ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PPOS ! (PARTICLES,DIMN)
     
!     Total, normal and tangetial forces      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FC ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FN ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FT ! (PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FNS2 ! (DIMN)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FTS2 ! (DIMN)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FNS1 ! (DIMN)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FTS1 ! (DIMN)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GRAV ! (DIMN)
     
!     Torque      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TOW ! (PARTICLES,DIMN)
     
!     Accumulated spring forces      
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: PFN ! (PARTICLES,DIMN,MAXNEIGHBORS)
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: PFT ! (PARTICLES,DIMN,MAXNEIGHBORS)
     
!     Wall position, velocity and normal vector
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_WALL_POS ! (NWALLS,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_WALL_VEL ! (NWALLS,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WALL_NORMAL ! (NWALLS,DIMN)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PN ! (PARTICLES, MAXNEIGHBORS)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PV ! (PARTICLES, MAXNEIGHBORS)
     
!     Particles in a computational cell (for volume fraction)
      INTEGER, DIMENSION(:), ALLOCATABLE :: PINC 
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PIJK ! (PARTCILES,5)=>I,J,K,IJK,M 

!     Volume averaged solids velocity in a cell      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_U_s
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_V_s
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_W_s

!     Averaged velocity obtained by avraging over all the particles
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_VEL_AVG
     
!     Drag exerted by the gas on solids
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: SOLID_DRAG
     
!     Neighbor search
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NEIGHBOURS ! (PARTICLES, MAXNEIGHBORS)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: LQUAD ! (MAXQUADS, NMQD)
      INTEGER, DIMENSION(:), ALLOCATABLE :: PQUAD ! (PARTICLES) 
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CQUAD ! (NWALLS,MAXQUADS)
     
!     Neighbor distances
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PN_DIST ! (PARTICLES,MAXNEIGHBORS)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PN_RLM ! (PARTICLES,MAXNEIGHBORS)
   
!     Granular temperature
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_THETA ! (PARTICLES,MMAX)
     
!     Cell faces
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XE ! (IMAX3)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: YN ! (JMAX3)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZT ! (KMAX3)

     
!********************************************************************************
!     COHESION      
     
!     Square-well potential parameters
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WELL_WIDTH ! (PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WELL_DEPTH ! (PARTICLES)

!     Well depth read in from mfix.dat
      DOUBLE PRECISION MASTER_WELL_DEPTH
      DOUBLE PRECISION MASTER_WALL_WELL_DEPTH

!     Ratio of well width to inner radius
      DOUBLE PRECISION RADIUS_RATIO
      DOUBLE PRECISION WALL_RADIUS_RATIO

!     Array of linked partners
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: LINKS ! (PARTICLES, MAXNEIGHBORS)

!     Array of agglomerated partners
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: AGGS ! (PARTICLES, MAXNEIGHBORS)

!     Does particle have at least one linked partner
      INTEGER, DIMENSION(:), ALLOCATABLE :: IS_LINKED ! (PARTICLES)

!     Switch to turn cohesion on and off
      LOGICAL USE_COHESION      

!     Switch to turn square well on and off
      LOGICAL SQUARE_WELL

!     Switch to run debuggin on and off
      INTEGER COHESION_DEBUG
      INTEGER COHESION_DEBUG_START

!     Specific particle to target in debugging
      INTEGER COH_DEBUG_PARTICLE

!     Specific time step to target in debugging
      INTEGER COH_DEBUG_STEP

!     1x2 array to hold the last two interacing particles
      INTEGER LAST_COLLISION(2)

!     Parameters to control writing of animation file
      INTEGER ANIMATION_WRITE_INTERVAL
      INTEGER ANIMATION_WRITE_COUNTER
      DOUBLE PRECISION TIME_FACTOR
      INTEGER MAX_ANIMATOR_STEP
      INTEGER MIN_ANIMATOR_STEP

!     Parameters to control writing of single particle log and collision log
      INTEGER MAX_LOG_STEP
      INTEGER MIN_LOG_STEP

!     Does particle have at least one aggloerated partner
      INTEGER, DIMENSION(:), ALLOCATABLE :: IS_AGGLOMERATED ! (PARTICLES)

!     Number of particle that is followed by single particle log
      INTEGER FOCUS_PARTICLE

!     Number of search grids in the x-direction
      INTEGER SEARCH_GRIDS(3)

!     Matrix of particles in each search grid
      INTEGER PART_IN_GRID(0:110, 0:110, 0:110, 0:55)

!     Matrix location of particle 
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PART_GRID ! (PARTICLES,4)

!     Maximum number of particles in a search grid
      INTEGER MAX_PART_IN_GRID

!     Size of search grids
      DOUBLE PRECISION SEARCH_GRID_SIZE(3)

!     Use new algorithm for inelastic collisions
      LOGICAL USE_COL_MW

!     Van der Waals constants
      LOGICAL VAN_DER_WAALS
      DOUBLE PRECISION HAMAKER_CONSTANT
      DOUBLE PRECISION VDW_INNER_CUTOFF ! (in cm)
      DOUBLE PRECISION VDW_OUTER_CUTOFF
      DOUBLE PRECISION WALL_HAMAKER_CONSTANT
      DOUBLE PRECISION WALL_VDW_INNER_CUTOFF
      DOUBLE PRECISION WALL_VDW_OUTER_CUTOFF
      DOUBLE PRECISION SURFACE_ENERGY
      DOUBLE PRECISION WALL_SURFACE_ENERGY

!     Parameters to control Rhodes (2001) cohesion model (10/16/03)
      LOGICAL RHODES_COHESION
      DOUBLE PRECISION RHODES_COHESION_FACTOR
      DOUBLE PRECISION RHODES_COHESION_FACTOR_WALL
      DOUBLE PRECISION RHODES_COHESION_LENGTH_SCALE
      DOUBLE PRECISION RHODES_COHESION_LENGTH_SCALE_WALL

!     Variables used in net force measurements
      LOGICAL RECORD_NET_FORCES
      INTEGER VERTICAL_NET_FORCE_BINS
      DOUBLE PRECISION NET_FORCE_COUNTER
      DOUBLE PRECISION VERTICAL_NET_FORCE_INCREMENT
      DOUBLE PRECISION NET_FORCE_TIME_INCREMENT
      DOUBLE PRECISION TIME_AT_PREVIOUS_NET_FORCE
      DOUBLE PRECISION NET_PART_TAN_FORCE(20,2) !Contributions to net
      DOUBLE PRECISION NET_WALL_TAN_FORCE(20,2) !  force in x- and y-
      DOUBLE PRECISION NET_PART_NORM_FORCE(20,2) !  directions from up to
      DOUBLE PRECISION NET_WALL_NORM_FORCE(20,2) !  20 different heights
      DOUBLE PRECISION NET_GRAVITY_FORCE(20,2)
      DOUBLE PRECISION NET_PART_COH_FORCE(20,2)
      DOUBLE PRECISION NET_WALL_COH_FORCE(20,2)
      DOUBLE PRECISION NET_FLUID_DRAG_FORCE(20,2)
      DOUBLE PRECISION NET_PRESSURE_FORCE(20,2)
      INTEGER APP_COH_DIST_INT
      INTEGER CAP_COH_DIST_INT
      INTEGER ESC_COH_DIST_INT

!     END COHESION 
!********************************************************************************
     
      END MODULE DISCRETELEMENT
