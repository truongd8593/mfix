!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C 
!     Module name: DES_INIT_NAMELIST                                      C
!     Purpose: DES - initialize the des-namelist                          C
!                                                                         C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!     Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!     Comments: Added some interpolation based inputs                     C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
      SUBROUTINE DES_INIT_NAMELIST 

      USE param1
      USE discretelement
      
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!     loop counters
      INTEGER :: LC, LCM, M, N, K, KKK
!     
!     Coefficient of restitution (old symbol)
      DOUBLE PRECISION :: E
!-----------------------------------------------
!     
!     
      
      INCLUDE 'des/desnamelist.inc'
      
      PARTICLES = UNDEFINED_I
      PARTICLES_FACTOR = 1.2D0
      MAXNEIGHBORS = 10
      MQUAD_FACTOR = 1.1D0
      PBP = 0.2*PARTICLES
      NFACTOR = 10
      
      DISCRETE_ELEMENT = .FALSE.
      DO_QUADTREE = .FALSE.
      DO_OCTREE = .FALSE.
      DO_NSQUARE = .FALSE.
      DO_NSEARCH = .FALSE.
      WALLFIXEDOVERLAP = .FALSE.
      WALLDTSPLIT = .FALSE.
      WALLREFLECT = .FALSE.
      DES_PERIODIC_WALLS = .FALSE.
      DES_PERIODIC_WALLS_X = .FALSE.
      DES_PERIODIC_WALLS_Y = .FALSE.
      DES_PERIODIC_WALLS_Z = .FALSE.
      INLET_OUTLET = .FALSE.
      INLET_OUTLET_X = .FALSE.
      INLET_OUTLET_Y = .FALSE.
      INLET_OUTLET_Z = .FALSE.
      DES_CONTINUUM_COUPLED = .FALSE.
      TSUJI_DRAG = .FALSE.
      PARTICLE_SLIDE = .FALSE.
      DEBUG_DES = .FALSE.
      NON_RECT_BC=.FALSE.

      KN = UNDEFINED
      KT = UNDEFINED
      KN_W = UNDEFINED
      KT_W = UNDEFINED
      ETA_DES_N = UNDEFINED
      ETA_DES_T = UNDEFINED
      ETA_N_W = UNDEFINED
      ETA_T_W = UNDEFINED
      MEW = UNDEFINED
      MEW_W = UNDEFINED
      DES_GAMMA = ZERO
      DES_F = ZERO
      DES_KE = UNDEFINED
      DES_PE = UNDEFINED

      S_TIME = UNDEFINED
      DES_SPX_TIME = UNDEFINED
      DTSOLID = UNDEFINED
      DTSOLID_FACTOR = 0.1D0
      P_TIME = UNDEFINED
      PTC = UNDEFINED

      N2CT = UNDEFINED
      NBSCT = UNDEFINED
      QUADCT = UNDEFINED
      OCTCT = UNDEFINED

      WX1 = UNDEFINED
      EX2 = UNDEFINED
      BY1 = UNDEFINED
      TY2 = UNDEFINED
      SZ1 = UNDEFINED
      NZ2 = UNDEFINED
      RADIUS_EQ = UNDEFINED
      NQUAD = UNDEFINED_I
      INIT_QUAD_COUNT = UNDEFINED_I 
      INQC = UNDEFINED_I 
      QLM = 1                   ! Number of levels to go up in the tree to move particle to a new quad 
      QLN = 1                   ! Number of levels to go up in the tree to do neighbor search
      DIMN = UNDEFINED_I
      IFI = 0
      NEIGHBOR_SEARCH_N = 25
      NEIGHBOR_SEARCH_RAD_RATIO = 1.0D0
      NEIGHBOR_SEARCH_DIST = UNDEFINED 
      DES_NEIGHBOR_SEARCH = UNDEFINED_I 
      USE_COHESION = .FALSE.
      FACTOR_RLM = 1.2 
      DES_INTERP_ON = .FALSE.
      intx_per = .false.
      inty_per = .false.
      intz_per = .false.
      NPC = 1
      pgrad(:) = zero 
      pvel_mean = zero 
      PVEL_StDev = zero 

      DES_EN_INPUT(:) = UNDEFINED
      DES_ET_INPUT(:) = UNDEFINED
      
      DES_EN_WALL_INPUT(:) = UNDEFINED
      DES_ET_WALL_INPUT(:) = UNDEFINED
      
      DEM_OUTPUT_DATA_TECPLOT = .FALSE.
      LID_VEL = 0.0d0

      GENER_PART_CONFIG = .FALSE.
      VOL_FRAC(:) = UNDEFINED
      DES_EPS_XSTART = UNDEFINED
      DES_EPS_YSTART = UNDEFINED 
      DES_EPS_ZSTART = UNDEFINED 	 
      RETURN
      END SUBROUTINE DES_INIT_NAMELIST
