!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C 
!     Module name: DES_INIT_NAMELIST                                      C
!     Purpose: DES - initialize the des-namelist                          C
!                                                                         C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!     Reviewer:                                          Date:            C
!                                                                         C
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
      
              PARTICLES = 100 
              PARTICLES_FACTOR = 1.2D0
              MAXNEIGHBORS = 10
              MQUAD_FACTOR = 1.1D0
              IF(DIMN.EQ.3) THEN
                 NMQD = 11
              ELSE
                 NMQD = 7
              END IF
              PBP = 0.2*PARTICLES
              NFACTOR = 500
      
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
                  E_RESTITUTION = UNDEFINED
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
                  QLM = 1 ! Number of levels to go up in the tree to move particle to a new quad 
                  QLN = 1 ! Number of levels to go up in the tree to do neighbor search
                  DIMN = UNDEFINED_I
                  IFI = 0
                  NEIGHBOR_SEARCH_N = 1
                  NEIGHBOR_SEARCH_RAD_RATIO = 1000 
                  NEIGHBOR_SEARCH_DIST = UNDEFINED 
                  DES_NEIGHBOR_SEARCH = UNDEFINED_I 
                  USE_COHESION = .FALSE.

                  RETURN
                  END SUBROUTINE DES_INIT_NAMELIST
