!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C 
!     Module name: DES_INIT_NAMELIST                                      C
!>     Purpose: DES - initialize the des-namelist                          
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
      USE des_bc
      
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------

      INCLUDE 'des/desnamelist.inc'

      PARTICLES = UNDEFINED_I
      FACTOR_RLM = 1.2 
      MN = 10
      NFACTOR = 10

      NON_RECT_BC=.FALSE.
      PARTICLES_FACTOR = 1.2D0

      DISCRETE_ELEMENT = .FALSE.
      DES_CONTINUUM_COUPLED = .FALSE.
      DES_INTERP_ON = .FALSE.
      TSUJI_DRAG = .FALSE.
      INT_METHOD = 0

      PARTICLE_SLIDE = .FALSE.
      DO_NSEARCH = .FALSE.

      DES_NEIGHBOR_SEARCH = 1
      NEIGHBOR_SEARCH_N = 25
      NEIGHBOR_SEARCH_RAD_RATIO = 1.0D0

      USE_COHESION = .FALSE.
      RADIUS_EQ = UNDEFINED     

      WALLDTSPLIT = .FALSE.
      WALLREFLECT = .FALSE.

      DES_PERIODIC_WALLS = .FALSE.
      DES_PERIODIC_WALLS_X = .FALSE.
      DES_PERIODIC_WALLS_Y = .FALSE.
      DES_PERIODIC_WALLS_Z = .FALSE.
      intx_per = .false.
      inty_per = .false.
      intz_per = .false.
      INLET_OUTLET = .FALSE.
      INLET_OUTLET_X = .FALSE.
      INLET_OUTLET_Y = .FALSE.
      INLET_OUTLET_Z = .FALSE.

      KN = UNDEFINED
      KT = UNDEFINED
      KT_FAC = UNDEFINED
      KN_W = UNDEFINED
      KT_W = UNDEFINED
      KT_W_FAC = UNDEFINED
      MEW = UNDEFINED
      MEW_W = UNDEFINED

      DES_EN_INPUT(:) = UNDEFINED
      DES_EN_WALL_INPUT(:) = UNDEFINED
      DES_ET_INPUT(:) = UNDEFINED
      DES_ET_WALL_INPUT(:) = UNDEFINED

      DES_ETAT_FAC = UNDEFINED
      DES_ETAT_W_FAC = UNDEFINED

      DES_GAMMA = ZERO
      DES_F = ZERO
      LID_VEL = 0.0d0

      S_TIME = UNDEFINED
      DES_SPX_TIME = UNDEFINED
      DTSOLID = UNDEFINED
      DTSOLID_FACTOR = 0.1D0
      P_TIME = UNDEFINED
      DEM_OUTPUT_DATA_TECPLOT = .FALSE.
      DEBUG_DES = .FALSE.
      IFI = 0

      DIMN = UNDEFINED_I  
      WX1 = UNDEFINED
      EX2 = UNDEFINED
      BY1 = UNDEFINED
      TY2 = UNDEFINED
      SZ1 = UNDEFINED
      NZ2 = UNDEFINED

      MQUAD_FACTOR = 1.1D0
      NQUAD = UNDEFINED_I
      INIT_QUAD_COUNT = UNDEFINED_I 
      INQC = UNDEFINED_I 
      QLM = 1       ! Number of levels to go up in the tree to move particle to a new quad 
      QLN = 1       ! Number of levels to go up in the tree to do neighbor search
      N2CT = UNDEFINED
      QUADCT = UNDEFINED
      OCTCT = UNDEFINED

      GENER_PART_CONFIG = .FALSE.
      VOL_FRAC(:) = UNDEFINED
      DES_EPS_XSTART = UNDEFINED
      DES_EPS_YSTART = UNDEFINED 
      DES_EPS_ZSTART = UNDEFINED 
      pvel_mean = zero 
      PVEL_StDev = zero 
      pgrad(:) = zero 

! J.Musser : mass inlet/outlet 
      MAX_PIS = UNDEFINED_I
      DES_BC_X_w(:) = UNDEFINED
      DES_BC_X_e(:) = UNDEFINED
      DES_BC_Y_s(:) = UNDEFINED
      DES_BC_Y_n(:) = UNDEFINED
      DES_BC_Z_b(:) = UNDEFINED
      DES_BC_Z_t(:) = UNDEFINED
      DES_BC_VOLFLOW_s(:) = UNDEFINED
      DES_BC_MASSFLOW_s(:) = UNDEFINED
      DES_BC_U_s(:) = UNDEFINED
      DES_BC_V_s(:) = UNDEFINED
      DES_BC_W_s(:) = UNDEFINED
      DES_BC_TYPE(:) = UNDEFINED_C

      COLL_MODEL = UNDEFINED_C
! T.Li : Hertzian collision model
      ew_young = undefined
      vw_poisson = undefined
      e_young(:) = undefined
      v_poisson(:) = undefined

      RETURN
      END SUBROUTINE DES_INIT_NAMELIST
