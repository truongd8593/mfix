!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C 
!     Module name: DES_INIT_NAMELIST                                      C
!     Purpose: DES - initialize the des-namelist                          
!                                                                         C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!     Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!     Comments: Added some interpolation based inputs                     C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   
      SUBROUTINE DES_INIT_NAMELIST 

      USE param1
      USE discretelement
      USE mfix_pic
      USE des_bc
      USE des_ic
      USE des_thermo
      USE des_rxns
      
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------

      INCLUDE 'des/desnamelist.inc'

      PARTICLES = UNDEFINED_I

      MN = 10
      NFACTOR = 10
      PARTICLES_FACTOR = 1.2D0

      DISCRETE_ELEMENT = .FALSE.
      DES_CONTINUUM_COUPLED = .FALSE.
      DES_INTERP_ON = .FALSE.
      DES_INTERP_MEAN_FIELDS = .TRUE.
      DES_REPORT_MASS_INTERP = .false. 
      TSUJI_DRAG = .FALSE.
      DES_INTG_METHOD = 'EULER'
      USE_COHESION = .FALSE.
      SQUARE_WELL = .FALSE.
      VAN_DER_WAALS = .FALSE.
      DES_CONTINUUM_HYBRID = .FALSE.

      DES_NEIGHBOR_SEARCH = 1
      DO_NSEARCH = .FALSE.
      NEIGHBOR_SEARCH_N = 25
      NEIGHBOR_SEARCH_RAD_RATIO = 1.0D0
      FACTOR_RLM = 1.2 

      DESGRIDSEARCH_IMAX = UNDEFINED_I
      DESGRIDSEARCH_JMAX = UNDEFINED_I
      DESGRIDSEARCH_KMAX = UNDEFINED_I      
      MQUAD_FACTOR = 1.1D0
      NQUAD = UNDEFINED_I
      INIT_QUAD_COUNT = UNDEFINED_I 
      INQC = UNDEFINED_I 
      QLM = 1       ! Number of levels to go up in the tree to move particle to a new quad 
      QLN = 1       ! Number of levels to go up in the tree to do neighbor search
      N2CT = UNDEFINED
      QUADCT = UNDEFINED
      OCTCT = UNDEFINED

! particle properties      
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

      DES_COLL_MODEL = UNDEFINED_C
! T.Li : Hertzian collision model
      ew_young = undefined
      vw_poisson = undefined
      e_young(:) = undefined
      v_poisson(:) = undefined

      DES_D_P0(:) = UNDEFINED
      DES_RO_S(:) = UNDEFINED
      DES_MMAX = UNDEFINED_I

      DIMN = UNDEFINED_I  
      WX1 = UNDEFINED
      EX2 = UNDEFINED
      BY1 = UNDEFINED
      TY2 = UNDEFINED
      SZ1 = UNDEFINED
      NZ2 = UNDEFINED

      DES_PERIODIC_WALLS = .FALSE.
      DES_PERIODIC_WALLS_X = .FALSE.
      DES_PERIODIC_WALLS_Y = .FALSE.
      DES_PERIODIC_WALLS_Z = .FALSE.

! lees edwards bc keywords      
      DES_LE_BC = .FALSE.
      DES_LE_REL_VEL = UNDEFINED
      DES_LE_SHEAR_DIR = UNDEFINED_C

      WALLDTSPLIT = .FALSE.
      WALLREFLECT = .FALSE.

! des wall boundaries: vibrational parameters      
      DES_GAMMA = ZERO
      DES_F = ZERO
      LID_VEL = 0.0d0

! des wall boundaries: wall velocities
      DES_BC_Uw_s(:,:) = ZERO
      DES_BC_Vw_s(:,:) = ZERO
      DES_BC_Ww_s(:,:) = ZERO

! J.Musser : des mass inlet/outlet (flow) boundaries
      MAX_PIS = UNDEFINED_I
      DES_BC_X_w(:) = UNDEFINED
      DES_BC_X_e(:) = UNDEFINED
      DES_BC_Y_s(:) = UNDEFINED
      DES_BC_Y_n(:) = UNDEFINED
      DES_BC_Z_b(:) = UNDEFINED
      DES_BC_Z_t(:) = UNDEFINED
      DES_BC_VOLFLOW_s(:,:) = UNDEFINED
      DES_BC_MASSFLOW_s(:,:) = UNDEFINED
      DES_BC_TYPE(:) = UNDEFINED_C
      DES_BC_ROP_s(:,:) = UNDEFINED
      DES_BC_T_s(:,:) = UNDEFINED
      DES_BC_X_s(:,:,:) = UNDEFINED
      FORCE_ORD_BC = .FALSE.

! for cohesion:
      COHESION_DEBUG = UNDEFINED_I
      
! for cohesion: squarewell
      MASTER_WELL_DEPTH = UNDEFINED
      MASTER_WALL_WELL_DEPTH = UNDEFINED
      RADIUS_RATIO = UNDEFINED
      WALL_RADIUS_RATIO = UNDEFINED
! for cohesion: van der waals
      HAMAKER_CONSTANT = UNDEFINED
      WALL_HAMAKER_CONSTANT = UNDEFINED
      VDW_OUTER_CUTOFF = UNDEFINED
      VDW_INNER_CUTOFF = UNDEFINED
      WALL_VDW_OUTER_CUTOFF = ZERO ! important to set to zero in case of no cohesion.
      WALL_VDW_INNER_CUTOFF = UNDEFINED
      Asperities = ZERO

! J.Musser : des particle initial conditions
      DES_IC_X_w(:) = UNDEFINED
      DES_IC_X_e(:) = UNDEFINED
      DES_IC_Y_s(:) = UNDEFINED
      DES_IC_Y_n(:) = UNDEFINED
      DES_IC_Z_b(:) = UNDEFINED
      DES_IC_Z_t(:) = UNDEFINED
      DES_IC_T_s(:,:) = UNDEFINED
      DES_IC_X_s(:,:,:) = UNDEFINED

      DTSOLID = UNDEFINED
      DTSOLID_FACTOR = 0.1D0

      PRINT_DES_DATA = .FALSE.
      DES_SPX_DT = LARGE_NUMBER
      DES_RES_DT = LARGE_NUMBER
      DES_OUTPUT_TYPE = UNDEFINED_C
      DEBUG_DES = .FALSE.
      FOCUS_PARTICLE = 0
      VTP_FINDEX = 0
      TECPLOT_FINDEX = 0      

! flag whether calculate bed height of a phase      
      DES_CALC_BEDHEIGHT = .FALSE.

! flag whether to invoke cluster identification algorithm
      DES_CALC_CLUSTER = .FALSE.
      CLUSTER_LENGTH_CUTOFF = UNDEFINED      
 
! flag whether to use automatic particle configuration
      GENER_PART_CONFIG = .FALSE.
      VOL_FRAC(:) = UNDEFINED
      DES_EPS_XSTART = UNDEFINED
      DES_EPS_YSTART = UNDEFINED 
      DES_EPS_ZSTART = UNDEFINED 
      pvel_mean = zero 
      PVEL_StDev = zero 

      MAX_DES_BC_CELL = 4
      MPPIC_SOLID_STRESS_SNIDER = .false. 
      MPPIC_COEFF_EN = UNDEFINED
      MPPIC_COEFF_EN2 = UNDEFINED 
      MPPIC_COEFF_EN_WALL = UNDEFINED
      MPPIC_COEFF_ET_WALL = UNDEFINED

      MPPIC_PDRAG_IMPLICIT = .false. 
      MPPIC_GRAV_TREATMENT = .true. 
      
      PSFAC_FRIC_PIC = 100 
      FRIC_EXP_PIC = 2.5
      FRIC_NON_SING_FAC = 1E-07
      CFL_PIC = 0.1 

! J.Musser : des energy equations
      DES_ENERGY_EQ = .FALSE.
      DES_CONV_EQ = .TRUE.
      DES_COND_EQ = .TRUE.
      DES_RADI_EQ = .FALSE.
      DES_COND_EQ_PFP = .TRUE.
      DES_COND_EQ_PP  = .TRUE.

      DES_CONV_CORR = 'RANZ_1952'

      DES_MIN_COND_DIST = UNDEFINED
      FLPC = 1.0d0/5.0d0

      DES_K_s0(:) = UNDEFINED
      DES_C_ps0(:) = UNDEFINED
      DES_Em(:) = UNDEFINED

! J.Musser : species equations (reactive chemistry)
      DES_SPECIES_EQ(:) = .FALSE.
      DES_NMAX_s(:) = UNDEFINED_I 
      DES_MW_s(:,:) = UNDEFINED
      DES_SPECIES_s(:,:) = UNDEFINED_C
      DES_SPECIES_ALIAS_s(:,:) = UNDEFINED_C
      REACTION_MODEL = 'VARIABLE_DENSITY'

      RETURN
      END SUBROUTINE DES_INIT_NAMELIST
