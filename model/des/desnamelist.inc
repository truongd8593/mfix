       NAMELIST / DES_INPUT_DATA /&
!
! USE discretelement
	DISCRETE_ELEMENT, DES_CONTINUUM_COUPLED, DES_INTERP_ON, &
        TSUJI_DRAG, DES_INTG_METHOD, &
	PARTICLES, DIMN, NFACTOR, &
	GENER_PART_CONFIG, VOL_FRAC, DES_EPS_XSTART,&
       	DES_EPS_YSTART, DES_EPS_ZSTART, &
        pgrad, pvel_mean, pvel_stdev, &
	MN, PARTICLES_FACTOR, &
	DES_NEIGHBOR_SEARCH, NEIGHBOR_SEARCH_N,&
        FACTOR_RLM, NEIGHBOR_SEARCH_RAD_RATIO,&
! 
! Used for quadtree or octree neighbor search methods
	MQUAD_FACTOR, QLM, QLN, INIT_QUAD_COUNT,&
! Particle properties	
        KN, KN_W, KT_FAC, KT_W_FAC, MEW, MEW_W, & 
        DES_EN_INPUT, DES_ET_INPUT, DES_ETAT_FAC, &
        DES_EN_WALL_INPUT, DES_ET_WALL_INPUT, DES_ETAT_W_FAC, &
!      
! Additional properties for Hertzian model: T.Li
        DES_COLL_MODEL, ew_young, vw_poisson, e_young, v_poisson,&
!	
! Not presently used
	DTSOLID_FACTOR,&
!
! Wall conditions
        DES_PERIODIC_WALLS, DES_PERIODIC_WALLS_X,&
        DES_PERIODIC_WALLS_Y, DES_PERIODIC_WALLS_Z,&
        WALLDTSPLIT,&
! Not presently used
        WALLREFLECT, NON_RECT_BC,& 
!
! Use either DES wall vibration parameters 
	DES_F, DES_GAMMA, LID_VEL,&
!
! or DES wall boundary condition information
        DES_BC_Uw_s, DES_BC_Vw_s, DES_BC_Ww_s,&
!	
! DES flow boundary condition information: J.Musser
        MAX_PIS, &
        DES_BC_X_w, DES_BC_X_e, &
        DES_BC_Y_s, DES_BC_Y_n, &
        DES_BC_Z_b, DES_BC_Z_t, &
        DES_BC_VOLFLOW_s, DES_BC_MASSFLOW_s, &
        DES_BC_U_s, DES_BC_V_s, DES_BC_W_s, DES_BC_TYPE, &
!
! Obsolete variables
        INLET_OUTLET, INLET_OUTLET_X, INLET_OUTLET_Y, INLET_OUTLET_Z,&
!
! DES output control information
        DEBUG_DES, FOCUS_PARTICLE,&
        PRINT_DES_DATA, P_TIME,&
        DEM_OUTPUT_DATA_TECPLOT, DES_SPX_TIME,&
!
! DES Cohesion Information: M.Weber
	MASTER_WELL_DEPTH, MASTER_WALL_WELL_DEPTH,&
        RADIUS_RATIO, WALL_RADIUS_RATIO, USE_COHESION, SQUARE_WELL,&
	MAX_LOG_STEP, VAN_DER_WAALS, HAMAKER_CONSTANT, VDW_INNER_CUTOFF,& 
        VDW_OUTER_CUTOFF, WALL_HAMAKER_CONSTANT,&
        WALL_VDW_INNER_CUTOFF, WALL_VDW_OUTER_CUTOFF,COHESION_DEBUG,&
        MAX_LOG_STEP, MIN_LOG_STEP, COHESION_DEBUG_START,&
!
! Obsolete Controls
        ANIMATION_WRITE_INTERVAL, ANIMATION_WRITE_COUNTER, TIME_FACTOR,&
        MAX_ANIMATOR_STEP, MIN_ANIMATOR_STEP
!


