       NAMELIST / DES_INPUT_DATA /&

!	    USE discretelement
! start Jay 
        MAXNEIGHBORS, MAXQUADS, NMQD, NWALLS, PBP, DES_SPX_TIME,&	
        DIMN, S_TIME, DTSOLID, DTSOLID_FACTOR, P_TIME, PTC,&
        PARTICLES, PARTICLES_FACTOR, MQUAD_FACTOR,&
        KN, KN_W, KT, KT_W, NEIGHBOR_SEARCH_N,&
        NEIGHBOR_SEARCH_RAD_RATIO, NEIGHBOR_SEARCH_DIST,&
        ETA_DES_N, ETA_N_W, ETA_DES_T, ETA_T_W,& 
        MEW, MEW_W, E_RESTITUTION, WALLCONTACT,&
        WX1, EX2, BY1, TY2, SZ1, NZ2, DES_GAMMA, DES_F,&
        DES_NEIGHBOR_SEARCH, MN, NQUAD,&
        QLM, QLN, INIT_QUAD_COUNT, INQC, RADIUS_EQ,&
        N2CT, NBSCT, QUADCT, OCTCT, DES_KE, DES_PE,&
        DISCRETE_ELEMENT, DES_CONTINUUM_COUPLED, PARTICLE_SLIDE,&
        DO_QUADTREE, DO_OCTREE, DO_NSQUARE, DO_NSEARCH,&
        WALLFIXEDOVERLAP, WALLDTSPLIT, WALLREFLECT, DES_PERIODIC_WALLS,& 
	DES_PERIODIC_WALLS_X, DES_PERIODIC_WALLS_Y, DES_PERIODIC_WALLS_Z,&
        INLET_OUTLET, INLET_OUTLET_X, INLET_OUTLET_Y, INLET_OUTLET_Z,&
        NON_RECT_BC,&
        TSUJI_DRAG, IFI, PRINT_DES_DATA, &
	MASTER_WELL_DEPTH, MASTER_WALL_WELL_DEPTH,&
        RADIUS_RATIO,WALL_RADIUS_RATIO, USE_COHESION, SQUARE_WELL,&
	MAX_LOG_STEP, VAN_DER_WAALS, HAMAKER_CONSTANT, VDW_INNER_CUTOFF,& 
        VDW_OUTER_CUTOFF, WALL_HAMAKER_CONSTANT,&
        WALL_VDW_INNER_CUTOFF, WALL_VDW_OUTER_CUTOFF,COHESION_DEBUG,FOCUS_PARTICLE,&
        MAX_LOG_STEP, MIN_LOG_STEP, COHESION_DEBUG_START,&
        ANIMATION_WRITE_INTERVAL, ANIMATION_WRITE_COUNTER, TIME_FACTOR,&
	 MAX_ANIMATOR_STEP, MIN_ANIMATOR_STEP, NFACTOR, DES_INTERP_ON 
! end Jay 

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DESNAMELIST.INC                                        C
!  Purpose: Variables needed for DES input                             C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      
