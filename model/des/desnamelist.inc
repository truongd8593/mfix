       NAMELIST / DES_INPUT_DATA /&

! USE discretelement
	DISCRETE_ELEMENT, DES_CONTINUUM_COUPLED, DES_INTERP_ON, &
	DES_ONEWAY_COUPLED, DES_INTERP_MEAN_FIELDS, & 
	 DES_REPORT_MASS_INTERP, USE_STL_DES, &
      DES_INTG_METHOD, DES_CONTINUUM_HYBRID,&
	 PARTICLES, DIMN, NFACTOR, USE_VDH_DEM_MODEL, &
	GENER_PART_CONFIG, MN, PARTICLES_FACTOR, &
	DES_NEIGHBOR_SEARCH, NEIGHBOR_SEARCH_N,&
        FACTOR_RLM, NEIGHBOR_SEARCH_RAD_RATIO,&
	
! Used for grid/cell based search methods
        DESGRIDSEARCH_IMAX, DESGRIDSEARCH_JMAX, DESGRIDSEARCH_KMAX, &

! Particle properties	
        KN, KN_W, KT_FAC, KT_W_FAC, MEW, MEW_W, & 
        DES_EN_INPUT, DES_ET_INPUT, DES_ETAT_FAC, &
        DES_EN_WALL_INPUT, DES_ET_WALL_INPUT, DES_ETAT_W_FAC, &
        DES_PSI_s, &
!      
! Additional properties for Hertzian model: 
        DES_COLL_MODEL, EW_YOUNG, VW_POISSON, E_YOUNG, V_POISSON,&
	
! Lees Edwards BC (for DEM only)
        DES_LE_BC, DES_LE_REL_VEL, DES_LE_SHEAR_DIR,&
	

! DES wall BC: wall translational velocity
! These need to be improved to work for arbitrary geoms
        DES_BC_Uw_s, DES_BC_Vw_s, DES_BC_Ww_s,&
	
! DES flow boundary condition information: J.Musser
        MAX_PIS, &
        DES_BC_X_w, DES_BC_X_e, &
        DES_BC_Y_s, DES_BC_Y_n, &
        DES_BC_Z_b, DES_BC_Z_t, &
        DES_BC_VOLFLOW_s, DES_BC_MASSFLOW_s, &
        DES_BC_ROP_s, DES_BC_TYPE, &
        DES_BC_T_s, DES_BC_X_s, &
       	FORCE_ORD_BC, &

! DES initial condition information: J.Musser
        DES_IC_X_w, DES_IC_X_e, &
        DES_IC_Y_s, DES_IC_Y_n, &
        DES_IC_Z_b, DES_IC_Z_t, &
        DES_IC_T_s, DES_IC_X_s, &

! DES output control information
        DEBUG_DES, FOCUS_PARTICLE,&
        PRINT_DES_DATA, DES_SPX_DT, DES_RES_DT, &
        DES_OUTPUT_TYPE,  &

! DES Cohesion Information: M.Weber
        USE_COHESION,VAN_DER_WAALS, HAMAKER_CONSTANT, & 
	WALL_HAMAKER_CONSTANT, &
        VDW_INNER_CUTOFF, VDW_OUTER_CUTOFF, ASPERITIES, &
        WALL_VDW_INNER_CUTOFF, WALL_VDW_OUTER_CUTOFF,&

! Des cluster identification
        DES_CALC_CLUSTER, CLUSTER_LENGTH_CUTOFF, &


! DES Energy Equations: J.Musser
        DES_ENERGY_EQ, DES_CONV_EQ, DES_COND_EQ, DES_RADI_EQ, &
        DES_COND_EQ_PFP, DES_COND_EQ_PP, DES_CONV_CORR, FLPC, & 
        DES_Em, DES_MIN_COND_DIST, &

! DES Species Equations (reactive chemistry): J. Musser
       DES_NMAX_s,  DES_MW_s, DES_SPECIES_EQ, &
       DES_SPECIES_s, DES_SPECIES_ALIAS_s, REACTION_MODEL, &

! MP-PIC related flags 
	 MPPIC, MPPIC_SOLID_STRESS_SNIDER, PSFAC_FRIC_PIC,& 
	 FRIC_EXP_PIC, FRIC_NON_SING_FAC,&
	 CFL_PIC, MPPIC_COEFF_EN1, MPPIC_COEFF_EN2, MPPIC_PDRAG_IMPLICIT, &
	 MPPIC_GRAV_TREATMENT, MPPIC_COEFF_EN_WALL, MPPIC_COEFF_ET_WALL, &

! des_cellwise_bc_data size 
	 MAX_DES_BC_CELL, & 

! flag to convert des outside box to facets needed if outside box not part of stl geom

	 des_convert_box_to_facets

