" Vim syntax file
" Language: mfix input file: mfix.dat
" Version:  1.0

if exists("b:current_syntax")
  finish
endif

syn case ignore



" list of MFIX keywords
" This list was generated from all *namelist.inc files
" and sorted alphabetically.

syn keyword mfixKeyword	ALPHA_MAX
syn keyword mfixKeyword	ASPERITIES
syn keyword mfixKeyword	AUTO_RESTART
syn keyword mfixKeyword	BAR_CHAR
syn keyword mfixKeyword	BAR_RESOLUTION
syn keyword mfixKeyword	BAR_WIDTH
syn keyword mfixKeyword	BATCH_WALLCLOCK
syn keyword mfixKeyword	BC_APPLY_TO_MPPIC
syn keyword mfixKeyword	BC_C_SCALAR
syn keyword mfixKeyword	BC_C_T_G
syn keyword mfixKeyword	BC_C_THETA_M
syn keyword mfixKeyword	BC_C_T_S
syn keyword mfixKeyword	BC_C_X_G
syn keyword mfixKeyword	BC_C_X_S
syn keyword mfixKeyword	BC_DT_0
syn keyword mfixKeyword	BC_DT_H
syn keyword mfixKeyword	BC_DT_L
syn keyword mfixKeyword	BC_EP_G
syn keyword mfixKeyword	BC_E_TURB_G
syn keyword mfixKeyword	BC_HW_G
syn keyword mfixKeyword	BC_HW_S
syn keyword mfixKeyword	BC_HW_SCALAR
syn keyword mfixKeyword	BC_HW_T_G
syn keyword mfixKeyword	BC_HW_THETA_M
syn keyword mfixKeyword	BC_HW_T_S
syn keyword mfixKeyword	BC_HW_X_G
syn keyword mfixKeyword	BC_HW_X_S
syn keyword mfixKeyword	BC_ID_Q
syn keyword mfixKeyword	BC_I_E
syn keyword mfixKeyword	BC_I_W
syn keyword mfixKeyword	BC_JET_G0
syn keyword mfixKeyword	BC_JET_GH
syn keyword mfixKeyword	BC_JET_GL
syn keyword mfixKeyword	BC_JJ_PS
syn keyword mfixKeyword	BC_J_N
syn keyword mfixKeyword	BC_J_S
syn keyword mfixKeyword	BC_K_B
syn keyword mfixKeyword	BC_K_T
syn keyword mfixKeyword	BC_K_TURB_G
syn keyword mfixKeyword	BC_MASSFLOW_G
syn keyword mfixKeyword	BC_MASSFLOW_S
syn keyword mfixKeyword	BC_P_G
syn keyword mfixKeyword	BC_RO_G
syn keyword mfixKeyword	BC_ROP_G
syn keyword mfixKeyword	BC_ROP_S
syn keyword mfixKeyword	BC_SCALAR
syn keyword mfixKeyword	BC_SCALARW
syn keyword mfixKeyword	BC_T_G
syn keyword mfixKeyword	BC_THETA_M
syn keyword mfixKeyword	BC_THETAW_M
syn keyword mfixKeyword	BC_T_S
syn keyword mfixKeyword	BC_TW_G
syn keyword mfixKeyword	BC_TW_S
syn keyword mfixKeyword	BC_TYPE
syn keyword mfixKeyword	BC_U_G
syn keyword mfixKeyword	BC_U_S
syn keyword mfixKeyword	BC_UW_G
syn keyword mfixKeyword	BC_UW_S
syn keyword mfixKeyword	BC_VELMAG_G
syn keyword mfixKeyword	BC_VELMAG_S
syn keyword mfixKeyword	BC_V_G
syn keyword mfixKeyword	BC_VOLFLOW_G
syn keyword mfixKeyword	BC_VOLFLOW_S
syn keyword mfixKeyword	BC_V_S
syn keyword mfixKeyword	BC_VW_G
syn keyword mfixKeyword	BC_VW_S
syn keyword mfixKeyword	BC_W_G
syn keyword mfixKeyword	BC_W_S
syn keyword mfixKeyword	BC_WW_G
syn keyword mfixKeyword	BC_WW_S
syn keyword mfixKeyword	BC_X_E
syn keyword mfixKeyword	BC_X_G
syn keyword mfixKeyword	BC_X_S
syn keyword mfixKeyword	BC_X_W
syn keyword mfixKeyword	BC_XW_G
syn keyword mfixKeyword	BC_XW_S
syn keyword mfixKeyword	BC_Y_N
syn keyword mfixKeyword	BC_Y_S
syn keyword mfixKeyword	BC_Z_B
syn keyword mfixKeyword	BC_Z_T
syn keyword mfixKeyword	BLENDING_STRESS
syn keyword mfixKeyword	BREAKAGE_EFF
syn keyword mfixKeyword	BWRITE_NETCDF
syn keyword mfixKeyword	C
syn keyword mfixKeyword	CAD_PROPAGATE_ORDER
syn keyword mfixKeyword	CALL_DI
syn keyword mfixKeyword	CALL_DQMOM
syn keyword mfixKeyword	CALL_GROW
syn keyword mfixKeyword	CALL_ISAT
syn keyword mfixKeyword	CALL_USR
syn keyword mfixKeyword	CARTESIAN_GRID
syn keyword mfixKeyword	C_E
syn keyword mfixKeyword	C_F
syn keyword mfixKeyword	C_FAC
syn keyword mfixKeyword	CFL_PIC
syn keyword mfixKeyword		CFL_PIC
syn keyword mfixKeyword	CG_SAFE_MODE
syn keyword mfixKeyword	CG_UR_FAC
syn keyword mfixKeyword	CHI_SCHEME
syn keyword mfixKeyword	CHK_BATCHQ_END
syn keyword mfixKeyword	CLIP_XMAX
syn keyword mfixKeyword	CLIP_XMIN
syn keyword mfixKeyword	CLIP_YMAX
syn keyword mfixKeyword	CLIP_YMIN
syn keyword mfixKeyword	CLIP_ZMAX
syn keyword mfixKeyword	CLIP_ZMIN
syn keyword mfixKeyword	CLOSE_PACKED
syn keyword mfixKeyword	C_NAME
syn keyword mfixKeyword	CN_ON
syn keyword mfixKeyword		COHESION_DEBUG
syn keyword mfixKeyword	COHESION_DEBUG_START
syn keyword mfixKeyword	CONSTANTNPC
syn keyword mfixKeyword	CONSTANTWT
syn keyword mfixKeyword	COORDINATES
syn keyword mfixKeyword	C_PG0
syn keyword mfixKeyword	C_PS0
syn keyword mfixKeyword	CPX
syn keyword mfixKeyword	CPY
syn keyword mfixKeyword	CPZ
syn keyword mfixKeyword	CYCLIC_X
syn keyword mfixKeyword	CYCLIC_X_PD
syn keyword mfixKeyword	CYCLIC_Y
syn keyword mfixKeyword	CYCLIC_Y_PD
syn keyword mfixKeyword	CYCLIC_Z
syn keyword mfixKeyword	CYCLIC_Z_PD
syn keyword mfixKeyword	DBGPRN_LAYOUT
syn keyword mfixKeyword	DEBUG_DES
syn keyword mfixKeyword	DEBUG_RESID
syn keyword mfixKeyword	DEF_COR
syn keyword mfixKeyword	DELP_X
syn keyword mfixKeyword	DELP_Y
syn keyword mfixKeyword	DELP_Z
syn keyword mfixKeyword	DES_BC_CORE_RHO
syn keyword mfixKeyword	DES_BC_MASSFLOW_S
syn keyword mfixKeyword	DES_BC_ROP_S
syn keyword mfixKeyword	DES_BC_T_S
syn keyword mfixKeyword	DES_BC_TYPE
syn keyword mfixKeyword	DES_BC_UW_S
syn keyword mfixKeyword	DES_BC_VOLFLOW_S
syn keyword mfixKeyword	DES_BC_VW_S
syn keyword mfixKeyword	DES_BC_WW_S
syn keyword mfixKeyword	DES_BC_X_E
syn keyword mfixKeyword	DES_BC_X_S
syn keyword mfixKeyword	DES_BC_X_W
syn keyword mfixKeyword	DES_BC_Y_N
syn keyword mfixKeyword	DES_BC_Y_S
syn keyword mfixKeyword	DES_BC_Z_B
syn keyword mfixKeyword	DES_BC_Z_T
syn keyword mfixKeyword	DES_CALC_BEDHEIGHT
syn keyword mfixKeyword	DES_COLL_MODEL
syn keyword mfixKeyword	DES_COND_EQ
syn keyword mfixKeyword	DES_COND_EQ_PFP
syn keyword mfixKeyword	DES_COND_EQ_PP
syn keyword mfixKeyword	DES_CONTINUUM_COUPLED
syn keyword mfixKeyword	DES_CONTINUUM_HYBRID
syn keyword mfixKeyword	DES_CONV_CORR
syn keyword mfixKeyword	DES_CONV_EQ
syn keyword mfixKeyword	DES_C_PS0
syn keyword mfixKeyword	DESCRIPTION
syn keyword mfixKeyword	DES_D_P0
syn keyword mfixKeyword	DES_EM
syn keyword mfixKeyword	DES_ENERGY_EQ
syn keyword mfixKeyword	DES_EN_INPUT
syn keyword mfixKeyword	DES_EN_WALL_INPUT
syn keyword mfixKeyword	DES_EPS_XSTART
syn keyword mfixKeyword		DES_EPS_YSTART
syn keyword mfixKeyword	DES_EPS_ZSTART
syn keyword mfixKeyword	DES_ETAT_FAC
syn keyword mfixKeyword	DES_ETAT_W_FAC
syn keyword mfixKeyword	DES_ET_INPUT
syn keyword mfixKeyword	DES_ET_WALL_INPUT
syn keyword mfixKeyword		DES_F
syn keyword mfixKeyword	DES_GAMMA
syn keyword mfixKeyword	DESGRIDSEARCH_IMAX
syn keyword mfixKeyword	DESGRIDSEARCH_JMAX
syn keyword mfixKeyword		DESGRIDSEARCH_KMAX
syn keyword mfixKeyword	DES_IC_CORE_RHO
syn keyword mfixKeyword	DES_IC_T_S
syn keyword mfixKeyword	DES_IC_X_E
syn keyword mfixKeyword	DES_IC_X_S
syn keyword mfixKeyword	DES_IC_X_W
syn keyword mfixKeyword	DES_IC_Y_N
syn keyword mfixKeyword	DES_IC_Y_S
syn keyword mfixKeyword	DES_IC_Z_B
syn keyword mfixKeyword	DES_IC_Z_T
syn keyword mfixKeyword	DES_INTERP_ON
syn keyword mfixKeyword	DES_INTG_METHOD
syn keyword mfixKeyword	DES_K_S0
syn keyword mfixKeyword	DES_LE_BC
syn keyword mfixKeyword	DES_LE_REL_VEL
syn keyword mfixKeyword	DES_LE_SHEAR_DIR
syn keyword mfixKeyword	DES_MIN_COND_DIST
syn keyword mfixKeyword	DES_MMAX
syn keyword mfixKeyword	DES_MW_S
syn keyword mfixKeyword		DES_NEIGHBOR_SEARCH
syn keyword mfixKeyword	DES_NMAX
syn keyword mfixKeyword		DES_ONEWAY_COUPLED
syn keyword mfixKeyword	DES_OUTPUT_TYPE
syn keyword mfixKeyword	DES_PERIODIC_WALLS
syn keyword mfixKeyword	DES_PERIODIC_WALLS_X
syn keyword mfixKeyword	DES_PERIODIC_WALLS_Y
syn keyword mfixKeyword	DES_PERIODIC_WALLS_Z
syn keyword mfixKeyword	DES_RADI_EQ
syn keyword mfixKeyword	DES_RES_DT
syn keyword mfixKeyword	DES_RO_S
syn keyword mfixKeyword	DES_SPECIES_EQ
syn keyword mfixKeyword	DES_SPECIES_NAME
syn keyword mfixKeyword	DES_SPX_DT
syn keyword mfixKeyword	DETECT_STALL
syn keyword mfixKeyword	DIF_G0
syn keyword mfixKeyword	DIF_S0
syn keyword mfixKeyword	DIMN
syn keyword mfixKeyword		DISCRETE_ELEMENT
syn keyword mfixKeyword	DISCRETIZE
syn keyword mfixKeyword	DO_TRANSPOSE
syn keyword mfixKeyword	D_P0
syn keyword mfixKeyword	DQUADRIC
syn keyword mfixKeyword	DRAG_C1
syn keyword mfixKeyword	DRAG_D1
syn keyword mfixKeyword	DRAG_TYPE
syn keyword mfixKeyword	DT
syn keyword mfixKeyword	DT_FAC
syn keyword mfixKeyword	DT_MAX
syn keyword mfixKeyword	DT_MIN
syn keyword mfixKeyword		DTSOLID_FACTOR
syn keyword mfixKeyword	DX
syn keyword mfixKeyword	DY
syn keyword mfixKeyword	DZ
syn keyword mfixKeyword	E
syn keyword mfixKeyword	ENABLE_DMP_LOG
syn keyword mfixKeyword	ENERGY_EQ
syn keyword mfixKeyword	EPS_F_MIN
syn keyword mfixKeyword	EP_S_MAX
syn keyword mfixKeyword	EP_STAR
syn keyword mfixKeyword	ERX
syn keyword mfixKeyword	ERY
syn keyword mfixKeyword	ERZ
syn keyword mfixKeyword	E_W
syn keyword mfixKeyword	EW_YOUNG
syn keyword mfixKeyword	E_YOUNG
syn keyword mfixKeyword	FAC_DIM_MAX_CUT_CELL
syn keyword mfixKeyword	FACTOR_RLM
syn keyword mfixKeyword	F_DASHBOARD
syn keyword mfixKeyword	FEDORS_LANDEL
syn keyword mfixKeyword	FIRST_DX
syn keyword mfixKeyword	FIRST_DY
syn keyword mfixKeyword	FIRST_DZ
syn keyword mfixKeyword	FLPC
syn keyword mfixKeyword	FLUID_IN_CLIPPED_REGION
syn keyword mfixKeyword	FLUX_G
syn keyword mfixKeyword	FOCUS_PARTICLE
syn keyword mfixKeyword		FORCE_ORD_BC
syn keyword mfixKeyword	FPFOI
syn keyword mfixKeyword	FRAME
syn keyword mfixKeyword		FRIC_EXP_PIC
syn keyword mfixKeyword	FRIC_NON_SING_FAC
syn keyword mfixKeyword	FRICTION
syn keyword mfixKeyword	FULL_LOG
syn keyword mfixKeyword		GENER_PART_CONFIG
syn keyword mfixKeyword	GRANULAR_ENERGY
syn keyword mfixKeyword	GRAVITY
syn keyword mfixKeyword	GROUP_Q
syn keyword mfixKeyword	GROUP_RELATION
syn keyword mfixKeyword	GROUP_SIZE
syn keyword mfixKeyword	HALF_ANGLE
syn keyword mfixKeyword	HAMAKER_CONSTANT
syn keyword mfixKeyword	IC_EP_G
syn keyword mfixKeyword	IC_E_TURB_G
syn keyword mfixKeyword	IC_GAMA_RG
syn keyword mfixKeyword	IC_GAMA_RS
syn keyword mfixKeyword	ICHECK_BICGS
syn keyword mfixKeyword	IC_I_E
syn keyword mfixKeyword	IC_I_W
syn keyword mfixKeyword	IC_J_N
syn keyword mfixKeyword	IC_J_S
syn keyword mfixKeyword	IC_K_B
syn keyword mfixKeyword	IC_K_T
syn keyword mfixKeyword	IC_K_TURB_G
syn keyword mfixKeyword	IC_L_SCALE
syn keyword mfixKeyword	IC_P_G
syn keyword mfixKeyword	IC_P_STAR
syn keyword mfixKeyword	IC_ROP_S
syn keyword mfixKeyword	IC_SCALAR
syn keyword mfixKeyword	IC_T_G
syn keyword mfixKeyword	IC_THETA_M
syn keyword mfixKeyword	IC_T_RG
syn keyword mfixKeyword	IC_T_RS
syn keyword mfixKeyword	IC_T_S
syn keyword mfixKeyword	IC_TYPE
syn keyword mfixKeyword	IC_U_G
syn keyword mfixKeyword	IC_U_S
syn keyword mfixKeyword	IC_V_G
syn keyword mfixKeyword	IC_V_S
syn keyword mfixKeyword	IC_W_G
syn keyword mfixKeyword	IC_W_S
syn keyword mfixKeyword	IC_X_E
syn keyword mfixKeyword	IC_X_G
syn keyword mfixKeyword	IC_X_S
syn keyword mfixKeyword	IC_X_W
syn keyword mfixKeyword	IC_Y_N
syn keyword mfixKeyword	IC_Y_S
syn keyword mfixKeyword	IC_Z_B
syn keyword mfixKeyword	IC_Z_T
syn keyword mfixKeyword	IMAX
syn keyword mfixKeyword	INIT_QUAD_COUNT
syn keyword mfixKeyword	ISATDT
syn keyword mfixKeyword	IS_I_E
syn keyword mfixKeyword	IS_I_W
syn keyword mfixKeyword	IS_J_N
syn keyword mfixKeyword	IS_J_S
syn keyword mfixKeyword	IS_K_B
syn keyword mfixKeyword	IS_K_T
syn keyword mfixKeyword	IS_PC
syn keyword mfixKeyword	IS_SERIAL
syn keyword mfixKeyword	IS_TYPE
syn keyword mfixKeyword	IS_VEL_S
syn keyword mfixKeyword	IS_X_E
syn keyword mfixKeyword	IS_X_W
syn keyword mfixKeyword	IS_Y_N
syn keyword mfixKeyword	IS_Y_S
syn keyword mfixKeyword	IS_Z_B
syn keyword mfixKeyword	IS_Z_T
syn keyword mfixKeyword	ITERMAX_INT
syn keyword mfixKeyword	JENKINS
syn keyword mfixKeyword	JMAX
syn keyword mfixKeyword	K_EPSILON
syn keyword mfixKeyword	K_G0
syn keyword mfixKeyword	KMAX
syn keyword mfixKeyword	KN
syn keyword mfixKeyword	KN_W
syn keyword mfixKeyword	K_S0
syn keyword mfixKeyword	KT_FAC
syn keyword mfixKeyword	KT_TYPE
syn keyword mfixKeyword	KT_W_FAC
syn keyword mfixKeyword	LAMBDA_X
syn keyword mfixKeyword	LAMBDA_Y
syn keyword mfixKeyword	LAMBDA_Z
syn keyword mfixKeyword	LAM_HYS
syn keyword mfixKeyword	LAST_DX
syn keyword mfixKeyword	LAST_DY
syn keyword mfixKeyword	LAST_DZ
syn keyword mfixKeyword	LEQ_IT
syn keyword mfixKeyword	LEQ_METHOD
syn keyword mfixKeyword	LEQ_PC
syn keyword mfixKeyword	LEQ_SWEEP
syn keyword mfixKeyword	LEQ_TOL
syn keyword mfixKeyword	LID_VEL
syn keyword mfixKeyword	L_SCALE0
syn keyword mfixKeyword	M_AM
syn keyword mfixKeyword		MASTER_WALL_WELL_DEPTH
syn keyword mfixKeyword	MASTER_WELL_DEPTH
syn keyword mfixKeyword		MAX_DES_BC_CELL
syn keyword mfixKeyword	MAX_INLET_VEL_FAC
syn keyword mfixKeyword	MAX_NIT
syn keyword mfixKeyword	MAX_PIS
syn keyword mfixKeyword	MEW
syn keyword mfixKeyword	MEW_W
syn keyword mfixKeyword	MINIMIZE_DOTPRODUCTS
syn keyword mfixKeyword	MMAX
syn keyword mfixKeyword	M_MAX
syn keyword mfixKeyword		MN
syn keyword mfixKeyword	MODEL_B
syn keyword mfixKeyword	MOMENTUM_X_EQ
syn keyword mfixKeyword	MOMENTUM_Y_EQ
syn keyword mfixKeyword	MOMENTUM_Z_EQ
syn keyword mfixKeyword		MPPIC
syn keyword mfixKeyword	MPPIC_COEFF_EN
syn keyword mfixKeyword	MPPIC_PDRAG_IMPLICIT
syn keyword mfixKeyword	MPPIC_SOLID_STRESS_SNIDER
syn keyword mfixKeyword		MQUAD_FACTOR
syn keyword mfixKeyword	MSH_SMALL_ANGLE
syn keyword mfixKeyword	MU_G0
syn keyword mfixKeyword	MU_GMAX
syn keyword mfixKeyword	MU_S0
syn keyword mfixKeyword	MW_AVG
syn keyword mfixKeyword	MW_G
syn keyword mfixKeyword	MW_S
syn keyword mfixKeyword	NAMELIST/QMOMK_INPUT_DATA/
syn keyword mfixKeyword	NCX
syn keyword mfixKeyword	NCY
syn keyword mfixKeyword	NCZ
syn keyword mfixKeyword	NEIGHBOR_SEARCH_N
syn keyword mfixKeyword	NEIGHBOR_SEARCH_RAD_RATIO
syn keyword mfixKeyword	NFACTOR
syn keyword mfixKeyword	N_GROUP
syn keyword mfixKeyword	NLOG
syn keyword mfixKeyword	NMAX
syn keyword mfixKeyword	NODESI
syn keyword mfixKeyword	NODESJ
syn keyword mfixKeyword	NODESK
syn keyword mfixKeyword	NO_I
syn keyword mfixKeyword	NO_J
syn keyword mfixKeyword	NO_K
syn keyword mfixKeyword	NORM_G
syn keyword mfixKeyword	NORM_S
syn keyword mfixKeyword		NPC_PIC
syn keyword mfixKeyword	N_QUADRIC
syn keyword mfixKeyword	NRR
syn keyword mfixKeyword	NSCALAR
syn keyword mfixKeyword	N_USR_DEF
syn keyword mfixKeyword	N_X
syn keyword mfixKeyword	N_Y
syn keyword mfixKeyword	N_Z
syn keyword mfixKeyword	OPT_PARALLEL
syn keyword mfixKeyword	OUT_DT
syn keyword mfixKeyword	OUT_MSH_VALUE
syn keyword mfixKeyword	OUT_STL_VALUE
syn keyword mfixKeyword		PARTICLES
syn keyword mfixKeyword	PARTICLES_FACTOR
syn keyword mfixKeyword	PG_OPTION
syn keyword mfixKeyword	PHASE4SCALAR
syn keyword mfixKeyword	PHI
syn keyword mfixKeyword	PHIP
syn keyword mfixKeyword	PHI_W
syn keyword mfixKeyword	PIECE_XMAX
syn keyword mfixKeyword	PIECE_XMIN
syn keyword mfixKeyword	PIECE_YMAX
syn keyword mfixKeyword	PIECE_YMIN
syn keyword mfixKeyword	PIECE_ZMAX
syn keyword mfixKeyword	PIECE_ZMIN
syn keyword mfixKeyword	P_REF
syn keyword mfixKeyword	PRINT_DES_DATA
syn keyword mfixKeyword	PRINT_PROGRESS_BAR
syn keyword mfixKeyword	PRINT_QMOMK_DATA
syn keyword mfixKeyword	PRINT_WARNINGS
syn keyword mfixKeyword	P_SCALE
syn keyword mfixKeyword	PSFAC_FRIC_PIC
syn keyword mfixKeyword	PVEL_MEAN
syn keyword mfixKeyword	PVEL_STDEV
syn keyword mfixKeyword	QLM
syn keyword mfixKeyword	QLN
syn keyword mfixKeyword	QMOMK
syn keyword mfixKeyword	QMOMK_CFL
syn keyword mfixKeyword	QMOMK_COLLISIONS
syn keyword mfixKeyword	QMOMK_COLLISIONS_ORDER
syn keyword mfixKeyword		QMOMK_COUPLED
syn keyword mfixKeyword	QMOMK_TYPE
syn keyword mfixKeyword	QMOMK_WALL_BC_TYPE
syn keyword mfixKeyword	QUADRIC_FORM
syn keyword mfixKeyword	QUADRIC_SCALE
syn keyword mfixKeyword	RADIUS
syn keyword mfixKeyword	RADIUS_RATIO
syn keyword mfixKeyword	RAY_DIR
syn keyword mfixKeyword	RDF_TYPE
syn keyword mfixKeyword	RDPC
syn keyword mfixKeyword	REACTION_MODEL
syn keyword mfixKeyword	RELATION_WITH_PREVIOUS
syn keyword mfixKeyword	REPORT_MASS_BALANCE_DT
syn keyword mfixKeyword	RES_DT
syn keyword mfixKeyword	RESID_STRING
syn keyword mfixKeyword	RO_G0
syn keyword mfixKeyword	RO_S
syn keyword mfixKeyword	R_P
syn keyword mfixKeyword	RUN_NAME
syn keyword mfixKeyword	RUN_TYPE
syn keyword mfixKeyword	SAVAGE
syn keyword mfixKeyword	SCALE_MSH
syn keyword mfixKeyword	SCALE_STL
syn keyword mfixKeyword	SCHAEFFER
syn keyword mfixKeyword	SEGREGATION_SLOPE_COEFFICIENT
syn keyword mfixKeyword	SET_CORNER_CELLS
syn keyword mfixKeyword	SHEAR
syn keyword mfixKeyword	SIGM_BLEND
syn keyword mfixKeyword	SIMONIN
syn keyword mfixKeyword	SOLVER_STATISTICS
syn keyword mfixKeyword	SPECIES_EQ
syn keyword mfixKeyword	SPECIES_NAME
syn keyword mfixKeyword	SPX_DT
syn keyword mfixKeyword	SQUARE_WELL
syn keyword mfixKeyword	STATWT_PIC
syn keyword mfixKeyword	STATWT_PIC
syn keyword mfixKeyword	STL_BC_ID
syn keyword mfixKeyword	STL_SMALL_ANGLE
syn keyword mfixKeyword	TANH_BLEND
syn keyword mfixKeyword	TERM_BUFFER
syn keyword mfixKeyword	THETA_X
syn keyword mfixKeyword	THETA_Y
syn keyword mfixKeyword	THETA_Z
syn keyword mfixKeyword	TIME
syn keyword mfixKeyword	TIME_DEPENDENT_FILENAME
syn keyword mfixKeyword	TOL_DELH
syn keyword mfixKeyword	TOL_DIVERGE
syn keyword mfixKeyword	TOL_F
syn keyword mfixKeyword	TOL_MERGE
syn keyword mfixKeyword	TOL_MSH
syn keyword mfixKeyword	TOL_POLY
syn keyword mfixKeyword	TOL_RESID
syn keyword mfixKeyword	TOL_RESID_SCALAR
syn keyword mfixKeyword	TOL_RESID_T
syn keyword mfixKeyword	TOL_RESID_TH
syn keyword mfixKeyword	TOL_RESID_X
syn keyword mfixKeyword	TOL_SMALL_AREA
syn keyword mfixKeyword	TOL_SMALL_CELL
syn keyword mfixKeyword	TOL_SNAP
syn keyword mfixKeyword	TOL_STL
syn keyword mfixKeyword	TSTOP
syn keyword mfixKeyword	TSUJI_DRAG
syn keyword mfixKeyword	T_X
syn keyword mfixKeyword	TX_MSH
syn keyword mfixKeyword	TX_STL
syn keyword mfixKeyword	T_Y
syn keyword mfixKeyword	TY_MSH
syn keyword mfixKeyword	TY_STL
syn keyword mfixKeyword	T_Z
syn keyword mfixKeyword	TZ_MSH
syn keyword mfixKeyword	TZ_STL
syn keyword mfixKeyword	U_G0
syn keyword mfixKeyword	UNITS
syn keyword mfixKeyword	UR_FAC
syn keyword mfixKeyword	UR_F_GS
syn keyword mfixKeyword	UR_KTH_SML
syn keyword mfixKeyword	U_S0
syn keyword mfixKeyword		USE_COHESION
syn keyword mfixKeyword	USE_DOLOOP
syn keyword mfixKeyword	USE_MSH
syn keyword mfixKeyword	USE_POLYGON
syn keyword mfixKeyword	USE_STL
syn keyword mfixKeyword	USR_DT
syn keyword mfixKeyword	USR_EXT
syn keyword mfixKeyword	USR_FORMAT
syn keyword mfixKeyword	USR_TYPE
syn keyword mfixKeyword	USR_VAR
syn keyword mfixKeyword	USR_X_E
syn keyword mfixKeyword	USR_X_W
syn keyword mfixKeyword	USR_Y_N
syn keyword mfixKeyword	USR_Y_S
syn keyword mfixKeyword	USR_Z_B
syn keyword mfixKeyword	USR_Z_T
syn keyword mfixKeyword		VAN_DER_WAALS
syn keyword mfixKeyword	VDW_INNER_CUTOFF
syn keyword mfixKeyword	VDW_OUTER_CUTOFF
syn keyword mfixKeyword	V_EX
syn keyword mfixKeyword	V_G0
syn keyword mfixKeyword	VOL_FRAC
syn keyword mfixKeyword	V_POISSON
syn keyword mfixKeyword	V_S0
syn keyword mfixKeyword	V_SH
syn keyword mfixKeyword	VTK_DT
syn keyword mfixKeyword	VTK_VAR
syn keyword mfixKeyword	VTU_DIR
syn keyword mfixKeyword	VW_POISSON
syn keyword mfixKeyword	WALLDTSPLIT
syn keyword mfixKeyword	WALL_HAMAKER_CONSTANT
syn keyword mfixKeyword	WALL_RADIUS_RATIO
syn keyword mfixKeyword	WALLREFLECT
syn keyword mfixKeyword	WALL_VDW_INNER_CUTOFF
syn keyword mfixKeyword	WALL_VDW_OUTER_CUTOFF
syn keyword mfixKeyword	W_G0
syn keyword mfixKeyword	WRITE_DASHBOARD
syn keyword mfixKeyword	WRITE_VTK_FILES
syn keyword mfixKeyword	W_S0
syn keyword mfixKeyword	XLENGTH
syn keyword mfixKeyword	XMIN
syn keyword mfixKeyword	YLENGTH
syn keyword mfixKeyword	YU_STANDISH
syn keyword mfixKeyword	ZLENGTH


" Detect type of text (comments, string of character, etc.)

syn match mfixBoolean "\.\s*\(true\|false\|t\|f\)\s*\."
syn match mfixComment      excludenl "!.*$" contains=@spell
syn region mfixString      start=+'+ end=+'+ contains=@spell
syn region mfixString      start=+"+ end=+"+ contains=@spell
syn region mfixBracket      start=+(+ end=+)+ contains=@spell

" Assign arbitrary color from predefined types

hi def link mfixKeyword    Keyword
hi def link mfixBoolean    Boolean
hi def link mfixComment    Comment
hi def link mfixString     Preproc
hi def link mfixBracket    Type 



let b:current_syntax = "mfix"

