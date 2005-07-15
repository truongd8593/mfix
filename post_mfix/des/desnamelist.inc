       NAMELIST / DES_INPUT_DATA /&

!	    USE discretelement
! start Dhanunjay 

             DO_NSQUARE, OUTOFBOX, DIMN, WALLDTSPLIT, COUPLED_FLOW,&
             DTSOLID, PARTICLES, N2CT, NBSCT, QUADCT, OCTCT, DES_GAMMA,&
             FACTOR, NEIGHBOURS, FC, CALLED, WALLFIXEDOVERLAP,& 
             FN, FT, TOW, DISCRETE_ELEMENT, DES_POS_NEW, DES_VEL_NEW,&
             KN, KT, ETA_DES_N, MEW, DES_POS_OLD, DES_VEL_OLD,&
             OMEGA_OLD, OMEGA_NEW, ETA_DES_T, RO_Sol, DES_NEIGHBOR_SEARCH,&
             DO_NBS, DO_OCTREE, PFN, PFT, PN, PV, TIME_ADJUST,&
             DES_WALL_POS, DES_WALL_VEL, DES_PERIODIC_WALLS_X, WX1, EX2,&
             BY1, TY2, SZ1, NZ2, MN, PMASS, MOI, GRAV, ROs,&
             DES_RADIUS, PR, RADIUS_EQ, DO_QUADTREE, INLET_OUTLET,&
	     LQUAD, CQUAD, NQUAD, WALLREFLECT, E_RESTITUTION, FTS1,&
	     WALL_NORMAL,KN_W, KT_W, ETA_N_W, ETA_T_W, MEW_W, FNS1,&
	     DES_KE, DES_PE, DES_F, ZONES, ZN1, ZN2, ZN3, ZN4,&
	     DES_CONTINUUM_COUPLED, TSUJI_DRAG, EQUIVALENT_RADIUS,&
	     KUIPERS_DRAG, DES_PERIODIC_WALLS_Y, DES_PERIODIC_WALLS_Z,&
	     DES_PERIODIC_WALLS, MASTER_WELL_DEPTH, MASTER_WALL_WELL_DEPTH,&
             RADIUS_RATIO,WALL_RADIUS_RATIO, USE_COHESION, SQUARE_WELL,&
             MAX_LOG_STEP, VAN_DER_WAALS, HAMAKER_CONSTANT, VDW_INNER_CUTOFF,& 
             VDW_OUTER_CUTOFF, WALL_HAMAKER_CONSTANT,&
             WALL_VDW_INNER_CUTOFF, WALL_VDW_OUTER_CUTOFF,COHESION_DEBUG,FOCUS_PARTICLE,&
             MAX_LOG_STEP, MIN_LOG_STEP, COHESION_DEBUG_START,&
             ANIMATION_WRITE_INTERVAL, ANIMATION_WRITE_COUNTER, TIME_FACTOR,&
             MAX_ANIMATOR_STEP, MIN_ANIMATOR_STEP
! end Dhanunjay 

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
