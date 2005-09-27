      NAMELIST / INPUT_DATA /&
!
! 1.1 Run control section
!
             RUN_NAME, DESCRIPTION, UNITS, RUN_TYPE, TIME, TSTOP,&
             DT, ENERGY_EQ, DEF_COR, C_FAC, FPFOI, SPECIES_EQ, &
             CALL_USR, MODEL_B, DT_MAX, DT_MIN, DT_FAC, DETECT_STALL,&
             MOMENTUM_X_EQ, MOMENTUM_Y_EQ,MOMENTUM_Z_EQ, DISCRETIZE,&
             GRANULAR_ENERGY,& 
	     NScalar, Phase4Scalar, &
	     K_Epsilon, Call_DQMOM,&
! sof 02-01-05
	     SIMONIN, AHMADI, JENKINS, &
             nRR, chi_scheme, &
! start anuj 4/20
             FRICTION, SAVAGE,&
! end anuj 4/20
!
! sof (02/16/2005)
             SCHAEFFER,&
! end sof

! start peter 7/15
             SHEAR, V_sh,&
! end peter 7/15
             DRAG_TYPE, &

!AE TIME 043001  added new logical flag var for 2nd order time into NAMELIST
             CN_ON,&
!
! 1.2 Output control section
!
             RES_DT, SPX_DT , OUT_DT, &
             USR_DT, USR_TYPE, USR_VAR, USR_FORMAT, USR_EXT,&
             USR_X_w , USR_X_e , USR_Y_s , USR_Y_n ,&
             USR_Z_b , USR_Z_t, NLOG, FULL_LOG, RESID_STRING,&
	     report_mass_balance_dt, &
!
! 1.3 Physical and numerical parameters section
!            (e and C_e are same)
             e, C_e, C_f, C, C_NAME, Phi, Phi_w, L_scale0, MU_gmax,&
             TOL_RESID, TOL_DIVERGE, TOL_RESID_Scalar, MAX_NIT, LEQ_IT, LEQ_METHOD, LEQ_SWEEP,&
	     LEQ_TOL,  V_ex,&
             P_ref, P_scale, UR_FAC, UR_F_gs, GRAVITY, TOL_RESID_T, TOL_RESID_X,&
	     TOL_RESID_Th, MAX_INLET_VEL_FAC, &
             NORM_g, NORM_s, PHIP, e_w, EPS_f_min, EPS_max,&
	     drag_c1, drag_d1, &
!
! start GERA 8/15/03
             SEGREGATION_SLOPE_COEFFICIENT, ep_s_max,&
! modified by sof, removed MAX_SOLID_1_PACKING and MAX_SOLID_2_PACKING
! end GERA 8/15/03     
!
! start sof 5/04/05 (for rearanging indices from coarsest to finest particles)
             M_MAX,&
! end sof 5/04/05      
!
! start sof 5/04/05 Adding automatic_restart
             AUTO_RESTART,&
! end sof 5/04/05    
!	     	     
! start rong
             AGGREGATION_EFF,BREAKAGE_EFF,&
! end rong 
!
! 1.4 Geometry and discretization section
!
             COORDINATES, IMAX, DX, XMIN, XLENGTH, JMAX, DY, YLENGTH,&
             KMAX, DZ, ZLENGTH, NO_I, NO_J, NO_K,&
             CYCLIC_X, CYCLIC_Y, CYCLIC_Z, CYCLIC_X_PD, CYCLIC_Y_PD,&
             CYCLIC_Z_PD, &
!
! 1.5 Gas phase section
!
             RO_g0, MU_g0, MW_AVG, MW_g, K_g0, C_pg0, DIF_g0,&
!
! 1.6 Solids phase section
!
             MU_s0, D_p0, RO_s, EP_star, MMAX, NMAX, MW_s, CLOSE_PACKED,&
             K_s0, C_ps0, DIF_s0,&

!
! 1.7 Initial conditions section
!
             IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t,&
             IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t,&
             IC_EP_g, IC_P_g, IC_ROP_s, IC_T_g, IC_T_s,&
             IC_U_g, IC_U_s, IC_V_g, IC_V_s, IC_W_g, IC_W_s,&
             IC_X_g, IC_X_s, IC_TYPE, IC_P_star, IC_L_scale,&
             IC_GAMA_Rg, IC_T_Rg, IC_GAMA_Rs, IC_T_Rs,&
             IC_Theta_m, IC_Scalar, IC_K_Turb_G, IC_E_Turb_G, &
!
! 1.8 Boundary conditions section
!
             BC_X_w, BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b, BC_Z_t,&
             BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t,&
             BC_EP_g, BC_P_g, BC_RO_g, BC_ROP_g, BC_ROP_s,&
             BC_T_g, BC_T_s,&
             BC_U_g, BC_U_s, BC_V_g, BC_V_s, BC_W_g, BC_W_s,&
             BC_TYPE, BC_VOLFLOW_g, BC_VOLFLOW_s, BC_MASSFLOW_g, &
             BC_MASSFLOW_s, BC_X_g, BC_X_s,&
             BC_DT_0, BC_Jet_g0, BC_DT_h, BC_Jet_gh, BC_DT_l, &
             BC_Jet_gl,&
             DELP_X, DELP_Y, DELP_Z, U_g0, U_s0, V_g0, V_s0, W_g0, W_s0,&
             BC_hw_g, BC_hw_s,&
             BC_Uw_g, BC_Uw_s, BC_Vw_g, BC_Vw_s, BC_Ww_g, BC_Ww_s,&
             BC_hw_T_g, BC_hw_T_s, BC_Tw_g, BC_Tw_s, BC_C_T_g,&
             BC_C_T_s,&
             BC_hw_X_g, BC_hw_X_s, BC_Xw_g, BC_Xw_s, BC_C_X_g,&
             BC_C_X_s,&
!start kapil anuj 01/19/98
             BC_JJ_PS,&
             BC_hw_Theta_m, BC_Thetaw_m, BC_C_Theta_m,&
             BC_Theta_m,&
!end   kapil anuj 01/19/98
             BC_hw_Scalar, BC_Scalarw, BC_C_Scalar,&
             BC_Scalar,&
	     BC_K_Turb_G, BC_E_Turb_G, &
	     Flux_g, &
!
! 1.9 Internal surfaces section
!
             IS_X_w, IS_X_e, IS_Y_s, IS_Y_n, IS_Z_b, IS_Z_t,&
             IS_I_w, IS_I_e, IS_J_s, IS_J_n, IS_K_b, IS_K_t,&
             IS_TYPE, IS_PC, IS_VEL_s,&
!
! 1.10 Chemical reactions
!
             SPECIES_NAME, &
!
!     CHEM & ISAT begin (nan xie)
!
	     CALL_DI, CALL_GROW, CALL_ISAT, ISATdt,&
!     CHEM & ISAT end (nan xie)
! 1.11 distributed parallel
!
             NODESI , NODESJ , NODESK, &
!
! 1.12 Preconditioner
!
	     LEQ_PC, &

! 1.13 Parallel related

	     IS_SERIAL, USE_DOLOOP
