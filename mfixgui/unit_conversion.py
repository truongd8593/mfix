#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division

""" Factors for converting CGS to SI models, etc
    converting CGS to SI models"""

cgs_to_SI = {
#<start generated code>
  'added_mass' : 1,
  'adjust_proc_domain_size' : 1,
# 'aggregation_eff': UNDEFINED,
# 'alpha_max': UNDEFINED,
# 'asperities': UNDEFINED,
  'auto_restart' : 1,
  'automatic_restart' : 1,
  'bar_char' : 1,
# 'bar_resolution': UNDEFINED,
  'bar_width' : 1,
# 'batch_wallclock': UNDEFINED,
# 'bc_c_scalar': UNDEFINED,
# 'bc_c_t_g': UNDEFINED,
# 'bc_c_t_s': UNDEFINED,
# 'bc_c_theta_m': UNDEFINED,
# 'bc_c_x_g': UNDEFINED,
# 'bc_c_x_s': UNDEFINED,
# 'bc_dt_0': UNDEFINED,
# 'bc_dt_h': UNDEFINED,
# 'bc_dt_l': UNDEFINED,
# 'bc_e_turb_g': UNDEFINED,
  'bc_ep_g' : 1,
  'bc_ep_s' : 1,
# 'bc_hw_g': UNDEFINED,
# 'bc_hw_s': UNDEFINED,
# 'bc_hw_scalar': UNDEFINED,
# 'bc_hw_t_g': UNDEFINED,
# 'bc_hw_t_s': UNDEFINED,
# 'bc_hw_theta_m': UNDEFINED,
# 'bc_hw_x_g': UNDEFINED,
# 'bc_hw_x_s': UNDEFINED,
  'bc_i_e' : 1,
  'bc_i_w' : 1,
  'bc_id_q' : 1,
  'bc_j_n' : 1,
  'bc_j_s' : 1,
# 'bc_jet_g0': UNDEFINED,
# 'bc_jet_gh': UNDEFINED,
# 'bc_jet_gl': UNDEFINED,
  'bc_jj_m' : 1,
  'bc_jj_ps' : 1,
  'bc_k_b' : 1,
  'bc_k_t' : 1,
# 'bc_k_turb_g': UNDEFINED,
# 'bc_massflow_g': UNDEFINED,
# 'bc_massflow_s': UNDEFINED,
# 'bc_p_g': UNDEFINED,
  'bc_pic_mi_const_npc' : 1,
# 'bc_pic_mi_const_statwt': UNDEFINED,
  'bc_po_apply_to_des' : 1,
# 'bc_rop_s': UNDEFINED,
# 'bc_scalar': UNDEFINED,
# 'bc_scalarw': UNDEFINED,
# 'bc_t_g': UNDEFINED,
# 'bc_t_s': UNDEFINED,
# 'bc_theta_m': UNDEFINED,
# 'bc_thetaw_m': UNDEFINED,
# 'bc_tw_g': UNDEFINED,
# 'bc_tw_s': UNDEFINED,
  'bc_type' : 1,
# 'bc_u_g': UNDEFINED,
# 'bc_u_s': UNDEFINED,
# 'bc_uw_g': UNDEFINED,
# 'bc_uw_s': UNDEFINED,
# 'bc_v_g': UNDEFINED,
# 'bc_v_s': UNDEFINED,
# 'bc_velmag_g': UNDEFINED,
# 'bc_velmag_s': UNDEFINED,
# 'bc_volflow_g': UNDEFINED,
# 'bc_volflow_s': UNDEFINED,
# 'bc_vw_g': UNDEFINED,
# 'bc_vw_s': UNDEFINED,
# 'bc_w_g': UNDEFINED,
# 'bc_w_s': UNDEFINED,
# 'bc_ww_g': UNDEFINED,
# 'bc_ww_s': UNDEFINED,
# 'bc_x_e': UNDEFINED,
  'bc_x_g' : 1,
  'bc_x_s' : 1,
# 'bc_x_w': UNDEFINED,
  'bc_xw_g' : 1,
  'bc_xw_s' : 1,
# 'bc_y_n': UNDEFINED,
# 'bc_y_s': UNDEFINED,
# 'bc_z_b': UNDEFINED,
# 'bc_z_t': UNDEFINED,
  'bdist_io' : 1,
# 'bend_r1': UNDEFINED,
# 'bend_r2': UNDEFINED,
# 'bend_theta1': UNDEFINED,
# 'bend_theta2': UNDEFINED,
  'blending_function' : 1,
# 'breakage_eff': UNDEFINED,
  'bstart_with_one_res' : 1,
  'bwrite_netcdf' : 1,
# 'c': UNDEFINED,
# 'c2c_r1': UNDEFINED,
# 'c2c_r2': UNDEFINED,
  'c2c_y1' : 1,
# 'c2c_y2': UNDEFINED,
# 'c_e': UNDEFINED,
# 'c_f': UNDEFINED,
# 'c_fac': UNDEFINED,
  'c_name' : 1,
  'c_pg0' : 4184.0,
  'c_ps0' : 4184.0,
  'cad_propagate_order' : 1,
  'call_dqmom' : 1,
  'call_usr' : 1,
  'call_usr_source' : 1,
  'cartesian_grid' : 1,
# 'cfl_pic': UNDEFINED,
  'cg_safe_mode' : 1,
# 'cg_ur_fac': UNDEFINED,
  'chi_scheme' : 1,
  'chk_batchq_end' : 1,
# 'clip_xmax': UNDEFINED,
# 'clip_xmin': UNDEFINED,
# 'clip_ymax': UNDEFINED,
# 'clip_ymin': UNDEFINED,
# 'clip_zmax': UNDEFINED,
# 'clip_zmin': UNDEFINED,
  'close_packed' : 1,
  'cn_on' : 1,
  'coordinates' : 1,
# 'cpx': UNDEFINED,
# 'cpy': UNDEFINED,
# 'cpz': UNDEFINED,
  'cyclic_x' : 1,
  'cyclic_x_pd' : 1,
  'cyclic_y' : 1,
  'cyclic_y_pd' : 1,
  'cyclic_z' : 1,
  'cyclic_z_pd' : 1,
  'cylindrical_2d' : 1,
  'd_p0' : 0.01,
  'dbgprn_layout' : 1,
  'debug_des' : 1,
  'debug_resid' : 1,
  'def_cor' : 1,
# 'delp_x': UNDEFINED,
# 'delp_y': UNDEFINED,
# 'delp_z': UNDEFINED,
  'des_coll_model' : 1,
  'des_conv_corr' : 1,
# 'des_diffuse_width': UNDEFINED,
# 'des_em': UNDEFINED,
# 'des_en_input': UNDEFINED,
# 'des_en_wall_input': UNDEFINED,
# 'des_et_input': UNDEFINED,
# 'des_et_wall_input': UNDEFINED,
# 'des_etat_fac': UNDEFINED,
# 'des_etat_w_fac': UNDEFINED,
  'des_explicitly_coupled' : 1,
  'des_interp_mean_fields' : 1,
  'des_interp_on' : 1,
  'des_interp_scheme' : 1,
# 'des_interp_width': UNDEFINED,
  'des_intg_method' : 1,
# 'des_min_cond_dist': UNDEFINED,
  'des_neighbor_search' : 1,
  'des_oneway_coupled' : 1,
  'des_output_type' : 1,
  'des_report_mass_interp' : 1,
  'des_usr_var_size' : 1,
  'description' : 1,
  'desgridsearch_imax' : 1,
  'desgridsearch_jmax' : 1,
  'desgridsearch_kmax' : 1,
  'detect_stall' : 1,
  'dif_g0' : 0.0001,
  'dif_s0' : 0.0001,
  'dil_factor_vsd' : 1,
  'dil_inert_x_vsd' : 1,
  'dim_facets_per_cell' : 1,
  'discretize' : 1,
  'do_transpose' : 1,
# 'dquadric': UNDEFINED,
# 'drag_c1': UNDEFINED,
# 'drag_d1': UNDEFINED,
  'drag_type' : 1,
# 'dt': UNDEFINED,
# 'dt_fac': UNDEFINED,
# 'dt_max': UNDEFINED,
# 'dt_min': UNDEFINED,
  'dwall_brute_force' : 1,
# 'dx': UNDEFINED,
# 'dy': UNDEFINED,
# 'dz': UNDEFINED,
# 'e_w': UNDEFINED,
  'e_young' : 0.1,
  'e_young_actual' : 0.1,
  'enable_dmp_log' : 1,
  'energy_eq' : 1,
  'ep_s_max' : 1,
  'ep_star' : 1,
  'eps_f_min' : 1,
# 'erx': UNDEFINED,
# 'ery': UNDEFINED,
# 'erz': UNDEFINED,
  'ew_young' : 0.1,
  'ew_young_actual' : 0.1,
  'f_dashboard' : 1,
  'fac_dim_max_cut_cell' : 1,
# 'factor_rlm': UNDEFINED,
  'fedors_landel' : 1,
# 'filter_size_ratio': UNDEFINED,
# 'first_dx': UNDEFINED,
# 'first_dy': UNDEFINED,
# 'first_dz': UNDEFINED,
# 'flpc': UNDEFINED,
  'fluid_in_clipped_region' : 1,
# 'flux_g': UNDEFINED,
  'focus_particle' : 1,
  'fpfoi' : 1,
  'frame' : 1,
# 'fric_exp_pic': UNDEFINED,
  'fric_non_sing_fac' : 1,
  'friction_model' : 1,
  'full_log' : 1,
  'gener_part_config' : 1,
  'gravity' : 0.01,
  'gravity_x' : 0.01,
  'gravity_y' : 0.01,
  'gravity_z' : 0.01,
  'group_q' : 1,
  'group_relation' : 1,
  'group_resid' : 1,
  'group_size' : 1,
# 'half_angle': UNDEFINED,
# 'hamaker_constant': UNDEFINED,
  'i_cyl_num' : 1,
  'i_cyl_transition' : 1,
  'ic_des_fit_to_region' : 1,
# 'ic_e_turb_g': UNDEFINED,
  'ic_ep_g' : 1,
  'ic_ep_s' : 1,
# 'ic_gama_rg': UNDEFINED,
# 'ic_gama_rs': UNDEFINED,
  'ic_i_e' : 1,
  'ic_i_w' : 1,
  'ic_j_n' : 1,
  'ic_j_s' : 1,
  'ic_k_b' : 1,
  'ic_k_t' : 1,
# 'ic_k_turb_g': UNDEFINED,
# 'ic_l_scale': UNDEFINED,
# 'ic_p_g': UNDEFINED,
# 'ic_p_star': UNDEFINED,
  'ic_pic_const_npc' : 1,
# 'ic_pic_const_statwt': UNDEFINED,
# 'ic_rop_s': UNDEFINED,
# 'ic_scalar': UNDEFINED,
# 'ic_t_g': UNDEFINED,
# 'ic_t_rg': UNDEFINED,
# 'ic_t_rs': UNDEFINED,
# 'ic_t_s': UNDEFINED,
# 'ic_theta_m': UNDEFINED,
  'ic_type' : 1,
# 'ic_u_g': UNDEFINED,
# 'ic_u_s': UNDEFINED,
# 'ic_v_g': UNDEFINED,
# 'ic_v_s': UNDEFINED,
# 'ic_w_g': UNDEFINED,
# 'ic_w_s': UNDEFINED,
# 'ic_x_e': UNDEFINED,
  'ic_x_g' : 1,
  'ic_x_s' : 1,
# 'ic_x_w': UNDEFINED,
# 'ic_y_n': UNDEFINED,
# 'ic_y_s': UNDEFINED,
# 'ic_z_b': UNDEFINED,
# 'ic_z_t': UNDEFINED,
  'icheck_bicgs' : 1,
  'imax' : 1,
  'inert_species' : 1,
  'is_i_e' : 1,
  'is_i_w' : 1,
  'is_j_n' : 1,
  'is_j_s' : 1,
  'is_k_b' : 1,
  'is_k_t' : 1,
# 'is_pc': UNDEFINED,
  'is_serial' : 1,
  'is_type' : 1,
# 'is_vel_s': UNDEFINED,
# 'is_x_e': UNDEFINED,
# 'is_x_w': UNDEFINED,
# 'is_y_n': UNDEFINED,
# 'is_y_s': UNDEFINED,
# 'is_z_b': UNDEFINED,
# 'is_z_t': UNDEFINED,
  'ishii' : 1,
  'iter_restart' : 1,
  'itermax_int' : 1,
  'jackson' : 1,
  'jenkins' : 1,
  'jmax' : 1,
  'k_g0' : 418.4,
  'k_s0' : 418.4,
  'kmax' : 1,
  'kn' : 0.001,
  'kn_w' : 0.001,
# 'kt_fac': UNDEFINED,
  'kt_type' : 1,
# 'kt_w_fac': UNDEFINED,
# 'lam_hys': UNDEFINED,
# 'lambda_x': UNDEFINED,
# 'lambda_y': UNDEFINED,
# 'lambda_z': UNDEFINED,
# 'last_dx': UNDEFINED,
# 'last_dy': UNDEFINED,
# 'last_dz': UNDEFINED,
  'leq_it' : 1,
  'leq_method' : 1,
  'leq_pc' : 1,
  'leq_sweep' : 1,
# 'leq_tol': UNDEFINED,
  'm_am' : 1,
# 'max_inlet_vel_fac': UNDEFINED,
  'max_nit' : 1,
# 'mew': UNDEFINED,
# 'mew_w': UNDEFINED,
  'minimize_des_facet_list' : 1,
  'minimize_send_recv' : 1,
  'mmax' : 1,
  'model_b' : 1,
  'momentum_x_eq' : 1,
  'momentum_y_eq' : 1,
  'momentum_z_eq' : 1,
# 'mppic_coeff_en1': UNDEFINED,
# 'mppic_coeff_en2': UNDEFINED,
# 'mppic_coeff_en_wall': UNDEFINED,
# 'mppic_coeff_et_wall': UNDEFINED,
  'mppic_grav_treatment' : 1,
  'mppic_pdrag_implicit' : 1,
  'mppic_solid_stress_snider' : 1,
  'mu_g0' : 0.1,
# 'mu_gmax': UNDEFINED,
# 'mu_s0': UNDEFINED,
  'mw_avg' : 1.0,
  'mw_g' : 1.0,
  'mw_s' : 1.0,
  'n_group' : 1,
  'n_quadric' : 1,
  'n_usr_def' : 1,
# 'n_x': UNDEFINED,
# 'n_y': UNDEFINED,
# 'n_z': UNDEFINED,
  'ncx' : 1,
  'ncy' : 1,
  'ncz' : 1,
  'neighbor_search_n' : 1,
# 'neighbor_search_rad_ratio': UNDEFINED,
  'nfactor' : 1,
  'nlog' : 1,
  'nmax' : 1,
  'nmax_g' : 1,
  'nmax_s' : 1,
  'no_k' : 1,
  'nodesi' : 1,
  'nodesi_report' : 1,
  'nodesj' : 1,
  'nodesj_report' : 1,
  'nodesk' : 1,
  'nodesk_report' : 1,
# 'norm_g': UNDEFINED,
  'norm_s' : 1,
  'nrr' : 1,
  'nscalar' : 1,
  'opt_parallel' : 1,
# 'out_dt': UNDEFINED,
# 'out_msh_value': UNDEFINED,
# 'out_stl_value': UNDEFINED,
# 'p_ref': UNDEFINED,
# 'p_scale': UNDEFINED,
  'particles' : 1,
  'persistent_mode' : 1,
  'pg_option' : 1,
  'phase4scalar' : 1,
# 'phi': UNDEFINED,
# 'phi_w': UNDEFINED,
# 'phip': UNDEFINED,
# 'phip0': UNDEFINED,
  'phip_out_jj' : 1,
  'pic_report_deletion_stats' : 1,
  'pic_report_min_epg' : 1,
  'pic_report_seeding_stats' : 1,
# 'piece_xmax': UNDEFINED,
# 'piece_xmin': UNDEFINED,
# 'piece_ymax': UNDEFINED,
# 'piece_ymin': UNDEFINED,
# 'piece_zmax': UNDEFINED,
# 'piece_zmin': UNDEFINED,
  'print_des_data' : 1,
  'print_progress_bar' : 1,
  'print_warnings' : 1,
  'ps_i_e' : 1,
  'ps_i_w' : 1,
  'ps_j_n' : 1,
  'ps_j_s' : 1,
  'ps_k_b' : 1,
  'ps_k_t' : 1,
# 'ps_massflow_g': UNDEFINED,
# 'ps_massflow_s': UNDEFINED,
# 'ps_t_g': UNDEFINED,
# 'ps_t_s': UNDEFINED,
# 'ps_u_g': UNDEFINED,
# 'ps_u_s': UNDEFINED,
# 'ps_v_g': UNDEFINED,
# 'ps_v_s': UNDEFINED,
# 'ps_w_g': UNDEFINED,
# 'ps_w_s': UNDEFINED,
# 'ps_x_e': UNDEFINED,
  'ps_x_g' : 1,
  'ps_x_s' : 1,
# 'ps_x_w': UNDEFINED,
# 'ps_y_n': UNDEFINED,
# 'ps_y_s': UNDEFINED,
# 'ps_z_b': UNDEFINED,
# 'ps_z_t': UNDEFINED,
  'psfac_fric_pic' : 1,
  'quadric_form' : 1,
# 'quadric_scale': UNDEFINED,
# 'r_p': UNDEFINED,
# 'radius': UNDEFINED,
  'ray_dir' : 1,
  'rdf_type' : 1,
  're_indexing' : 1,
# 'reactor1_r1': UNDEFINED,
# 'reactor1_r2': UNDEFINED,
# 'reactor1_rr1': UNDEFINED,
# 'reactor1_rr2': UNDEFINED,
# 'reactor1_theta1': UNDEFINED,
# 'reactor1_theta2': UNDEFINED,
# 'reactor1_y1': UNDEFINED,
# 'reactor1_y2': UNDEFINED,
# 'reactor1_yr1': UNDEFINED,
# 'reactor1_yr2': UNDEFINED,
  'relation_with_previous' : 1,
  'report_best_domain_size' : 1,
# 'report_mass_balance_dt': UNDEFINED,
  'report_neg_density' : 1,
  'report_neg_specificheat' : 1,
# 'res_backup_dt': UNDEFINED,
  'res_backups' : 1,
# 'res_dt': UNDEFINED,
  'resid_string' : 1,
  'ro_g0' : 1000.0,
  'ro_s0' : 1000.0,
  'ro_xs0' : 1000.0,
  'run_name' : 1,
  'run_type' : 1,
# 'scale_msh': UNDEFINED,
# 'scale_stl': UNDEFINED,
# 'segregation_slope_coefficient': UNDEFINED,
  'set_corner_cells' : 1,
  'shear' : 1,
  'solids_model' : 1,
  'solver_statistics' : 1,
  'species_alias_g' : 1,
  'species_alias_s' : 1,
  'species_eq' : 1,
  'species_g' : 1,
  'species_name' : 1,
  'species_s' : 1,
  'spx_dt' : 1,
  'stiff_chem_max_steps' : 1,
  'stiff_chemistry' : 1,
  'stl_bc_id' : 1,
# 'stl_small_angle': UNDEFINED,
  'subgrid_type' : 1,
  'subgrid_wall' : 1,
# 't_x': UNDEFINED,
# 't_y': UNDEFINED,
# 't_z': UNDEFINED,
# 'term_buffer': UNDEFINED,
# 'theta_x': UNDEFINED,
# 'theta_y': UNDEFINED,
# 'theta_z': UNDEFINED,
# 'time': UNDEFINED,
  'time_dependent_filename' : 1,
  'tol_delh' : 1,
# 'tol_diverge': UNDEFINED,
# 'tol_f': UNDEFINED,
  'tol_merge' : 1,
# 'tol_msh': UNDEFINED,
# 'tol_poly': UNDEFINED,
# 'tol_resid': UNDEFINED,
# 'tol_resid_k_epsilon': UNDEFINED,
# 'tol_resid_scalar': UNDEFINED,
# 'tol_resid_t': UNDEFINED,
# 'tol_resid_th': UNDEFINED,
# 'tol_resid_x': UNDEFINED,
  'tol_small_area' : 1,
  'tol_small_cell' : 1,
  'tol_snap' : 1,
# 'tol_stl': UNDEFINED,
# 'tol_stl_dp': UNDEFINED,
# 'torus_r1': UNDEFINED,
# 'torus_r2': UNDEFINED,
# 'tstop': UNDEFINED,
  'turbulence_model' : 1,
# 'tx_msh': UNDEFINED,
# 'tx_stl': UNDEFINED,
# 'ty_msh': UNDEFINED,
# 'ty_stl': UNDEFINED,
# 'tz_msh': UNDEFINED,
# 'tz_stl': UNDEFINED,
# 'ucoil_r1': UNDEFINED,
# 'ucoil_r2': UNDEFINED,
  'ucoil_y1' : 1,
# 'ucoil_y2': UNDEFINED,
  'units' : 1,
# 'ur_f_gs': UNDEFINED,
# 'ur_fac': UNDEFINED,
# 'ur_kth_sml': UNDEFINED,
  'use_cohesion' : 1,
  'use_doloop' : 1,
  'use_msh' : 1,
  'use_polygon' : 1,
  'use_rrates' : 1,
  'use_stl' : 1,
  'use_vdh_dem_model' : 1,
  'usr_cpg' : 1,
  'usr_cps' : 1,
  'usr_difg' : 1,
  'usr_difs' : 1,
# 'usr_dt': UNDEFINED,
  'usr_ext' : 1,
  'usr_fgs' : 1,
  'usr_format' : 1,
  'usr_fss' : 1,
  'usr_gama' : 1,
  'usr_i_e' : 1,
  'usr_i_w' : 1,
  'usr_j_n' : 1,
  'usr_j_s' : 1,
  'usr_k_b' : 1,
  'usr_k_t' : 1,
  'usr_kg' : 1,
  'usr_ks' : 1,
  'usr_mug' : 1,
  'usr_mus' : 1,
  'usr_rog' : 1,
  'usr_ros' : 1,
  'usr_type' : 1,
  'usr_var' : 1,
# 'usr_x_e': UNDEFINED,
# 'usr_x_w': UNDEFINED,
# 'usr_y_n': UNDEFINED,
# 'usr_y_s': UNDEFINED,
# 'usr_z_b': UNDEFINED,
# 'usr_z_t': UNDEFINED,
# 'v_ex': UNDEFINED,
# 'v_poisson': UNDEFINED,
# 'v_poisson_actual': UNDEFINED,
# 'v_sh': UNDEFINED,
  'van_der_waals' : 1,
# 'vdw_inner_cutoff': UNDEFINED,
# 'vdw_outer_cutoff': UNDEFINED,
  'vtk_bc_id' : 1,
  'vtk_cutcell_only' : 1,
  'vtk_data' : 1,
  'vtk_dbg_file' : 1,
  'vtk_debug' : 1,
# 'vtk_dt': UNDEFINED,
  'vtk_dwall' : 1,
  'vtk_e_turb_g' : 1,
  'vtk_ep_g' : 1,
  'vtk_facet_count_des' : 1,
  'vtk_ijk' : 1,
  'vtk_k_turb_g' : 1,
  'vtk_lambda_2' : 1,
  'vtk_nb_facet_des' : 1,
  'vtk_normal' : 1,
  'vtk_p_g' : 1,
  'vtk_p_star' : 1,
  'vtk_part_angular_vel' : 1,
  'vtk_part_cohesion' : 1,
  'vtk_part_diameter' : 1,
  'vtk_part_orientation' : 1,
  'vtk_part_rank' : 1,
  'vtk_part_temp' : 1,
  'vtk_part_usr_var' : 1,
  'vtk_part_vel' : 1,
  'vtk_part_x_s' : 1,
  'vtk_partition' : 1,
  'vtk_rop_s' : 1,
  'vtk_rrate' : 1,
  'vtk_scalar' : 1,
  'vtk_select_mode' : 1,
# 'vtk_slice_tol': UNDEFINED,
  'vtk_t_g' : 1,
  'vtk_t_s' : 1,
  'vtk_theta_m' : 1,
  'vtk_u_g' : 1,
  'vtk_u_s' : 1,
  'vtk_v_g' : 1,
  'vtk_v_s' : 1,
  'vtk_var' : 1,
  'vtk_vel_g' : 1,
  'vtk_vel_s' : 1,
  'vtk_vorticity' : 1,
  'vtk_w_g' : 1,
  'vtk_w_s' : 1,
# 'vtk_x_e': UNDEFINED,
  'vtk_x_g' : 1,
  'vtk_x_s' : 1,
  'vtk_x_w' : 1,
# 'vtk_y_n': UNDEFINED,
  'vtk_y_s' : 1,
  'vtk_z_b' : 1,
# 'vtk_z_t': UNDEFINED,
  'vtp_dir' : 1,
  'vtu_dir' : 1,
# 'vw_poisson': UNDEFINED,
# 'vw_poisson_actual': UNDEFINED,
# 'wall_hamaker_constant': UNDEFINED,
# 'wall_vdw_inner_cutoff': UNDEFINED,
# 'wall_vdw_outer_cutoff': UNDEFINED,
  'write_dashboard' : 1,
  'write_vtk_files' : 1,
  'x_s0' : 1,
# 'xlength': UNDEFINED,
# 'xmin': UNDEFINED,
# 'ylength': UNDEFINED,
  'yu_standish' : 1,
# 'zlength': UNDEFINED,
#<end generated code>
}

if (__name__ == '__main__'):
    import os
    import re
    from tools.general import SCRIPT_DIRECTORY
    from tools.namelistparser import buildKeywordDoc

    doc = buildKeywordDoc(os.path.join(SCRIPT_DIRECTORY, os.pardir, 'model'))
    keys = doc.keys()
    keys.sort()

    pat = re.compile(r'\[(.*) in SI.*\]')


    unit_to_SI = { 'm':    1e-2,  # from cm
                   'kg':   1e-3,  # from g
                   'J':    4.184, # from cal
                   'kmol': 1e-3,  # from mol
                   'N':    1e-5,  # from dyne
                   'Pa':   0.1,   # from barye
                   'K':    1,
                   's':    1
    }

    # Try to get units from keyword doc. Add to predefined factors above

    for key in keys[:]:
        if key in cgs_to_SI:
            continue
        if 'gravity' in key:
            factor = 1e-2 # cm -> m
            cgs_to_SI[key] = factor
            continue
        entry = doc[key]
        if not isinstance(entry, dict):
            keys.remove(key)
            del doc[key]
            continue
        dtype = entry['dtype']
        if dtype !=  'DP':
            cgs_to_SI[key] = 1
            continue
        desc = entry['description']
        dlower = desc.lower()

        if 'fraction' in dlower:
            cgs_to_SI[key] = 1
            continue

        match = pat.search(desc)

        if match:
            SI_unit = match.group(1)
            if SI_unit.startswith('(') and SI_unit.endswith(')'):
                SI_unit = SI_unit[1:-1]
            if '/' in SI_unit:
                num, denom = SI_unit.split('/')
            else:
                num, denom = SI_unit, ''
            if num.startswith('(') and num.endswith(')'):
                num = num[1:-1]
            if denom.startswith('(') and denom.endswith(')'):
                denom = denom[1:-1]
            factor = 1.0
            for term in num.split('.'):
                if '^' in term:
                    u, pow = term.split('^')
                    pow = int(pow)
                else:
                    u = term
                    pow = 1
                factor *= unit_to_SI[u]**pow
            if denom:
                for term in denom.split('.'):
                    if '^' in term:
                        u, pow = term.split('^')
                        pow = int(pow)
                    else:
                        u = term
                        pow = 1
                    factor /= unit_to_SI[u]**pow
            cgs_to_SI[key] = round(factor,10)
    infile = open(__file__, 'r')
    outfile = open(__file__ + '.tmp', 'w')
    skip = False
    for line in infile:
        line = line.decode('utf-8')
        if line.startswith('#<start'):
            outfile.write(line)
            for key in keys:
                factor = cgs_to_SI.get(key)
                if factor is None:
                    outfile.write("# '%s': UNDEFINED,\n" % key)
                else:
                    outfile.write("  '%s' : %s,\n" % (key, factor))
            skip = True
        elif line.startswith('#<end'):
            outfile.write(line)
            skip = False
            continue
        elif not skip:
            outfile.write(line.encode('utf-8'))

    # We could replace the input file, but let's not go crazy
    print("Generated %s" % (__file__ + '.tmp'))










"""
Simulation Units
Simulations can be setup using the International System of Units (SI) or the centimetergram-second system (CGS). Although the majority of units are consistent with the specified systems, there are exceptions. The following table provides the SI and CGS units employed by MFIX for various quantities.


Quantity                      MFIX SI unit              MFIX CGS unit

length, position              meter(m)                  centimeter (cm)
mass                          kilogram (kg)             gram (g)
time                          second (s)                second (s)
thermal temperature           Kelvin (K)                Kelvin (K)
energy†                       Joule (J)                 calorie (cal)
amount of substance‡          kilomole (kmol)           mole (mol)
force                         Newton (1 N = 1 kg·m·s-2) dyne (1 dyn = 1 g·cm·s-2)
pressure                      Pascal (1 Pa = 1 N·m-2)   barye (1 Ba = 1 dyn·cm-2)
dynamic viscosity             Pa·s                      poise (1 P = 1 g·cm-1·s-1)
kinematic viscosity           m2·s-1                    Stokes (1 St = 1 cm2·s-1)
gas constant                  J·K-1·kmol-1              erg·K-1·mol-1
enthalpy                      J                         cal
specific heat                 J·kg-1·K-1                cal·g-1·K-1
thermal conductivity          J·m-1·K-1·s-1             cal·cm-1·K-1·s-1

† The CGS unit for energy is the ergon (1 erg = 1 dyne·cm). This is reflected in MFIX
through the gas constant. However, all thermochemical properties related to the energy
equations are based in calories for CGS units. Entries in the Burcat database are
always specified in terms of calories regardless of the simulation units. MFIX converts
the entries to joules after reading the database when SI units are used.
‡ The SI unit for the amount of a substance is the mole (mol). These units are needed
when specifying reaction rates:
• amount per time per volume for Eulerian model reactions
• amount per time for Lagrangian model reactions
"""
