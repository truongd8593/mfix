.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
mfix.exe : \
    adjust_a_u_g.$(OBJ_EXT) \
    adjust_a_u_s.$(OBJ_EXT) \
    adjust_a_v_g.$(OBJ_EXT) \
    adjust_a_v_s.$(OBJ_EXT) \
    adjust_a_w_g.$(OBJ_EXT) \
    adjust_a_w_s.$(OBJ_EXT) \
    adjust_dt.$(OBJ_EXT) \
    adjust_eps.$(OBJ_EXT) \
    adjust_leq.$(OBJ_EXT) \
    adjust_rop.$(OBJ_EXT) \
    adjust_theta.$(OBJ_EXT) \
    allocate_arrays.$(OBJ_EXT) \
    ambm_mod.$(OBJ_EXT) \
    b_m_p_star.$(OBJ_EXT) \
    bc_mod.$(OBJ_EXT) \
    bc_phi.$(OBJ_EXT) \
    bc_theta.$(OBJ_EXT) \
    bound_x.$(OBJ_EXT) \
    boundfunijk_mod.$(OBJ_EXT) \
    cal_d.$(OBJ_EXT) \
    calc_cell.$(OBJ_EXT) \
    calc_coeff.$(OBJ_EXT) \
    calc_d.$(OBJ_EXT) \
    calc_dif_g.$(OBJ_EXT) \
    calc_dif_s.$(OBJ_EXT) \
    calc_drag.$(OBJ_EXT) \
    calc_e.$(OBJ_EXT) \
    calc_gama.$(OBJ_EXT) \
    calc_grbdry.$(OBJ_EXT) \
    calc_k_cp.$(OBJ_EXT) \
    calc_k_g.$(OBJ_EXT) \
    calc_k_s.$(OBJ_EXT) \
    calc_mu_g.$(OBJ_EXT) \
    calc_mu_s.$(OBJ_EXT) \
    calc_mw.$(OBJ_EXT) \
    calc_outflow.$(OBJ_EXT) \
    calc_p_star.$(OBJ_EXT) \
    calc_resid.$(OBJ_EXT) \
    calc_s_ddot_s.$(OBJ_EXT) \
    calc_trd_g.$(OBJ_EXT) \
    calc_trd_s.$(OBJ_EXT) \
    calc_u_friction.$(OBJ_EXT) \
    calc_vol_fr.$(OBJ_EXT) \
    calc_xsi.$(OBJ_EXT) \
    check_ab_m.$(OBJ_EXT) \
    check_convergence.$(OBJ_EXT) \
    check_data_01.$(OBJ_EXT) \
    check_data_02.$(OBJ_EXT) \
    check_data_03.$(OBJ_EXT) \
    check_data_04.$(OBJ_EXT) \
    check_data_05.$(OBJ_EXT) \
    check_data_06.$(OBJ_EXT) \
    check_data_07.$(OBJ_EXT) \
    check_data_08.$(OBJ_EXT) \
    check_data_09.$(OBJ_EXT) \
    check_data_20.$(OBJ_EXT) \
    check_data_30.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    check_plane.$(OBJ_EXT) \
    coeff_mod.$(OBJ_EXT) \
    compare.$(OBJ_EXT) \
    constant_mod.$(OBJ_EXT) \
    cont_mod.$(OBJ_EXT) \
    conv_dif_phi.$(OBJ_EXT) \
    conv_dif_u_g.$(OBJ_EXT) \
    conv_dif_u_s.$(OBJ_EXT) \
    conv_dif_v_g.$(OBJ_EXT) \
    conv_dif_v_s.$(OBJ_EXT) \
    conv_dif_w_g.$(OBJ_EXT) \
    conv_dif_w_s.$(OBJ_EXT) \
    conv_pp_g.$(OBJ_EXT) \
    conv_rop_g.$(OBJ_EXT) \
    conv_rop_s.$(OBJ_EXT) \
    conv_source_epp.$(OBJ_EXT) \
    copy_a.$(OBJ_EXT) \
    corner.$(OBJ_EXT) \
    corner_mod.$(OBJ_EXT) \
    correct_0.$(OBJ_EXT) \
    correct_1.$(OBJ_EXT) \
    dgtsl.$(OBJ_EXT) \
    dif_u_is.$(OBJ_EXT) \
    dif_v_is.$(OBJ_EXT) \
    dif_w_is.$(OBJ_EXT) \
    discretize.$(OBJ_EXT) \
    display_resid.$(OBJ_EXT) \
    drag_gs.$(OBJ_EXT) \
    drag_mod.$(OBJ_EXT) \
    drag_ss.$(OBJ_EXT) \
    energy_mod.$(OBJ_EXT) \
    eosg.$(OBJ_EXT) \
    equal.$(OBJ_EXT) \
    error_routine.$(OBJ_EXT) \
    exchange.$(OBJ_EXT) \
    exit.$(OBJ_EXT) \
    fldvar_mod.$(OBJ_EXT) \
    flow_to_vel.$(OBJ_EXT) \
    function_mod.$(OBJ_EXT) \
    funits_mod.$(OBJ_EXT) \
    g_0.$(OBJ_EXT) \
    geometry_mod.$(OBJ_EXT) \
    get_bc_area.$(OBJ_EXT) \
    get_data.$(OBJ_EXT) \
    get_eq.$(OBJ_EXT) \
    get_flow_bc.$(OBJ_EXT) \
    get_hloss.$(OBJ_EXT) \
    get_is.$(OBJ_EXT) \
    get_philoss.$(OBJ_EXT) \
    get_smass.$(OBJ_EXT) \
    get_stats.$(OBJ_EXT) \
    get_walls_bc.$(OBJ_EXT) \
    ic_mod.$(OBJ_EXT) \
    in_bin_512.$(OBJ_EXT) \
    in_bin_512i.$(OBJ_EXT) \
    indices_mod.$(OBJ_EXT) \
    init_ab_m.$(OBJ_EXT) \
    init_fvars.$(OBJ_EXT) \
    init_namelist.$(OBJ_EXT) \
    init_resid.$(OBJ_EXT) \
    is_mod.$(OBJ_EXT) \
    iterate.$(OBJ_EXT) \
    leq_bicgs.$(OBJ_EXT) \
    leq_gmres.$(OBJ_EXT) \
    leq_sor.$(OBJ_EXT) \
    leqsol_mod.$(OBJ_EXT) \
    line_too_big.$(OBJ_EXT) \
    location.$(OBJ_EXT) \
    location_check.$(OBJ_EXT) \
    machine.$(OBJ_EXT) \
    machine_mod.$(OBJ_EXT) \
    make_upper_case.$(OBJ_EXT) \
    mark_phase_4_cor.$(OBJ_EXT) \
    matrix_mod.$(OBJ_EXT) \
    mfix.$(OBJ_EXT) \
    mod_bc_i.$(OBJ_EXT) \
    mod_bc_j.$(OBJ_EXT) \
    mod_bc_k.$(OBJ_EXT) \
    open_file.$(OBJ_EXT) \
    open_files.$(OBJ_EXT) \
    out_array.$(OBJ_EXT) \
    out_array_c.$(OBJ_EXT) \
    out_array_k.$(OBJ_EXT) \
    out_array_kc.$(OBJ_EXT) \
    out_bin_512.$(OBJ_EXT) \
    out_bin_512i.$(OBJ_EXT) \
    out_bin_512r.$(OBJ_EXT) \
    out_bin_r.$(OBJ_EXT) \
    output_mod.$(OBJ_EXT) \
    parallel_mod.$(OBJ_EXT) \
    param1_mod.$(OBJ_EXT) \
    param_mod.$(OBJ_EXT) \
    parse_line.$(OBJ_EXT) \
    parse_mod.$(OBJ_EXT) \
    parse_resid_string.$(OBJ_EXT) \
    parse_rxn.$(OBJ_EXT) \
    partial_elim.$(OBJ_EXT) \
    pgcor_mod.$(OBJ_EXT) \
    physical_prop.$(OBJ_EXT) \
    physprop_mod.$(OBJ_EXT) \
    pscor_mod.$(OBJ_EXT) \
    read_namelist.$(OBJ_EXT) \
    read_res0.$(OBJ_EXT) \
    read_res1.$(OBJ_EXT) \
    remove_comment.$(OBJ_EXT) \
    reset_new.$(OBJ_EXT) \
    residual_mod.$(OBJ_EXT) \
    rrates.$(OBJ_EXT) \
    rrates0.$(OBJ_EXT) \
    run_mod.$(OBJ_EXT) \
    rxns_mod.$(OBJ_EXT) \
    scalar_prop.$(OBJ_EXT) \
    scalars_mod.$(OBJ_EXT) \
    scales_mod.$(OBJ_EXT) \
    seek_comment.$(OBJ_EXT) \
    seek_end.$(OBJ_EXT) \
    set_bc0.$(OBJ_EXT) \
    set_bc1.$(OBJ_EXT) \
    set_constants.$(OBJ_EXT) \
    set_constprop.$(OBJ_EXT) \
    set_flags.$(OBJ_EXT) \
    set_fluidbed_p.$(OBJ_EXT) \
    set_geometry.$(OBJ_EXT) \
    set_geometry1.$(OBJ_EXT) \
    set_ic.$(OBJ_EXT) \
    set_increments.$(OBJ_EXT) \
    set_index1.$(OBJ_EXT) \
    set_index1a.$(OBJ_EXT) \
    set_l_scale.$(OBJ_EXT) \
    set_max2.$(OBJ_EXT) \
    set_mw_mix_g.$(OBJ_EXT) \
    set_outflow.$(OBJ_EXT) \
    set_ro_g.$(OBJ_EXT) \
    set_wall_bc.$(OBJ_EXT) \
    shift_dxyz.$(OBJ_EXT) \
    solve_continuity.$(OBJ_EXT) \
    solve_energy_eq.$(OBJ_EXT) \
    solve_epp.$(OBJ_EXT) \
    solve_granular_energy.$(OBJ_EXT) \
    solve_lin_eq.$(OBJ_EXT) \
    solve_pp_g.$(OBJ_EXT) \
    solve_scalar_eq.$(OBJ_EXT) \
    solve_species_eq.$(OBJ_EXT) \
    solve_vel_star.$(OBJ_EXT) \
    source_granular_energy.$(OBJ_EXT) \
    source_phi.$(OBJ_EXT) \
    source_pp_g.$(OBJ_EXT) \
    source_rop_g.$(OBJ_EXT) \
    source_rop_s.$(OBJ_EXT) \
    source_u_g.$(OBJ_EXT) \
    source_u_s.$(OBJ_EXT) \
    source_v_g.$(OBJ_EXT) \
    source_v_s.$(OBJ_EXT) \
    source_w_g.$(OBJ_EXT) \
    source_w_s.$(OBJ_EXT) \
    tau_g_mod.$(OBJ_EXT) \
    tau_s_mod.$(OBJ_EXT) \
    tau_u_g.$(OBJ_EXT) \
    tau_u_s.$(OBJ_EXT) \
    tau_v_g.$(OBJ_EXT) \
    tau_v_s.$(OBJ_EXT) \
    tau_w_g.$(OBJ_EXT) \
    tau_w_s.$(OBJ_EXT) \
    test_lin_eq.$(OBJ_EXT) \
    time_cpu_mod.$(OBJ_EXT) \
    time_march.$(OBJ_EXT) \
    tmp_array1_mod.$(OBJ_EXT) \
    tmp_array_mod.$(OBJ_EXT) \
    toleranc_mod.$(OBJ_EXT) \
    trace_mod.$(OBJ_EXT) \
    transfer.$(OBJ_EXT) \
    transport_prop.$(OBJ_EXT) \
    undef_2_0.$(OBJ_EXT) \
    under_relax.$(OBJ_EXT) \
    update_old.$(OBJ_EXT) \
    ur_facs_mod.$(OBJ_EXT) \
    usr0.$(OBJ_EXT) \
    usr1.$(OBJ_EXT) \
    usr2.$(OBJ_EXT) \
    usr3.$(OBJ_EXT) \
    usr_init_namelist.$(OBJ_EXT) \
    usr_mod.$(OBJ_EXT) \
    usr_write_out0.$(OBJ_EXT) \
    usr_write_out1.$(OBJ_EXT) \
    vavg_u_g.$(OBJ_EXT) \
    vavg_u_s.$(OBJ_EXT) \
    vavg_v_g.$(OBJ_EXT) \
    vavg_v_s.$(OBJ_EXT) \
    vavg_w_g.$(OBJ_EXT) \
    vavg_w_s.$(OBJ_EXT) \
    vf_gs_x.$(OBJ_EXT) \
    vf_gs_y.$(OBJ_EXT) \
    vf_gs_z.$(OBJ_EXT) \
    visc_g_mod.$(OBJ_EXT) \
    visc_s_mod.$(OBJ_EXT) \
    vshear_mod.$(OBJ_EXT) \
    write_ab_m.$(OBJ_EXT) \
    write_error.$(OBJ_EXT) \
    write_header.$(OBJ_EXT) \
    write_out0.$(OBJ_EXT) \
    write_out1.$(OBJ_EXT) \
    write_out3.$(OBJ_EXT) \
    write_res0.$(OBJ_EXT) \
    write_res1.$(OBJ_EXT) \
    write_spx0.$(OBJ_EXT) \
    write_spx1.$(OBJ_EXT) \
    write_table.$(OBJ_EXT) \
    write_usr0.$(OBJ_EXT) \
    write_usr1.$(OBJ_EXT) \
    xerbla.$(OBJ_EXT) \
    xsi_array_mod.$(OBJ_EXT) \
    zero_array.$(OBJ_EXT) \
    zero_norm_vel.$(OBJ_EXT) \
    compar_mod.$(OBJ_EXT) \
    dbg_util_mod.$(OBJ_EXT) \
    debug_mod.$(OBJ_EXT) \
    gridmap_mod.$(OBJ_EXT) \
    mpi_mod.$(OBJ_EXT) \
    mpi_utility_mod.$(OBJ_EXT) \
    parallel_mpi_mod.$(OBJ_EXT) \
    sendrecv_mod.$(OBJ_EXT) \
    blas90.a 
	$(LINK_CMD) $(LINK_FLAGS) \
    adjust_a_u_g.$(OBJ_EXT) \
    adjust_a_u_s.$(OBJ_EXT) \
    adjust_a_v_g.$(OBJ_EXT) \
    adjust_a_v_s.$(OBJ_EXT) \
    adjust_a_w_g.$(OBJ_EXT) \
    adjust_a_w_s.$(OBJ_EXT) \
    adjust_dt.$(OBJ_EXT) \
    adjust_eps.$(OBJ_EXT) \
    adjust_leq.$(OBJ_EXT) \
    adjust_rop.$(OBJ_EXT) \
    adjust_theta.$(OBJ_EXT) \
    allocate_arrays.$(OBJ_EXT) \
    ambm_mod.$(OBJ_EXT) \
    b_m_p_star.$(OBJ_EXT) \
    bc_mod.$(OBJ_EXT) \
    bc_phi.$(OBJ_EXT) \
    bc_theta.$(OBJ_EXT) \
    bound_x.$(OBJ_EXT) \
    boundfunijk_mod.$(OBJ_EXT) \
    cal_d.$(OBJ_EXT) \
    calc_cell.$(OBJ_EXT) \
    calc_coeff.$(OBJ_EXT) \
    calc_d.$(OBJ_EXT) \
    calc_dif_g.$(OBJ_EXT) \
    calc_dif_s.$(OBJ_EXT) \
    calc_drag.$(OBJ_EXT) \
    calc_e.$(OBJ_EXT) \
    calc_gama.$(OBJ_EXT) \
    calc_grbdry.$(OBJ_EXT) \
    calc_k_cp.$(OBJ_EXT) \
    calc_k_g.$(OBJ_EXT) \
    calc_k_s.$(OBJ_EXT) \
    calc_mu_g.$(OBJ_EXT) \
    calc_mu_s.$(OBJ_EXT) \
    calc_mw.$(OBJ_EXT) \
    calc_outflow.$(OBJ_EXT) \
    calc_p_star.$(OBJ_EXT) \
    calc_resid.$(OBJ_EXT) \
    calc_s_ddot_s.$(OBJ_EXT) \
    calc_trd_g.$(OBJ_EXT) \
    calc_trd_s.$(OBJ_EXT) \
    calc_u_friction.$(OBJ_EXT) \
    calc_vol_fr.$(OBJ_EXT) \
    calc_xsi.$(OBJ_EXT) \
    check_ab_m.$(OBJ_EXT) \
    check_convergence.$(OBJ_EXT) \
    check_data_01.$(OBJ_EXT) \
    check_data_02.$(OBJ_EXT) \
    check_data_03.$(OBJ_EXT) \
    check_data_04.$(OBJ_EXT) \
    check_data_05.$(OBJ_EXT) \
    check_data_06.$(OBJ_EXT) \
    check_data_07.$(OBJ_EXT) \
    check_data_08.$(OBJ_EXT) \
    check_data_09.$(OBJ_EXT) \
    check_data_20.$(OBJ_EXT) \
    check_data_30.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    check_plane.$(OBJ_EXT) \
    coeff_mod.$(OBJ_EXT) \
    compare.$(OBJ_EXT) \
    constant_mod.$(OBJ_EXT) \
    cont_mod.$(OBJ_EXT) \
    conv_dif_phi.$(OBJ_EXT) \
    conv_dif_u_g.$(OBJ_EXT) \
    conv_dif_u_s.$(OBJ_EXT) \
    conv_dif_v_g.$(OBJ_EXT) \
    conv_dif_v_s.$(OBJ_EXT) \
    conv_dif_w_g.$(OBJ_EXT) \
    conv_dif_w_s.$(OBJ_EXT) \
    conv_pp_g.$(OBJ_EXT) \
    conv_rop_g.$(OBJ_EXT) \
    conv_rop_s.$(OBJ_EXT) \
    conv_source_epp.$(OBJ_EXT) \
    copy_a.$(OBJ_EXT) \
    corner.$(OBJ_EXT) \
    corner_mod.$(OBJ_EXT) \
    correct_0.$(OBJ_EXT) \
    correct_1.$(OBJ_EXT) \
    dgtsl.$(OBJ_EXT) \
    dif_u_is.$(OBJ_EXT) \
    dif_v_is.$(OBJ_EXT) \
    dif_w_is.$(OBJ_EXT) \
    discretize.$(OBJ_EXT) \
    display_resid.$(OBJ_EXT) \
    drag_gs.$(OBJ_EXT) \
    drag_mod.$(OBJ_EXT) \
    drag_ss.$(OBJ_EXT) \
    energy_mod.$(OBJ_EXT) \
    eosg.$(OBJ_EXT) \
    equal.$(OBJ_EXT) \
    error_routine.$(OBJ_EXT) \
    exchange.$(OBJ_EXT) \
    exit.$(OBJ_EXT) \
    fldvar_mod.$(OBJ_EXT) \
    flow_to_vel.$(OBJ_EXT) \
    function_mod.$(OBJ_EXT) \
    funits_mod.$(OBJ_EXT) \
    g_0.$(OBJ_EXT) \
    geometry_mod.$(OBJ_EXT) \
    get_bc_area.$(OBJ_EXT) \
    get_data.$(OBJ_EXT) \
    get_eq.$(OBJ_EXT) \
    get_flow_bc.$(OBJ_EXT) \
    get_hloss.$(OBJ_EXT) \
    get_is.$(OBJ_EXT) \
    get_philoss.$(OBJ_EXT) \
    get_smass.$(OBJ_EXT) \
    get_stats.$(OBJ_EXT) \
    get_walls_bc.$(OBJ_EXT) \
    ic_mod.$(OBJ_EXT) \
    in_bin_512.$(OBJ_EXT) \
    in_bin_512i.$(OBJ_EXT) \
    indices_mod.$(OBJ_EXT) \
    init_ab_m.$(OBJ_EXT) \
    init_fvars.$(OBJ_EXT) \
    init_namelist.$(OBJ_EXT) \
    init_resid.$(OBJ_EXT) \
    is_mod.$(OBJ_EXT) \
    iterate.$(OBJ_EXT) \
    leq_bicgs.$(OBJ_EXT) \
    leq_gmres.$(OBJ_EXT) \
    leq_sor.$(OBJ_EXT) \
    leqsol_mod.$(OBJ_EXT) \
    line_too_big.$(OBJ_EXT) \
    location.$(OBJ_EXT) \
    location_check.$(OBJ_EXT) \
    machine.$(OBJ_EXT) \
    machine_mod.$(OBJ_EXT) \
    make_upper_case.$(OBJ_EXT) \
    mark_phase_4_cor.$(OBJ_EXT) \
    matrix_mod.$(OBJ_EXT) \
    mfix.$(OBJ_EXT) \
    mod_bc_i.$(OBJ_EXT) \
    mod_bc_j.$(OBJ_EXT) \
    mod_bc_k.$(OBJ_EXT) \
    open_file.$(OBJ_EXT) \
    open_files.$(OBJ_EXT) \
    out_array.$(OBJ_EXT) \
    out_array_c.$(OBJ_EXT) \
    out_array_k.$(OBJ_EXT) \
    out_array_kc.$(OBJ_EXT) \
    out_bin_512.$(OBJ_EXT) \
    out_bin_512i.$(OBJ_EXT) \
    out_bin_512r.$(OBJ_EXT) \
    out_bin_r.$(OBJ_EXT) \
    output_mod.$(OBJ_EXT) \
    parallel_mod.$(OBJ_EXT) \
    param1_mod.$(OBJ_EXT) \
    param_mod.$(OBJ_EXT) \
    parse_line.$(OBJ_EXT) \
    parse_mod.$(OBJ_EXT) \
    parse_resid_string.$(OBJ_EXT) \
    parse_rxn.$(OBJ_EXT) \
    partial_elim.$(OBJ_EXT) \
    pgcor_mod.$(OBJ_EXT) \
    physical_prop.$(OBJ_EXT) \
    physprop_mod.$(OBJ_EXT) \
    pscor_mod.$(OBJ_EXT) \
    read_namelist.$(OBJ_EXT) \
    read_res0.$(OBJ_EXT) \
    read_res1.$(OBJ_EXT) \
    remove_comment.$(OBJ_EXT) \
    reset_new.$(OBJ_EXT) \
    residual_mod.$(OBJ_EXT) \
    rrates.$(OBJ_EXT) \
    rrates0.$(OBJ_EXT) \
    run_mod.$(OBJ_EXT) \
    rxns_mod.$(OBJ_EXT) \
    scalar_prop.$(OBJ_EXT) \
    scalars_mod.$(OBJ_EXT) \
    scales_mod.$(OBJ_EXT) \
    seek_comment.$(OBJ_EXT) \
    seek_end.$(OBJ_EXT) \
    set_bc0.$(OBJ_EXT) \
    set_bc1.$(OBJ_EXT) \
    set_constants.$(OBJ_EXT) \
    set_constprop.$(OBJ_EXT) \
    set_flags.$(OBJ_EXT) \
    set_fluidbed_p.$(OBJ_EXT) \
    set_geometry.$(OBJ_EXT) \
    set_geometry1.$(OBJ_EXT) \
    set_ic.$(OBJ_EXT) \
    set_increments.$(OBJ_EXT) \
    set_index1.$(OBJ_EXT) \
    set_index1a.$(OBJ_EXT) \
    set_l_scale.$(OBJ_EXT) \
    set_max2.$(OBJ_EXT) \
    set_mw_mix_g.$(OBJ_EXT) \
    set_outflow.$(OBJ_EXT) \
    set_ro_g.$(OBJ_EXT) \
    set_wall_bc.$(OBJ_EXT) \
    shift_dxyz.$(OBJ_EXT) \
    solve_continuity.$(OBJ_EXT) \
    solve_energy_eq.$(OBJ_EXT) \
    solve_epp.$(OBJ_EXT) \
    solve_granular_energy.$(OBJ_EXT) \
    solve_lin_eq.$(OBJ_EXT) \
    solve_pp_g.$(OBJ_EXT) \
    solve_scalar_eq.$(OBJ_EXT) \
    solve_species_eq.$(OBJ_EXT) \
    solve_vel_star.$(OBJ_EXT) \
    source_granular_energy.$(OBJ_EXT) \
    source_phi.$(OBJ_EXT) \
    source_pp_g.$(OBJ_EXT) \
    source_rop_g.$(OBJ_EXT) \
    source_rop_s.$(OBJ_EXT) \
    source_u_g.$(OBJ_EXT) \
    source_u_s.$(OBJ_EXT) \
    source_v_g.$(OBJ_EXT) \
    source_v_s.$(OBJ_EXT) \
    source_w_g.$(OBJ_EXT) \
    source_w_s.$(OBJ_EXT) \
    tau_g_mod.$(OBJ_EXT) \
    tau_s_mod.$(OBJ_EXT) \
    tau_u_g.$(OBJ_EXT) \
    tau_u_s.$(OBJ_EXT) \
    tau_v_g.$(OBJ_EXT) \
    tau_v_s.$(OBJ_EXT) \
    tau_w_g.$(OBJ_EXT) \
    tau_w_s.$(OBJ_EXT) \
    test_lin_eq.$(OBJ_EXT) \
    time_cpu_mod.$(OBJ_EXT) \
    time_march.$(OBJ_EXT) \
    tmp_array1_mod.$(OBJ_EXT) \
    tmp_array_mod.$(OBJ_EXT) \
    toleranc_mod.$(OBJ_EXT) \
    trace_mod.$(OBJ_EXT) \
    transfer.$(OBJ_EXT) \
    transport_prop.$(OBJ_EXT) \
    undef_2_0.$(OBJ_EXT) \
    under_relax.$(OBJ_EXT) \
    update_old.$(OBJ_EXT) \
    ur_facs_mod.$(OBJ_EXT) \
    usr0.$(OBJ_EXT) \
    usr1.$(OBJ_EXT) \
    usr2.$(OBJ_EXT) \
    usr3.$(OBJ_EXT) \
    usr_init_namelist.$(OBJ_EXT) \
    usr_mod.$(OBJ_EXT) \
    usr_write_out0.$(OBJ_EXT) \
    usr_write_out1.$(OBJ_EXT) \
    vavg_u_g.$(OBJ_EXT) \
    vavg_u_s.$(OBJ_EXT) \
    vavg_v_g.$(OBJ_EXT) \
    vavg_v_s.$(OBJ_EXT) \
    vavg_w_g.$(OBJ_EXT) \
    vavg_w_s.$(OBJ_EXT) \
    vf_gs_x.$(OBJ_EXT) \
    vf_gs_y.$(OBJ_EXT) \
    vf_gs_z.$(OBJ_EXT) \
    visc_g_mod.$(OBJ_EXT) \
    visc_s_mod.$(OBJ_EXT) \
    vshear_mod.$(OBJ_EXT) \
    write_ab_m.$(OBJ_EXT) \
    write_error.$(OBJ_EXT) \
    write_header.$(OBJ_EXT) \
    write_out0.$(OBJ_EXT) \
    write_out1.$(OBJ_EXT) \
    write_out3.$(OBJ_EXT) \
    write_res0.$(OBJ_EXT) \
    write_res1.$(OBJ_EXT) \
    write_spx0.$(OBJ_EXT) \
    write_spx1.$(OBJ_EXT) \
    write_table.$(OBJ_EXT) \
    write_usr0.$(OBJ_EXT) \
    write_usr1.$(OBJ_EXT) \
    xerbla.$(OBJ_EXT) \
    xsi_array_mod.$(OBJ_EXT) \
    zero_array.$(OBJ_EXT) \
    zero_norm_vel.$(OBJ_EXT) \
    compar_mod.$(OBJ_EXT) \
    dbg_util_mod.$(OBJ_EXT) \
    debug_mod.$(OBJ_EXT) \
    gridmap_mod.$(OBJ_EXT) \
    mpi_mod.$(OBJ_EXT) \
    mpi_utility_mod.$(OBJ_EXT) \
    parallel_mpi_mod.$(OBJ_EXT) \
    sendrecv_mod.$(OBJ_EXT) \
  -o mfix.exe $(LIB_FLAGS)
  
blas90.a : BLAS.o
	ar cr blas90.a BLAS.o
BLAS.o : BLAS.F
	$(FORTRAN_CMD) $(FORT_FLAGS) BLAS.F
AMBM.mod : ambm_mod.f \
            PARAM.mod \
            PARAM1.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ambm_mod.f 
BC.mod : bc_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) bc_mod.f 
BOUNDFUNIJK.mod : boundfunijk_mod.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk_mod.f 
COEFF.mod : coeff_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) coeff_mod.f 
CONSTANT.mod : constant_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) constant_mod.f 
CONT.mod : cont_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) cont_mod.f 
CORNER.mod : corner_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) corner_mod.f 
DRAG.mod : drag_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) drag_mod.f 
ENERGY.mod : energy_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) energy_mod.f 
FLDVAR.mod : fldvar_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) fldvar_mod.f 
FUNCTION.mod : function_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) function_mod.f 
FUNITS.mod : funits_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) funits_mod.f 
GEOMETRY.mod : geometry_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) geometry_mod.f 
IC.mod : ic_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ic_mod.f 
INDICES.mod : indices_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) indices_mod.f 
IS.mod : is_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) is_mod.f 
LEQSOL.mod : leqsol_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) leqsol_mod.f 
MACHINE.mod : machine_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) machine_mod.f 
MATRIX.mod : matrix_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) matrix_mod.f 
OUTPUT.mod : output_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) output_mod.f 
PARALLEL.mod : parallel_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parallel_mod.f 
PARAM1.mod : param1_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) param1_mod.f 
PARAM.mod : param_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) param_mod.f 
PARSE.mod : parse_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_mod.f 
PGCOR.mod : pgcor_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pgcor_mod.f 
PHYSPROP.mod : physprop_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) physprop_mod.f 
PSCOR.mod : pscor_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pscor_mod.f 
RESIDUAL.mod : residual_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) residual_mod.f 
RUN.mod : run_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) run_mod.f 
RXNS.mod : rxns_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) rxns_mod.f 
SCALARS.mod : scalars_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scalars_mod.f 
SCALES.mod : scales_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scales_mod.f 
TAU_G.mod : tau_g_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_g_mod.f 
TAU_S.mod : tau_s_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_s_mod.f 
TIME_CPU.mod : time_cpu_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) time_cpu_mod.f 
TMP_ARRAY1.mod : tmp_array1_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array1_mod.f 
TMP_ARRAY.mod : tmp_array_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array_mod.f 
TOLERANC.mod : toleranc_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) toleranc_mod.f 
TRACE.mod : trace_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) trace_mod.f 
UR_FACS.mod : ur_facs_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ur_facs_mod.f 
USR.mod : usr_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_mod.f 
VISC_G.mod : visc_g_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_g_mod.f 
VISC_S.mod : visc_s_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_s_mod.f 
VSHEAR.mod : vshear_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) vshear_mod.f 
XSI_ARRAY.mod : xsi_array_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) xsi_array_mod.f 
COMPAR.mod : ./dmp_modules/compar_mod.f \
            MPI.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/compar_mod.f 
DBG_UTIL.mod : ./dmp_modules/dbg_util_mod.f \
            COMPAR.mod \
            GEOMETRY.mod \
            PARALLEL_MPI.mod \
            INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/dbg_util_mod.f 
DEBUG.mod : ./dmp_modules/debug_mod.f \
            DBG_UTIL.mod \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/debug_mod.f 
GRIDMAP.mod : ./dmp_modules/gridmap_mod.f \
            MPI_UTILITY.mod \
            PARALLEL_MPI.mod \
            GEOMETRY.mod \
            SENDRECV.mod \
            COMPAR.mod \
            INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/gridmap_mod.f 
MPI.mod : ./dmp_modules/mpi_mod.f \
            mpif.h                                                      
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_mod.f 
MPI_UTILITY.mod : ./dmp_modules/mpi_utility_mod.f \
            GEOMETRY.mod \
            COMPAR.mod \
            PARALLEL_MPI.mod \
            DEBUG.mod \
            INDICES.mod \
            FUNITS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_utility_mod.f 
PARALLEL_MPI.mod : ./dmp_modules/parallel_mpi_mod.f \
            GEOMETRY.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/parallel_mpi_mod.f 
SENDRECV.mod : ./dmp_modules/sendrecv_mod.f \
            PARALLEL_MPI.mod \
            DEBUG.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            MPI.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv_mod.f 
adjust_a_u_g.$(OBJ_EXT) : adjust_a_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            USR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_u_s.$(OBJ_EXT) : adjust_a_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_v_g.$(OBJ_EXT) : adjust_a_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_v_s.$(OBJ_EXT) : adjust_a_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_w_g.$(OBJ_EXT) : adjust_a_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_w_s.$(OBJ_EXT) : adjust_a_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_dt.$(OBJ_EXT) : adjust_dt.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
adjust_eps.$(OBJ_EXT) : adjust_eps.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
adjust_leq.$(OBJ_EXT) : adjust_leq.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            LEQSOL.mod 
adjust_rop.$(OBJ_EXT) : adjust_rop.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
adjust_theta.$(OBJ_EXT) : adjust_theta.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            COMPAR.mod \
            function.inc                                                
allocate_arrays.$(OBJ_EXT) : allocate_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            AMBM.mod \
            COEFF.mod \
            CONT.mod \
            DRAG.mod \
            ENERGY.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            PHYSPROP.mod \
            PSCOR.mod \
            RESIDUAL.mod \
            RXNS.mod \
            SCALARS.mod \
            TAU_G.mod \
            TAU_S.mod \
            TMP_ARRAY.mod \
            TMP_ARRAY1.mod \
            TRACE.mod \
            VISC_G.mod \
            VISC_S.mod \
            XSI_ARRAY.mod \
            VSHEAR.mod 
b_m_p_star.$(OBJ_EXT) : b_m_p_star.f \
            PARAM.mod \
            PARAM1.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            RXNS.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
bc_phi.$(OBJ_EXT) : bc_phi.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            GEOMETRY.mod \
            OUTPUT.mod \
            INDICES.mod \
            BC.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
bc_theta.$(OBJ_EXT) : bc_theta.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            GEOMETRY.mod \
            OUTPUT.mod \
            INDICES.mod \
            BC.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
bound_x.$(OBJ_EXT) : bound_x.f \
            PARAM.mod \
            PARAM1.mod 
cal_d.$(OBJ_EXT) : cal_d.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            VSHEAR.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_cell.$(OBJ_EXT) : calc_cell.f \
            PARAM.mod \
            PARAM1.mod 
calc_coeff.$(OBJ_EXT) : calc_coeff.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RXNS.mod \
            COMPAR.mod 
calc_d.$(OBJ_EXT) : calc_d.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            SCALES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_dif_g.$(OBJ_EXT) : calc_dif_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            CONSTANT.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
calc_dif_s.$(OBJ_EXT) : calc_dif_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
calc_drag.$(OBJ_EXT) : calc_drag.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            DRAG.mod \
            COMPAR.mod 
calc_e.$(OBJ_EXT) : calc_e.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            CONSTANT.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_gama.$(OBJ_EXT) : calc_gama.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            ENERGY.mod \
            RXNS.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
calc_grbdry.$(OBJ_EXT) : calc_grbdry.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_k_cp.$(OBJ_EXT) : calc_k_cp.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            PSCOR.mod \
            GEOMETRY.mod \
            CONSTANT.mod \
            RUN.mod \
            VISC_S.mod \
            TRACE.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
calc_k_g.$(OBJ_EXT) : calc_k_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            CONSTANT.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
calc_k_s.$(OBJ_EXT) : calc_k_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
calc_mu_g.$(OBJ_EXT) : calc_mu_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            VISC_G.mod \
            VISC_S.mod \
            INDICES.mod \
            CONSTANT.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                
calc_mu_s.$(OBJ_EXT) : calc_mu_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            DRAG.mod \
            RUN.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            VISC_G.mod \
            VISC_S.mod \
            TRACE.mod \
            INDICES.mod \
            CONSTANT.mod \
            VSHEAR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            s_pr1.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            s_pr2.inc                                                   
calc_mw.$(OBJ_EXT) : calc_mw.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod 
calc_outflow.$(OBJ_EXT) : calc_outflow.f \
            PARAM.mod \
            PARAM1.mod \
            BC.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
calc_p_star.$(OBJ_EXT) : calc_p_star.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            PGCOR.mod \
            PSCOR.mod \
            UR_FACS.mod \
            RESIDUAL.mod \
            COMPAR.mod \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                   
calc_resid.$(OBJ_EXT) : calc_resid.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
calc_s_ddot_s.$(OBJ_EXT) : calc_s_ddot_s.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
calc_trd_g.$(OBJ_EXT) : calc_trd_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
calc_trd_s.$(OBJ_EXT) : calc_trd_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
calc_u_friction.$(OBJ_EXT) : calc_u_friction.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            GEOMETRY.mod \
            INDICES.mod \
            BC.mod \
            COMPAR.mod \
            RUN.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_vol_fr.$(OBJ_EXT) : calc_vol_fr.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            PGCOR.mod \
            PSCOR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS2) calc_vol_fr.f 
calc_xsi.$(OBJ_EXT) : calc_xsi.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            VSHEAR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
check_ab_m.$(OBJ_EXT) : check_ab_m.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
check_convergence.$(OBJ_EXT) : check_convergence.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            RESIDUAL.mod \
            TOLERANC.mod \
            MPI_UTILITY.mod 
check_data_01.$(OBJ_EXT) : check_data_01.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            PHYSPROP.mod \
            INDICES.mod \
            SCALARS.mod 
check_data_02.$(OBJ_EXT) : check_data_02.f \
            PARAM.mod \
            PARAM1.mod \
            OUTPUT.mod \
            LEQSOL.mod \
            GEOMETRY.mod 
check_data_03.$(OBJ_EXT) : check_data_03.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
check_data_04.$(OBJ_EXT) : check_data_04.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            FUNITS.mod 
check_data_05.$(OBJ_EXT) : check_data_05.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FUNITS.mod \
            RUN.mod \
            INDICES.mod 
check_data_06.$(OBJ_EXT) : check_data_06.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            IC.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            INDICES.mod \
            FUNITS.mod \
            SCALARS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            function.inc                                                
check_data_07.$(OBJ_EXT) : check_data_07.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            BC.mod \
            INDICES.mod \
            FUNITS.mod \
            SCALARS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
check_data_08.$(OBJ_EXT) : check_data_08.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            IS.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
check_data_09.$(OBJ_EXT) : check_data_09.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            RXNS.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod 
check_data_20.$(OBJ_EXT) : check_data_20.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            VISC_G.mod \
            RXNS.mod \
            SCALARS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
check_data_30.$(OBJ_EXT) : check_data_30.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            RXNS.mod \
            VISC_S.mod \
            VISC_G.mod \
            GEOMETRY.mod \
            RUN.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
check_one_axis.$(OBJ_EXT) : check_one_axis.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod 
check_plane.$(OBJ_EXT) : check_plane.f \
            FUNITS.mod \
            COMPAR.mod 
compare.$(OBJ_EXT) : compare.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
conv_dif_phi.$(OBJ_EXT) : conv_dif_phi.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            SENDRECV.mod \
            XSI_ARRAY.mod \
            MPI_UTILITY.mod \
            INDICES.mod \
            PARALLEL.mod \
            MATRIX.mod \
            TOLERANC.mod \
            VSHEAR.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            OUTPUT.mod \
            IS.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
conv_dif_u_g.$(OBJ_EXT) : conv_dif_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            VISC_G.mod \
            COMPAR.mod \
            TOLERANC.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            VSHEAR.mod \
            XSI_ARRAY.mod \
            TMP_ARRAY.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_dif_u_s.$(OBJ_EXT) : conv_dif_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            PHYSPROP.mod \
            VISC_S.mod \
            COMPAR.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            XSI_ARRAY.mod \
            TMP_ARRAY.mod \
            VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_dif_v_g.$(OBJ_EXT) : conv_dif_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            VISC_G.mod \
            COMPAR.mod \
            TOLERANC.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            XSI_ARRAY.mod \
            VSHEAR.mod \
            TMP_ARRAY.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_dif_v_s.$(OBJ_EXT) : conv_dif_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            PHYSPROP.mod \
            VISC_S.mod \
            COMPAR.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            XSI_ARRAY.mod \
            TMP_ARRAY.mod \
            VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_dif_w_g.$(OBJ_EXT) : conv_dif_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            VISC_G.mod \
            COMPAR.mod \
            TOLERANC.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            XSI_ARRAY.mod \
            TMP_ARRAY.mod \
            VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_dif_w_s.$(OBJ_EXT) : conv_dif_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            PHYSPROP.mod \
            VISC_S.mod \
            COMPAR.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            XSI_ARRAY.mod \
            TMP_ARRAY.mod \
            VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_pp_g.$(OBJ_EXT) : conv_pp_g.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            PARALLEL.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            COMPAR.mod \
            XSI_ARRAY.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_rop_g.$(OBJ_EXT) : conv_rop_g.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            COMPAR.mod \
            PARALLEL.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            XSI_ARRAY.mod \
            function.inc                                                
conv_rop_s.$(OBJ_EXT) : conv_rop_s.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            COMPAR.mod \
            PARALLEL.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            PSCOR.mod \
            XSI_ARRAY.mod \
            function.inc                                                
conv_source_epp.$(OBJ_EXT) : conv_source_epp.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            SENDRECV.mod \
            XSI_ARRAY.mod \
            PARALLEL.mod \
            MATRIX.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            RXNS.mod \
            INDICES.mod \
            PGCOR.mod \
            PSCOR.mod \
            VSHEAR.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
copy_a.$(OBJ_EXT) : copy_a.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            function.inc                                                
corner.$(OBJ_EXT) : corner.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            MATRIX.mod \
            CORNER.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
correct_0.$(OBJ_EXT) : correct_0.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PGCOR.mod \
            UR_FACS.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            function.inc                                                
correct_1.$(OBJ_EXT) : correct_1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            GEOMETRY.mod \
            PSCOR.mod \
            UR_FACS.mod \
            CONSTANT.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
dgtsl.$(OBJ_EXT) : dgtsl.f 
dif_u_is.$(OBJ_EXT) : dif_u_is.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            OUTPUT.mod \
            INDICES.mod \
            IS.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
dif_v_is.$(OBJ_EXT) : dif_v_is.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            OUTPUT.mod \
            INDICES.mod \
            IS.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
dif_w_is.$(OBJ_EXT) : dif_w_is.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            OUTPUT.mod \
            INDICES.mod \
            IS.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
discretize.$(OBJ_EXT) : discretize.f \
            PARAM.mod \
            PARAM1.mod 
display_resid.$(OBJ_EXT) : display_resid.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RESIDUAL.mod \
            FLDVAR.mod \
            COMPAR.mod 
drag_gs.$(OBJ_EXT) : drag_gs.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            CONSTANT.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
drag_ss.$(OBJ_EXT) : drag_ss.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
eosg.$(OBJ_EXT) : eosg.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            SCALES.mod \
            sc_p_g1.inc                                                  \
            sc_p_g2.inc                                                 
equal.$(OBJ_EXT) : equal.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            PHYSPROP.mod 
error_routine.$(OBJ_EXT) : error_routine.f \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
exchange.$(OBJ_EXT) : exchange.f \
            PARAM.mod \
            PARAM1.mod \
            COMPAR.mod 
exit.$(OBJ_EXT) : exit.f \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
flow_to_vel.$(OBJ_EXT) : flow_to_vel.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            BC.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod 
g_0.$(OBJ_EXT) : g_0.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
get_bc_area.$(OBJ_EXT) : get_bc_area.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            BC.mod \
            COMPAR.mod 
get_data.$(OBJ_EXT) : get_data.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            FUNITS.mod \
            COMPAR.mod \
            GRIDMAP.mod 
get_eq.$(OBJ_EXT) : get_eq.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod 
get_flow_bc.$(OBJ_EXT) : get_flow_bc.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            BC.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
get_hloss.$(OBJ_EXT) : get_hloss.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            BC.mod \
            INDICES.mod \
            ENERGY.mod 
get_is.$(OBJ_EXT) : get_is.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            IS.mod \
            INDICES.mod \
            FUNITS.mod 
get_philoss.$(OBJ_EXT) : get_philoss.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            BC.mod \
            INDICES.mod \
            ENERGY.mod \
            COMPAR.mod \
            function.inc                                                
get_smass.$(OBJ_EXT) : get_smass.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
get_stats.$(OBJ_EXT) : get_stats.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            FUNITS.mod \
            RESIDUAL.mod \
            RUN.mod \
            COMPAR.mod \
            function.inc                                                
get_walls_bc.$(OBJ_EXT) : get_walls_bc.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            BC.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
in_bin_512.$(OBJ_EXT) : in_bin_512.f \
            MACHINE.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
in_bin_512i.$(OBJ_EXT) : in_bin_512i.f \
            MACHINE.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
init_ab_m.$(OBJ_EXT) : init_ab_m.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            PARALLEL.mod \
            COMPAR.mod 
init_fvars.$(OBJ_EXT) : init_fvars.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            SCALARS.mod \
            COMPAR.mod 
init_namelist.$(OBJ_EXT) : init_namelist.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            IC.mod \
            BC.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            IS.mod \
            TOLERANC.mod \
            SCALES.mod \
            UR_FACS.mod \
            LEQSOL.mod \
            RESIDUAL.mod \
            RXNS.mod \
            SCALARS.mod \
            COMPAR.mod \
            namelist.inc                                                
init_resid.$(OBJ_EXT) : init_resid.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RESIDUAL.mod 
iterate.$(OBJ_EXT) : iterate.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            FUNITS.mod \
            TIME_CPU.mod \
            PSCOR.mod \
            COEFF.mod \
            LEQSOL.mod \
            VISC_G.mod \
            PGCOR.mod \
            CONT.mod \
            SCALARS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
leq_bicgs.$(OBJ_EXT) : leq_bicgs.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            PARALLEL.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            FUNITS.mod \
            function.inc                                                
leq_gmres.$(OBJ_EXT) : leq_gmres.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            DEBUG.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PARALLEL.mod \
            FUNITS.mod \
            GRIDMAP.mod \
            function.inc                                                
leq_sor.$(OBJ_EXT) : leq_sor.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
line_too_big.$(OBJ_EXT) : line_too_big.f 
location.$(OBJ_EXT) : location.f \
            PARAM.mod \
            PARAM1.mod 
location_check.$(OBJ_EXT) : location_check.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod \
            GEOMETRY.mod 
machine.$(OBJ_EXT) : machine.f \
            MACHINE.mod \
            PARAM.mod \
            RUN.mod \
            COMPAR.mod \
            FUNITS.mod 
make_upper_case.$(OBJ_EXT) : make_upper_case.f 
mark_phase_4_cor.$(OBJ_EXT) : mark_phase_4_cor.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS2) mark_phase_4_cor.f 
mfix.$(OBJ_EXT) : mfix.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            TIME_CPU.mod \
            FUNITS.mod \
            OUTPUT.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PARALLEL_MPI.mod \
            function.inc                                                
mod_bc_i.$(OBJ_EXT) : mod_bc_i.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
mod_bc_j.$(OBJ_EXT) : mod_bc_j.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
mod_bc_k.$(OBJ_EXT) : mod_bc_k.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
open_file.$(OBJ_EXT) : open_file.f \
            COMPAR.mod 
open_files.$(OBJ_EXT) : open_files.f \
            MACHINE.mod \
            FUNITS.mod \
            COMPAR.mod 
out_array.$(OBJ_EXT) : out_array.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
out_array_c.$(OBJ_EXT) : out_array_c.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
out_array_k.$(OBJ_EXT) : out_array_k.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
out_array_kc.$(OBJ_EXT) : out_array_kc.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
out_bin_512.$(OBJ_EXT) : out_bin_512.f \
            MACHINE.mod 
out_bin_512i.$(OBJ_EXT) : out_bin_512i.f \
            MACHINE.mod 
out_bin_512r.$(OBJ_EXT) : out_bin_512r.f \
            MACHINE.mod 
out_bin_r.$(OBJ_EXT) : out_bin_r.f \
            PARAM.mod 
parse_line.$(OBJ_EXT) : parse_line.f \
            PARAM.mod \
            PARAM1.mod \
            PARSE.mod \
            COMPAR.mod 
parse_resid_string.$(OBJ_EXT) : parse_resid_string.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RESIDUAL.mod \
            FUNITS.mod \
            COMPAR.mod 
parse_rxn.$(OBJ_EXT) : parse_rxn.f \
            PARAM.mod \
            PARAM1.mod \
            PARSE.mod \
            RXNS.mod 
partial_elim.$(OBJ_EXT) : partial_elim.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            RUN.mod \
            function.inc                                                
physical_prop.$(OBJ_EXT) : physical_prop.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            TOLERANC.mod \
            CONSTANT.mod \
            COMPAR.mod \
            FUNITS.mod \
            cp_fun1.inc                                                  \
            function.inc                                                 \
            cp_fun2.inc                                                 
read_namelist.$(OBJ_EXT) : read_namelist.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            IC.mod \
            IS.mod \
            BC.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            TOLERANC.mod \
            FUNITS.mod \
            SCALES.mod \
            UR_FACS.mod \
            LEQSOL.mod \
            RESIDUAL.mod \
            RXNS.mod \
            SCALARS.mod \
            COMPAR.mod \
            usrnlst.inc                                                  \
            namelist.inc                                                
read_res0.$(OBJ_EXT) : read_res0.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            IC.mod \
            BC.mod \
            IS.mod \
            CONSTANT.mod \
            FUNITS.mod \
            OUTPUT.mod \
            SCALES.mod \
            UR_FACS.mod \
            TOLERANC.mod \
            LEQSOL.mod \
            SCALARS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
read_res1.$(OBJ_EXT) : read_res1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            ENERGY.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod 
remove_comment.$(OBJ_EXT) : remove_comment.f 
reset_new.$(OBJ_EXT) : reset_new.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            SCALARS.mod 
rrates.$(OBJ_EXT) : rrates.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RXNS.mod \
            ENERGY.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
rrates0.$(OBJ_EXT) : rrates0.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RXNS.mod \
            ENERGY.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            FUNITS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
scalar_prop.$(OBJ_EXT) : scalar_prop.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            SCALARS.mod \
            TOLERANC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
seek_comment.$(OBJ_EXT) : seek_comment.f 
seek_end.$(OBJ_EXT) : seek_end.f 
set_bc0.$(OBJ_EXT) : set_bc0.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PHYSPROP.mod \
            BC.mod \
            FLDVAR.mod \
            INDICES.mod \
            RUN.mod \
            FUNITS.mod \
            SCALES.mod \
            SCALARS.mod \
            BOUNDFUNIJK.mod \
            sc_p_g1.inc                                                  \
            function.inc                                                 \
            sc_p_g2.inc                                                 
set_bc1.$(OBJ_EXT) : set_bc1.f \
            PARAM.mod \
            PARAM1.mod \
            BC.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
set_constants.$(OBJ_EXT) : set_constants.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            VISC_S.mod \
            ENERGY.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            RUN.mod \
            FUNITS.mod \
            DRAG.mod 
set_constprop.$(OBJ_EXT) : set_constprop.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            VISC_S.mod \
            VISC_G.mod \
            ENERGY.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            RUN.mod \
            FUNITS.mod \
            DRAG.mod \
            COMPAR.mod \
            function.inc                                                
set_flags.$(OBJ_EXT) : set_flags.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            BC.mod \
            IS.mod \
            INDICES.mod \
            PHYSPROP.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
set_fluidbed_p.$(OBJ_EXT) : set_fluidbed_p.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            BC.mod \
            IC.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            FUNITS.mod \
            SCALES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            sc_p_g1.inc                                                  \
            b_force1.inc                                                 \
            function.inc                                                 \
            b_force2.inc                                                 \
            sc_p_g2.inc                                                 
set_geometry.$(OBJ_EXT) : set_geometry.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            GEOMETRY.mod \
            COMPAR.mod 
set_geometry1.$(OBJ_EXT) : set_geometry1.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
set_ic.$(OBJ_EXT) : set_ic.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            IC.mod \
            FLDVAR.mod \
            VISC_G.mod \
            INDICES.mod \
            SCALES.mod \
            ENERGY.mod \
            SCALARS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            sc_p_g1.inc                                                  \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            sc_p_g2.inc                                                 
set_increments.$(OBJ_EXT) : set_increments.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            FUNITS.mod \
            function.inc                                                
set_index1.$(OBJ_EXT) : set_index1.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            CONSTANT.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
set_index1a.$(OBJ_EXT) : set_index1a.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            BOUNDFUNIJK.mod \
            function.inc                                                
set_l_scale.$(OBJ_EXT) : set_l_scale.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            VISC_G.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod 
set_max2.$(OBJ_EXT) : set_max2.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod 
set_mw_mix_g.$(OBJ_EXT) : set_mw_mix_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
set_outflow.$(OBJ_EXT) : set_outflow.f \
            PARAM.mod \
            PARAM1.mod \
            BC.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            SCALARS.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
set_ro_g.$(OBJ_EXT) : set_ro_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
set_wall_bc.$(OBJ_EXT) : set_wall_bc.f \
            PARAM.mod \
            PARAM1.mod \
            BC.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
shift_dxyz.$(OBJ_EXT) : shift_dxyz.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod 
solve_continuity.$(OBJ_EXT) : solve_continuity.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            RESIDUAL.mod \
            CONT.mod \
            LEQSOL.mod \
            AMBM.mod 
solve_energy_eq.$(OBJ_EXT) : solve_energy_eq.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            BC.mod \
            ENERGY.mod \
            RXNS.mod \
            AMBM.mod \
            TMP_ARRAY.mod \
            TMP_ARRAY1.mod \
            COMPAR.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
solve_epp.$(OBJ_EXT) : solve_epp.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PSCOR.mod \
            RESIDUAL.mod \
            LEQSOL.mod \
            PHYSPROP.mod \
            AMBM.mod 
solve_granular_energy.$(OBJ_EXT) : solve_granular_energy.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            BC.mod \
            ENERGY.mod \
            RXNS.mod \
            AMBM.mod \
            TMP_ARRAY.mod \
            COMPAR.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
solve_lin_eq.$(OBJ_EXT) : solve_lin_eq.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod 
solve_pp_g.$(OBJ_EXT) : solve_pp_g.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            PGCOR.mod \
            RESIDUAL.mod \
            LEQSOL.mod \
            RUN.mod \
            AMBM.mod 
solve_scalar_eq.$(OBJ_EXT) : solve_scalar_eq.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            BC.mod \
            ENERGY.mod \
            RXNS.mod \
            SCALARS.mod \
            AMBM.mod \
            TMP_ARRAY.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
solve_species_eq.$(OBJ_EXT) : solve_species_eq.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            BC.mod \
            ENERGY.mod \
            RXNS.mod \
            AMBM.mod \
            TMP_ARRAY.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
solve_vel_star.$(OBJ_EXT) : solve_vel_star.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            AMBM.mod \
            TMP_ARRAY1.mod \
            COMPAR.mod 
source_granular_energy.$(OBJ_EXT) : source_granular_energy.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            DRAG.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            VISC_G.mod \
            VISC_S.mod \
            TRACE.mod \
            INDICES.mod \
            CONSTANT.mod \
            COMPAR.mod \
            s_pr1.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            s_pr2.inc                                                   
source_phi.$(OBJ_EXT) : source_phi.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
source_pp_g.$(OBJ_EXT) : source_pp_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RXNS.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            BC.mod \
            VSHEAR.mod \
            XSI_ARRAY.mod \
            COMPAR.mod \
            function.inc                                                
source_rop_g.$(OBJ_EXT) : source_rop_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            RXNS.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            COMPAR.mod \
            function.inc                                                
source_rop_s.$(OBJ_EXT) : source_rop_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            RXNS.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            PSCOR.mod \
            COMPAR.mod \
            function.inc                                                
source_u_g.$(OBJ_EXT) : source_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_u_s.$(OBJ_EXT) : source_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_v_g.$(OBJ_EXT) : source_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            VSHEAR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_v_s.$(OBJ_EXT) : source_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            VSHEAR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_w_g.$(OBJ_EXT) : source_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_w_s.$(OBJ_EXT) : source_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
tau_u_g.$(OBJ_EXT) : tau_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_u_s.$(OBJ_EXT) : tau_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            VSHEAR.mod \
            SENDRECV.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_v_g.$(OBJ_EXT) : tau_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            SENDRECV.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_v_s.$(OBJ_EXT) : tau_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            SENDRECV.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_w_g.$(OBJ_EXT) : tau_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            SENDRECV.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_w_s.$(OBJ_EXT) : tau_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            SENDRECV.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
test_lin_eq.$(OBJ_EXT) : test_lin_eq.f \
            PARAM.mod \
            PARAM1.mod 
time_march.$(OBJ_EXT) : time_march.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PGCOR.mod \
            PSCOR.mod \
            CONT.mod \
            COEFF.mod \
            TAU_G.mod \
            TAU_S.mod \
            VISC_G.mod \
            VISC_S.mod \
            FUNITS.mod \
            VSHEAR.mod \
            SCALARS.mod \
            DRAG.mod \
            COMPAR.mod 
transfer.$(OBJ_EXT) : transfer.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod 
transport_prop.$(OBJ_EXT) : transport_prop.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            TOLERANC.mod \
            COMPAR.mod 
undef_2_0.$(OBJ_EXT) : undef_2_0.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod 
under_relax.$(OBJ_EXT) : under_relax.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
update_old.$(OBJ_EXT) : update_old.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            TRACE.mod \
            VISC_S.mod \
            SCALARS.mod 
usr0.$(OBJ_EXT) : usr0.f \
            USR.mod 
usr1.$(OBJ_EXT) : usr1.f \
            USR.mod 
usr2.$(OBJ_EXT) : usr2.f \
            USR.mod 
usr3.$(OBJ_EXT) : usr3.f \
            USR.mod 
usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f \
            usrnlst.inc                                                 
usr_write_out0.$(OBJ_EXT) : usr_write_out0.f 
usr_write_out1.$(OBJ_EXT) : usr_write_out1.f 
vavg_u_g.$(OBJ_EXT) : vavg_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
vavg_u_s.$(OBJ_EXT) : vavg_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vavg_v_g.$(OBJ_EXT) : vavg_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
vavg_v_s.$(OBJ_EXT) : vavg_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vavg_w_g.$(OBJ_EXT) : vavg_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
vavg_w_s.$(OBJ_EXT) : vavg_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vf_gs_x.$(OBJ_EXT) : vf_gs_x.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
vf_gs_y.$(OBJ_EXT) : vf_gs_y.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
vf_gs_z.$(OBJ_EXT) : vf_gs_z.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
write_ab_m.$(OBJ_EXT) : write_ab_m.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            INDICES.mod \
            function.inc                                                
write_error.$(OBJ_EXT) : write_error.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod 
write_header.$(OBJ_EXT) : write_header.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            FUNITS.mod \
            COMPAR.mod 
write_out0.$(OBJ_EXT) : write_out0.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            IC.mod \
            BC.mod \
            IS.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            FUNITS.mod \
            TOLERANC.mod \
            SCALES.mod \
            SCALARS.mod \
            UR_FACS.mod \
            LEQSOL.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            function.inc                                                
write_out1.$(OBJ_EXT) : write_out1.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
write_out3.$(OBJ_EXT) : write_out3.f \
            FUNITS.mod \
            COMPAR.mod 
write_res0.$(OBJ_EXT) : write_res0.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            IC.mod \
            IS.mod \
            BC.mod \
            CONSTANT.mod \
            FUNITS.mod \
            OUTPUT.mod \
            SCALES.mod \
            SCALARS.mod \
            UR_FACS.mod \
            LEQSOL.mod \
            TOLERANC.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod 
write_res1.$(OBJ_EXT) : write_res1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            OUTPUT.mod \
            ENERGY.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod 
write_spx0.$(OBJ_EXT) : write_spx0.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
write_spx1.$(OBJ_EXT) : write_spx1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            FUNITS.mod \
            SCALARS.mod \
            OUTPUT.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod 
write_table.$(OBJ_EXT) : write_table.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod 
write_usr0.$(OBJ_EXT) : write_usr0.f 
write_usr1.$(OBJ_EXT) : write_usr1.f 
xerbla.$(OBJ_EXT) : xerbla.f \
            COMPAR.mod 
zero_array.$(OBJ_EXT) : zero_array.f \
            PARAM.mod \
            PARAM1.mod 
zero_norm_vel.$(OBJ_EXT) : zero_norm_vel.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            INDICES.mod \
            IS.mod \
            COMPAR.mod \
            function.inc                                                
