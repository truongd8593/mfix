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
    dsdgmr_a_mod.$(OBJ_EXT) \
    dslucs_a_mod.$(OBJ_EXT) \
    dslugm_a_mod.$(OBJ_EXT) \
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
    igcg_a_mod.$(OBJ_EXT) \
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
    leq_dsdgmr.$(OBJ_EXT) \
    leq_dslucs.$(OBJ_EXT) \
    leq_dslugm.$(OBJ_EXT) \
    leq_gmres.$(OBJ_EXT) \
    leq_igcg.$(OBJ_EXT) \
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
    xsi_array_mod.$(OBJ_EXT) \
    zero_array.$(OBJ_EXT) \
    zero_norm_vel.$(OBJ_EXT) \
    compar_mod.$(OBJ_EXT) \
    debug_mod.$(OBJ_EXT) \
    gridmap_mod.$(OBJ_EXT) \
    mpi_mod.$(OBJ_EXT) \
    mpi_utility_mod.$(OBJ_EXT) \
    parallel_mpi_mod.$(OBJ_EXT) \
    sendrecv_mod.$(OBJ_EXT) \
    igcg.$(OBJ_EXT) \
    igcg_i_mod.$(OBJ_EXT) \
    dcgs.$(OBJ_EXT) \
    dgmres.$(OBJ_EXT) \
    dlaputil.$(OBJ_EXT) \
    dmset.$(OBJ_EXT) \
    dmvops.$(OBJ_EXT) \
    mach.$(OBJ_EXT) \
    solblk_mod.$(OBJ_EXT) \
    xersla.$(OBJ_EXT) \
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
    dsdgmr_a_mod.$(OBJ_EXT) \
    dslucs_a_mod.$(OBJ_EXT) \
    dslugm_a_mod.$(OBJ_EXT) \
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
    igcg_a_mod.$(OBJ_EXT) \
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
    leq_dsdgmr.$(OBJ_EXT) \
    leq_dslucs.$(OBJ_EXT) \
    leq_dslugm.$(OBJ_EXT) \
    leq_gmres.$(OBJ_EXT) \
    leq_igcg.$(OBJ_EXT) \
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
    xsi_array_mod.$(OBJ_EXT) \
    zero_array.$(OBJ_EXT) \
    zero_norm_vel.$(OBJ_EXT) \
    compar_mod.$(OBJ_EXT) \
    debug_mod.$(OBJ_EXT) \
    gridmap_mod.$(OBJ_EXT) \
    mpi_mod.$(OBJ_EXT) \
    mpi_utility_mod.$(OBJ_EXT) \
    parallel_mpi_mod.$(OBJ_EXT) \
    sendrecv_mod.$(OBJ_EXT) \
    igcg.$(OBJ_EXT) \
    igcg_i_mod.$(OBJ_EXT) \
    dcgs.$(OBJ_EXT) \
    dgmres.$(OBJ_EXT) \
    dlaputil.$(OBJ_EXT) \
    dmset.$(OBJ_EXT) \
    dmvops.$(OBJ_EXT) \
    mach.$(OBJ_EXT) \
    solblk_mod.$(OBJ_EXT) \
    xersla.$(OBJ_EXT) \
  -o mfix.exe $(LIB_FLAGS)
  
blas90.a : BLAS.o
	ar cr blas90.a BLAS.o
BLAS.o : BLAS.F
	$(FORTRAN_CMD) $(FORT_FLAGS) BLAS.F
ambm.mod : ambm_mod.f \
            param.mod \
            param1.mod \
            compar.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ambm_mod.f 
bc.mod : bc_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) bc_mod.f 
boundfunijk.mod : boundfunijk_mod.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk_mod.f 
coeff.mod : coeff_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) coeff_mod.f 
constant.mod : constant_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) constant_mod.f 
cont.mod : cont_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) cont_mod.f 
corner.mod : corner_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) corner_mod.f 
drag.mod : drag_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) drag_mod.f 
dsdgmr_a.mod : dsdgmr_a_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) dsdgmr_a_mod.f 
dslucs_a.mod : dslucs_a_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) dslucs_a_mod.f 
dslugm_a.mod : dslugm_a_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) dslugm_a_mod.f 
energy.mod : energy_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) energy_mod.f 
fldvar.mod : fldvar_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) fldvar_mod.f 
function.mod : function_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) function_mod.f 
funits.mod : funits_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) funits_mod.f 
geometry.mod : geometry_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) geometry_mod.f 
ic.mod : ic_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ic_mod.f 
igcg_a.mod : igcg_a_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) igcg_a_mod.f 
indices.mod : indices_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) indices_mod.f 
is.mod : is_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) is_mod.f 
leqsol.mod : leqsol_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) leqsol_mod.f 
machine.mod : machine_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) machine_mod.f 
matrix.mod : matrix_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) matrix_mod.f 
output.mod : output_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) output_mod.f 
parallel.mod : parallel_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parallel_mod.f 
param1.mod : param1_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) param1_mod.f 
param.mod : param_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) param_mod.f 
parse.mod : parse_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_mod.f 
pgcor.mod : pgcor_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pgcor_mod.f 
physprop.mod : physprop_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) physprop_mod.f 
pscor.mod : pscor_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pscor_mod.f 
residual.mod : residual_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) residual_mod.f 
run.mod : run_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) run_mod.f 
rxns.mod : rxns_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) rxns_mod.f 
scales.mod : scales_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scales_mod.f 
tau_g.mod : tau_g_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_g_mod.f 
tau_s.mod : tau_s_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_s_mod.f 
time_cpu.mod : time_cpu_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) time_cpu_mod.f 
tmp_array1.mod : tmp_array1_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array1_mod.f 
tmp_array.mod : tmp_array_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array_mod.f 
toleranc.mod : toleranc_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) toleranc_mod.f 
trace.mod : trace_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) trace_mod.f 
ur_facs.mod : ur_facs_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ur_facs_mod.f 
usr.mod : usr_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_mod.f 
visc_g.mod : visc_g_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_g_mod.f 
visc_s.mod : visc_s_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_s_mod.f 
xsi_array.mod : xsi_array_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) xsi_array_mod.f 
compar.mod : ./dmp_modules/compar_mod.f \
            mpi.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/compar_mod.f 
debug.mod : ./dmp_modules/debug_mod.f \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/debug_mod.f 
gridmap.mod : ./dmp_modules/gridmap_mod.f \
            parallel_mpi.mod \
            geometry.mod \
            sendrecv.mod \
            compar.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/gridmap_mod.f 
mpi.mod : ./dmp_modules/mpi_mod.f \
            /usr/include/mpif.h                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_mod.f 
mpi_utility.mod : ./dmp_modules/mpi_utility_mod.f \
            geometry.mod \
            compar.mod \
            parallel_mpi.mod \
            debug.mod \
            indices.mod \
            funits.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_utility_mod.f 
parallel_mpi.mod : ./dmp_modules/parallel_mpi_mod.f \
            geometry.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/parallel_mpi_mod.f 
sendrecv.mod : ./dmp_modules/sendrecv_mod.f \
            parallel_mpi.mod \
            debug.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            mpi.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv_mod.f 
igcg_i.mod : ./leq_igcg/igcg_i_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./leq_igcg/igcg_i_mod.f 
solblk.mod : ./leq_slap/solblk_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./leq_slap/solblk_mod.f 
adjust_a_u_g.$(OBJ_EXT) : adjust_a_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_u_s.$(OBJ_EXT) : adjust_a_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_v_g.$(OBJ_EXT) : adjust_a_v_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_v_s.$(OBJ_EXT) : adjust_a_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_w_g.$(OBJ_EXT) : adjust_a_w_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_w_s.$(OBJ_EXT) : adjust_a_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_dt.$(OBJ_EXT) : adjust_dt.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            compar.mod \
            mpi_utility.mod 
adjust_eps.$(OBJ_EXT) : adjust_eps.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            compar.mod \
            function.inc                                                
adjust_leq.$(OBJ_EXT) : adjust_leq.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            leqsol.mod 
adjust_rop.$(OBJ_EXT) : adjust_rop.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
adjust_theta.$(OBJ_EXT) : adjust_theta.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            compar.mod \
            function.inc                                                
allocate_arrays.$(OBJ_EXT) : allocate_arrays.f \
            param.mod \
            param1.mod \
            ambm.mod \
            coeff.mod \
            cont.mod \
            drag.mod \
            dsdgmr_a.mod \
            energy.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            physprop.mod \
            pscor.mod \
            residual.mod \
            rxns.mod \
            tau_g.mod \
            tau_s.mod \
            tmp_array.mod \
            tmp_array1.mod \
            trace.mod \
            visc_g.mod \
            visc_s.mod \
            xsi_array.mod \
            dslucs_a.mod \
            dslugm_a.mod \
            igcg_i.mod \
            igcg_a.mod 
b_m_p_star.$(OBJ_EXT) : b_m_p_star.f \
            param.mod \
            param1.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            rxns.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
bc_phi.$(OBJ_EXT) : bc_phi.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            geometry.mod \
            output.mod \
            indices.mod \
            bc.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
bc_theta.$(OBJ_EXT) : bc_theta.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            geometry.mod \
            output.mod \
            indices.mod \
            bc.mod \
            compar.mod \
            mpi_utility.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
bound_x.$(OBJ_EXT) : bound_x.f \
            param.mod \
            param1.mod 
calc_cell.$(OBJ_EXT) : calc_cell.f \
            param.mod \
            param1.mod 
calc_coeff.$(OBJ_EXT) : calc_coeff.f \
            param.mod \
            param1.mod \
            physprop.mod \
            rxns.mod 
calc_d.$(OBJ_EXT) : calc_d.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            scales.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_dif_g.$(OBJ_EXT) : calc_dif_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            constant.mod \
            compar.mod \
            function.inc                                                
calc_dif_s.$(OBJ_EXT) : calc_dif_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            constant.mod \
            toleranc.mod \
            compar.mod \
            function.inc                                                
calc_drag.$(OBJ_EXT) : calc_drag.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            drag.mod 
calc_e.$(OBJ_EXT) : calc_e.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            constant.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_gama.$(OBJ_EXT) : calc_gama.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            energy.mod \
            rxns.mod \
            indices.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
calc_grbdry.$(OBJ_EXT) : calc_grbdry.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_k_cp.$(OBJ_EXT) : calc_k_cp.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            pscor.mod \
            geometry.mod \
            constant.mod \
            run.mod \
            visc_s.mod \
            trace.mod \
            compar.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
calc_k_g.$(OBJ_EXT) : calc_k_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            constant.mod \
            compar.mod \
            function.inc                                                
calc_k_s.$(OBJ_EXT) : calc_k_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            constant.mod \
            toleranc.mod \
            compar.mod \
            function.inc                                                
calc_mu_g.$(OBJ_EXT) : calc_mu_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            visc_g.mod \
            visc_s.mod \
            indices.mod \
            constant.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                
calc_mu_s.$(OBJ_EXT) : calc_mu_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            drag.mod \
            run.mod \
            geometry.mod \
            fldvar.mod \
            visc_g.mod \
            visc_s.mod \
            trace.mod \
            indices.mod \
            constant.mod \
            compar.mod \
            s_pr1.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            s_pr2.inc                                                   
calc_mw.$(OBJ_EXT) : calc_mw.f \
            param.mod \
            param1.mod \
            toleranc.mod 
calc_outflow.$(OBJ_EXT) : calc_outflow.f \
            param.mod \
            param1.mod \
            bc.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
calc_p_star.$(OBJ_EXT) : calc_p_star.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            pgcor.mod \
            pscor.mod \
            ur_facs.mod \
            residual.mod \
            compar.mod \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                   
calc_resid.$(OBJ_EXT) : calc_resid.f \
            param.mod \
            param1.mod \
            matrix.mod \
            parallel.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
calc_s_ddot_s.$(OBJ_EXT) : calc_s_ddot_s.f \
            param.mod \
            param1.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
calc_trd_g.$(OBJ_EXT) : calc_trd_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
calc_trd_s.$(OBJ_EXT) : calc_trd_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            function.inc                                                
calc_u_friction.$(OBJ_EXT) : calc_u_friction.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            geometry.mod \
            indices.mod \
            bc.mod \
            compar.mod \
            run.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_vol_fr.$(OBJ_EXT) : calc_vol_fr.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            pgcor.mod \
            pscor.mod \
            compar.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS2) calc_vol_fr.f 
calc_xsi.$(OBJ_EXT) : calc_xsi.f \
            param.mod \
            param1.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
check_ab_m.$(OBJ_EXT) : check_ab_m.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
check_convergence.$(OBJ_EXT) : check_convergence.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            residual.mod \
            toleranc.mod 
check_data_01.$(OBJ_EXT) : check_data_01.f \
            param.mod \
            param1.mod \
            run.mod \
            physprop.mod \
            indices.mod 
check_data_02.$(OBJ_EXT) : check_data_02.f \
            param.mod \
            param1.mod \
            output.mod \
            leqsol.mod \
            geometry.mod 
check_data_03.$(OBJ_EXT) : check_data_03.f \
            param.mod \
            param1.mod \
            geometry.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod 
check_data_04.$(OBJ_EXT) : check_data_04.f \
            param.mod \
            param1.mod \
            run.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            funits.mod 
check_data_05.$(OBJ_EXT) : check_data_05.f \
            param.mod \
            param1.mod \
            physprop.mod \
            funits.mod \
            run.mod \
            indices.mod 
check_data_06.$(OBJ_EXT) : check_data_06.f \
            param.mod \
            param1.mod \
            geometry.mod \
            ic.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod \
            function.inc                                                
check_data_07.$(OBJ_EXT) : check_data_07.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            bc.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
check_data_08.$(OBJ_EXT) : check_data_08.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            is.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
check_data_09.$(OBJ_EXT) : check_data_09.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            rxns.mod \
            indices.mod \
            funits.mod \
            compar.mod 
check_data_20.$(OBJ_EXT) : check_data_20.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            constant.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            visc_g.mod \
            rxns.mod \
            compar.mod \
            function.inc                                                
check_data_30.$(OBJ_EXT) : check_data_30.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            fldvar.mod \
            rxns.mod \
            visc_s.mod \
            geometry.mod \
            run.mod \
            constant.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
check_one_axis.$(OBJ_EXT) : check_one_axis.f \
            param.mod \
            param1.mod \
            funits.mod 
check_plane.$(OBJ_EXT) : check_plane.f \
            funits.mod \
            compar.mod 
compare.$(OBJ_EXT) : compare.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
conv_dif_phi.$(OBJ_EXT) : conv_dif_phi.f \
            param.mod \
            param1.mod \
            run.mod \
            parallel.mod \
            matrix.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            xsi_array.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            output.mod \
            is.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            b_force2.inc                                                
conv_dif_u_g.$(OBJ_EXT) : conv_dif_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            visc_g.mod \
            toleranc.mod \
            physprop.mod \
            fldvar.mod \
            output.mod \
            compar.mod \
            xsi_array.mod \
            tmp_array.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_dif_u_s.$(OBJ_EXT) : conv_dif_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            physprop.mod \
            visc_s.mod \
            toleranc.mod \
            fldvar.mod \
            output.mod \
            compar.mod \
            xsi_array.mod \
            tmp_array.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_dif_v_g.$(OBJ_EXT) : conv_dif_v_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            visc_g.mod \
            toleranc.mod \
            physprop.mod \
            fldvar.mod \
            output.mod \
            compar.mod \
            xsi_array.mod \
            tmp_array.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_dif_v_s.$(OBJ_EXT) : conv_dif_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            physprop.mod \
            visc_s.mod \
            toleranc.mod \
            fldvar.mod \
            output.mod \
            compar.mod \
            xsi_array.mod \
            tmp_array.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_dif_w_g.$(OBJ_EXT) : conv_dif_w_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            visc_g.mod \
            toleranc.mod \
            physprop.mod \
            fldvar.mod \
            output.mod \
            compar.mod \
            xsi_array.mod \
            tmp_array.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_dif_w_s.$(OBJ_EXT) : conv_dif_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            physprop.mod \
            visc_s.mod \
            toleranc.mod \
            fldvar.mod \
            output.mod \
            compar.mod \
            xsi_array.mod \
            tmp_array.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_pp_g.$(OBJ_EXT) : conv_pp_g.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            parallel.mod \
            matrix.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            compar.mod \
            xsi_array.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
conv_rop_g.$(OBJ_EXT) : conv_rop_g.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            parallel.mod \
            matrix.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            compar.mod \
            xsi_array.mod \
            function.inc                                                
conv_rop_s.$(OBJ_EXT) : conv_rop_s.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            parallel.mod \
            matrix.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            pscor.mod \
            compar.mod \
            xsi_array.mod \
            function.inc                                                
conv_source_epp.$(OBJ_EXT) : conv_source_epp.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            parallel.mod \
            matrix.mod \
            constant.mod \
            physprop.mod \
            rxns.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            pscor.mod \
            compar.mod \
            xsi_array.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
copy_a.$(OBJ_EXT) : copy_a.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            physprop.mod \
            function.inc                                                
corner.$(OBJ_EXT) : corner.f \
            param.mod \
            param1.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            matrix.mod \
            corner.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
correct_0.$(OBJ_EXT) : correct_0.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            pgcor.mod \
            ur_facs.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            function.inc                                                
correct_1.$(OBJ_EXT) : correct_1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            geometry.mod \
            pscor.mod \
            ur_facs.mod \
            constant.mod \
            compar.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
dgtsl.$(OBJ_EXT) : dgtsl.f 
dif_u_is.$(OBJ_EXT) : dif_u_is.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            output.mod \
            indices.mod \
            is.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
dif_v_is.$(OBJ_EXT) : dif_v_is.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            output.mod \
            indices.mod \
            is.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
dif_w_is.$(OBJ_EXT) : dif_w_is.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            output.mod \
            indices.mod \
            is.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
discretize.$(OBJ_EXT) : discretize.f \
            param.mod \
            param1.mod 
display_resid.$(OBJ_EXT) : display_resid.f \
            param.mod \
            param1.mod \
            physprop.mod \
            residual.mod \
            fldvar.mod \
            compar.mod 
drag_gs.$(OBJ_EXT) : drag_gs.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
drag_ss.$(OBJ_EXT) : drag_ss.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
eosg.$(OBJ_EXT) : eosg.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            scales.mod \
            sc_p_g1.inc                                                  \
            sc_p_g2.inc                                                 
equal.$(OBJ_EXT) : equal.f \
            param.mod \
            param1.mod \
            indices.mod \
            physprop.mod 
error_routine.$(OBJ_EXT) : error_routine.f \
            funits.mod \
            compar.mod \
            mpi_utility.mod 
exchange.$(OBJ_EXT) : exchange.f \
            param.mod \
            param1.mod 
exit.$(OBJ_EXT) : exit.f \
            compar.mod \
            mpi_utility.mod 
flow_to_vel.$(OBJ_EXT) : flow_to_vel.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            bc.mod \
            indices.mod \
            funits.mod \
            compar.mod 
g_0.$(OBJ_EXT) : g_0.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
get_bc_area.$(OBJ_EXT) : get_bc_area.f \
            param.mod \
            param1.mod \
            geometry.mod \
            bc.mod 
get_data.$(OBJ_EXT) : get_data.f \
            param.mod \
            param1.mod \
            run.mod \
            funits.mod \
            compar.mod \
            gridmap.mod 
get_eq.$(OBJ_EXT) : get_eq.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod 
get_flow_bc.$(OBJ_EXT) : get_flow_bc.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            bc.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
get_hloss.$(OBJ_EXT) : get_hloss.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            bc.mod \
            indices.mod \
            energy.mod 
get_is.$(OBJ_EXT) : get_is.f \
            param.mod \
            param1.mod \
            geometry.mod \
            is.mod \
            indices.mod \
            funits.mod 
get_philoss.$(OBJ_EXT) : get_philoss.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            bc.mod \
            indices.mod \
            energy.mod \
            compar.mod \
            function.inc                                                
get_smass.$(OBJ_EXT) : get_smass.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
get_stats.$(OBJ_EXT) : get_stats.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            funits.mod \
            residual.mod \
            run.mod \
            compar.mod \
            function.inc                                                
get_walls_bc.$(OBJ_EXT) : get_walls_bc.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            bc.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
in_bin_512.$(OBJ_EXT) : in_bin_512.f \
            machine.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
in_bin_512i.$(OBJ_EXT) : in_bin_512i.f \
            machine.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
init_ab_m.$(OBJ_EXT) : init_ab_m.f \
            param.mod \
            param1.mod \
            matrix.mod \
            parallel.mod 
init_fvars.$(OBJ_EXT) : init_fvars.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            indices.mod 
init_namelist.$(OBJ_EXT) : init_namelist.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            geometry.mod \
            ic.mod \
            bc.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            is.mod \
            toleranc.mod \
            scales.mod \
            ur_facs.mod \
            leqsol.mod \
            residual.mod \
            rxns.mod \
            namelist.inc                                                
init_resid.$(OBJ_EXT) : init_resid.f \
            param.mod \
            param1.mod \
            physprop.mod \
            residual.mod 
iterate.$(OBJ_EXT) : iterate.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            output.mod \
            indices.mod \
            funits.mod \
            time_cpu.mod \
            pscor.mod \
            coeff.mod \
            leqsol.mod \
            visc_g.mod \
            pgcor.mod \
            cont.mod \
            compar.mod 
leq_bicgs.$(OBJ_EXT) : leq_bicgs.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            parallel.mod \
            mpi_utility.mod \
            sendrecv.mod \
            function.inc                                                
leq_dsdgmr.$(OBJ_EXT) : leq_dsdgmr.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            dsdgmr_a.mod \
            compar.mod \
            function.inc                                                
leq_dslucs.$(OBJ_EXT) : leq_dslucs.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            dslucs_a.mod \
            compar.mod \
            function.inc                                                
leq_dslugm.$(OBJ_EXT) : leq_dslugm.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            dslugm_a.mod \
            compar.mod \
            function.inc                                                
leq_gmres.$(OBJ_EXT) : leq_gmres.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            debug.mod \
            compar.mod \
            mpi_utility.mod \
            parallel.mod \
            funits.mod \
            gridmap.mod \
            function.inc                                                
leq_igcg.$(OBJ_EXT) : leq_igcg.f \
            param.mod \
            param1.mod \
            geometry.mod \
            matrix.mod \
            igcg_a.mod 
leq_sor.$(OBJ_EXT) : leq_sor.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
line_too_big.$(OBJ_EXT) : line_too_big.f 
location.$(OBJ_EXT) : location.f \
            param.mod \
            param1.mod 
location_check.$(OBJ_EXT) : location_check.f \
            param.mod \
            param1.mod \
            funits.mod \
            geometry.mod 
machine.$(OBJ_EXT) : machine.f \
            machine.mod \
            param.mod \
            run.mod \
            funits.mod 
make_upper_case.$(OBJ_EXT) : make_upper_case.f 
mark_phase_4_cor.$(OBJ_EXT) : mark_phase_4_cor.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            fldvar.mod \
            physprop.mod \
            constant.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS2) mark_phase_4_cor.f 
mfix.$(OBJ_EXT) : mfix.f \
            param.mod \
            param1.mod \
            run.mod \
            time_cpu.mod \
            funits.mod \
            output.mod \
            compar.mod \
            mpi_utility.mod \
            parallel_mpi.mod 
mod_bc_i.$(OBJ_EXT) : mod_bc_i.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
mod_bc_j.$(OBJ_EXT) : mod_bc_j.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
mod_bc_k.$(OBJ_EXT) : mod_bc_k.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
open_file.$(OBJ_EXT) : open_file.f \
            compar.mod 
open_files.$(OBJ_EXT) : open_files.f \
            machine.mod \
            funits.mod \
            compar.mod 
out_array.$(OBJ_EXT) : out_array.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
out_array_c.$(OBJ_EXT) : out_array_c.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
out_array_k.$(OBJ_EXT) : out_array_k.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
out_array_kc.$(OBJ_EXT) : out_array_kc.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
out_bin_512.$(OBJ_EXT) : out_bin_512.f \
            machine.mod 
out_bin_512i.$(OBJ_EXT) : out_bin_512i.f \
            machine.mod 
out_bin_512r.$(OBJ_EXT) : out_bin_512r.f \
            machine.mod 
out_bin_r.$(OBJ_EXT) : out_bin_r.f \
            param.mod 
parse_line.$(OBJ_EXT) : parse_line.f \
            param.mod \
            param1.mod \
            parse.mod \
            compar.mod 
parse_resid_string.$(OBJ_EXT) : parse_resid_string.f \
            param.mod \
            param1.mod \
            physprop.mod \
            residual.mod 
parse_rxn.$(OBJ_EXT) : parse_rxn.f \
            param.mod \
            param1.mod \
            parse.mod \
            rxns.mod 
partial_elim.$(OBJ_EXT) : partial_elim.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            matrix.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            run.mod \
            function.inc                                                
physical_prop.$(OBJ_EXT) : physical_prop.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            toleranc.mod \
            constant.mod \
            compar.mod \
            cp_fun1.inc                                                  \
            function.inc                                                 \
            cp_fun2.inc                                                 
read_namelist.$(OBJ_EXT) : read_namelist.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            geometry.mod \
            ic.mod \
            is.mod \
            bc.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            toleranc.mod \
            funits.mod \
            scales.mod \
            ur_facs.mod \
            leqsol.mod \
            residual.mod \
            rxns.mod \
            compar.mod \
            usrnlst.inc                                                  \
            namelist.inc                                                
read_res0.$(OBJ_EXT) : read_res0.f \
            param.mod \
            param1.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            ic.mod \
            bc.mod \
            is.mod \
            constant.mod \
            funits.mod \
            output.mod \
            scales.mod \
            ur_facs.mod \
            toleranc.mod \
            leqsol.mod \
            tmp_array.mod \
            compar.mod \
            mpi_utility.mod 
read_res1.$(OBJ_EXT) : read_res1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            tmp_array.mod \
            compar.mod \
            mpi_utility.mod 
remove_comment.$(OBJ_EXT) : remove_comment.f 
reset_new.$(OBJ_EXT) : reset_new.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod 
rrates.$(OBJ_EXT) : rrates.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            rxns.mod \
            energy.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
rrates0.$(OBJ_EXT) : rrates0.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            rxns.mod \
            energy.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            funits.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
seek_comment.$(OBJ_EXT) : seek_comment.f 
seek_end.$(OBJ_EXT) : seek_end.f 
set_bc0.$(OBJ_EXT) : set_bc0.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod \
            physprop.mod \
            bc.mod \
            fldvar.mod \
            indices.mod \
            run.mod \
            funits.mod \
            scales.mod \
            boundfunijk.mod \
            sc_p_g1.inc                                                  \
            function.inc                                                 \
            sc_p_g2.inc                                                 
set_bc1.$(OBJ_EXT) : set_bc1.f \
            param.mod \
            param1.mod \
            bc.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
set_constants.$(OBJ_EXT) : set_constants.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            visc_s.mod \
            energy.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            run.mod \
            funits.mod \
            drag.mod 
set_constprop.$(OBJ_EXT) : set_constprop.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            visc_s.mod \
            visc_g.mod \
            energy.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            run.mod \
            funits.mod \
            drag.mod \
            compar.mod \
            function.inc                                                
set_flags.$(OBJ_EXT) : set_flags.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            bc.mod \
            is.mod \
            indices.mod \
            physprop.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
set_fluidbed_p.$(OBJ_EXT) : set_fluidbed_p.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            bc.mod \
            ic.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            funits.mod \
            scales.mod \
            compar.mod \
            sc_p_g1.inc                                                  \
            b_force1.inc                                                 \
            function.inc                                                 \
            b_force2.inc                                                 \
            sc_p_g2.inc                                                 
set_geometry.$(OBJ_EXT) : set_geometry.f \
            param.mod \
            param1.mod \
            run.mod \
            geometry.mod 
set_geometry1.$(OBJ_EXT) : set_geometry1.f \
            param.mod \
            param1.mod \
            parallel.mod \
            run.mod \
            geometry.mod \
            indices.mod 
set_ic.$(OBJ_EXT) : set_ic.f \
            param.mod \
            param1.mod \
            geometry.mod \
            constant.mod \
            physprop.mod \
            ic.mod \
            fldvar.mod \
            visc_g.mod \
            indices.mod \
            scales.mod \
            energy.mod \
            compar.mod \
            sc_p_g1.inc                                                  \
            b_force1.inc                                                 \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            b_force2.inc                                                 \
            sc_p_g2.inc                                                 
set_increments.$(OBJ_EXT) : set_increments.f \
            param.mod \
            param1.mod \
            indices.mod \
            geometry.mod \
            compar.mod \
            physprop.mod \
            fldvar.mod \
            funits.mod \
            function.inc                                                
set_index1.$(OBJ_EXT) : set_index1.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            constant.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
set_index1a.$(OBJ_EXT) : set_index1a.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            boundfunijk.mod \
            function.inc                                                
set_l_scale.$(OBJ_EXT) : set_l_scale.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            visc_g.mod \
            geometry.mod \
            indices.mod 
set_max2.$(OBJ_EXT) : set_max2.f \
            param.mod \
            param1.mod \
            geometry.mod 
set_mw_mix_g.$(OBJ_EXT) : set_mw_mix_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
set_outflow.$(OBJ_EXT) : set_outflow.f \
            param.mod \
            param1.mod \
            bc.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
set_ro_g.$(OBJ_EXT) : set_ro_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
set_wall_bc.$(OBJ_EXT) : set_wall_bc.f \
            param.mod \
            param1.mod \
            bc.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
shift_dxyz.$(OBJ_EXT) : shift_dxyz.f \
            param.mod \
            param1.mod \
            geometry.mod 
solve_continuity.$(OBJ_EXT) : solve_continuity.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            residual.mod \
            cont.mod \
            leqsol.mod \
            ambm.mod 
solve_energy_eq.$(OBJ_EXT) : solve_energy_eq.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            output.mod \
            indices.mod \
            drag.mod \
            residual.mod \
            ur_facs.mod \
            pgcor.mod \
            pscor.mod \
            leqsol.mod \
            bc.mod \
            energy.mod \
            rxns.mod \
            ambm.mod \
            tmp_array.mod \
            tmp_array1.mod \
            compar.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
solve_epp.$(OBJ_EXT) : solve_epp.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            pscor.mod \
            residual.mod \
            leqsol.mod \
            physprop.mod \
            ambm.mod 
solve_granular_energy.$(OBJ_EXT) : solve_granular_energy.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            output.mod \
            indices.mod \
            drag.mod \
            residual.mod \
            ur_facs.mod \
            pgcor.mod \
            pscor.mod \
            leqsol.mod \
            bc.mod \
            energy.mod \
            rxns.mod \
            ambm.mod \
            tmp_array.mod \
            compar.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
solve_lin_eq.$(OBJ_EXT) : solve_lin_eq.f \
            param.mod \
            param1.mod \
            geometry.mod 
solve_pp_g.$(OBJ_EXT) : solve_pp_g.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            pgcor.mod \
            residual.mod \
            leqsol.mod \
            run.mod \
            ambm.mod 
solve_species_eq.$(OBJ_EXT) : solve_species_eq.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            output.mod \
            indices.mod \
            drag.mod \
            residual.mod \
            ur_facs.mod \
            pgcor.mod \
            pscor.mod \
            leqsol.mod \
            bc.mod \
            energy.mod \
            rxns.mod \
            ambm.mod \
            tmp_array.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
solve_vel_star.$(OBJ_EXT) : solve_vel_star.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            output.mod \
            indices.mod \
            drag.mod \
            residual.mod \
            ur_facs.mod \
            pgcor.mod \
            pscor.mod \
            leqsol.mod \
            ambm.mod \
            tmp_array1.mod 
source_granular_energy.$(OBJ_EXT) : source_granular_energy.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            drag.mod \
            geometry.mod \
            fldvar.mod \
            visc_g.mod \
            visc_s.mod \
            trace.mod \
            indices.mod \
            constant.mod \
            compar.mod \
            s_pr1.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            s_pr2.inc                                                   
source_phi.$(OBJ_EXT) : source_phi.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
source_pp_g.$(OBJ_EXT) : source_pp_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            physprop.mod \
            fldvar.mod \
            rxns.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            bc.mod \
            xsi_array.mod \
            compar.mod \
            function.inc                                                
source_rop_g.$(OBJ_EXT) : source_rop_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            rxns.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            compar.mod \
            function.inc                                                
source_rop_s.$(OBJ_EXT) : source_rop_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            rxns.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            pscor.mod \
            compar.mod \
            function.inc                                                
source_u_g.$(OBJ_EXT) : source_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            compar.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_u_s.$(OBJ_EXT) : source_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            bc.mod \
            compar.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_v_g.$(OBJ_EXT) : source_v_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            compar.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_v_s.$(OBJ_EXT) : source_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            bc.mod \
            compar.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_w_g.$(OBJ_EXT) : source_w_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            compar.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_w_s.$(OBJ_EXT) : source_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            bc.mod \
            compar.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
tau_u_g.$(OBJ_EXT) : tau_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
tau_u_s.$(OBJ_EXT) : tau_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
tau_v_g.$(OBJ_EXT) : tau_v_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
tau_v_s.$(OBJ_EXT) : tau_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
tau_w_g.$(OBJ_EXT) : tau_w_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
tau_w_s.$(OBJ_EXT) : tau_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
test_lin_eq.$(OBJ_EXT) : test_lin_eq.f \
            param.mod \
            param1.mod 
time_march.$(OBJ_EXT) : time_march.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            pgcor.mod \
            pscor.mod \
            cont.mod \
            coeff.mod \
            tau_g.mod \
            tau_s.mod \
            visc_g.mod \
            visc_s.mod \
            funits.mod \
            compar.mod 
transfer.$(OBJ_EXT) : transfer.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod 
transport_prop.$(OBJ_EXT) : transport_prop.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            toleranc.mod 
undef_2_0.$(OBJ_EXT) : undef_2_0.f \
            param.mod \
            param1.mod \
            geometry.mod 
under_relax.$(OBJ_EXT) : under_relax.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
update_old.$(OBJ_EXT) : update_old.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            trace.mod \
            visc_s.mod 
usr0.$(OBJ_EXT) : usr0.f 
usr1.$(OBJ_EXT) : usr1.f 
usr2.$(OBJ_EXT) : usr2.f 
usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f \
            usrnlst.inc                                                 
usr_write_out0.$(OBJ_EXT) : usr_write_out0.f 
usr_write_out1.$(OBJ_EXT) : usr_write_out1.f 
vavg_u_g.$(OBJ_EXT) : vavg_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
vavg_u_s.$(OBJ_EXT) : vavg_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vavg_v_g.$(OBJ_EXT) : vavg_v_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
vavg_v_s.$(OBJ_EXT) : vavg_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vavg_w_g.$(OBJ_EXT) : vavg_w_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
vavg_w_s.$(OBJ_EXT) : vavg_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vf_gs_x.$(OBJ_EXT) : vf_gs_x.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
vf_gs_y.$(OBJ_EXT) : vf_gs_y.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
vf_gs_z.$(OBJ_EXT) : vf_gs_z.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
write_ab_m.$(OBJ_EXT) : write_ab_m.f \
            param.mod \
            param1.mod \
            matrix.mod \
            compar.mod \
            mpi_utility.mod 
write_error.$(OBJ_EXT) : write_error.f \
            param.mod \
            param1.mod \
            funits.mod 
write_header.$(OBJ_EXT) : write_header.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            funits.mod \
            compar.mod 
write_out0.$(OBJ_EXT) : write_out0.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            geometry.mod \
            ic.mod \
            bc.mod \
            is.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            funits.mod \
            toleranc.mod \
            scales.mod \
            ur_facs.mod \
            leqsol.mod \
            compar.mod \
            mpi_utility.mod 
write_out1.$(OBJ_EXT) : write_out1.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod 
write_out3.$(OBJ_EXT) : write_out3.f \
            funits.mod \
            compar.mod 
write_res0.$(OBJ_EXT) : write_res0.f \
            param.mod \
            param1.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            ic.mod \
            is.mod \
            bc.mod \
            constant.mod \
            funits.mod \
            output.mod \
            scales.mod \
            ur_facs.mod \
            leqsol.mod \
            toleranc.mod \
            compar.mod \
            mpi_utility.mod 
write_res1.$(OBJ_EXT) : write_res1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            output.mod \
            compar.mod \
            mpi_utility.mod 
write_spx0.$(OBJ_EXT) : write_spx0.f \
            param.mod \
            param1.mod \
            run.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod 
write_spx1.$(OBJ_EXT) : write_spx1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            output.mod \
            compar.mod \
            mpi_utility.mod 
write_table.$(OBJ_EXT) : write_table.f \
            param.mod \
            param1.mod \
            funits.mod 
write_usr0.$(OBJ_EXT) : write_usr0.f 
write_usr1.$(OBJ_EXT) : write_usr1.f 
zero_array.$(OBJ_EXT) : zero_array.f \
            param.mod \
            param1.mod 
zero_norm_vel.$(OBJ_EXT) : zero_norm_vel.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            physprop.mod \
            fldvar.mod \
            indices.mod \
            is.mod \
            compar.mod \
            function.inc                                                
igcg.$(OBJ_EXT) : ./leq_igcg/igcg.f \
            param.mod \
            igcg_i.mod \
            geometry.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./leq_igcg/igcg.f 
dcgs.$(OBJ_EXT) : ./leq_slap/dcgs.f \
            solblk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./leq_slap/dcgs.f 
dgmres.$(OBJ_EXT) : ./leq_slap/dgmres.f \
            solblk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./leq_slap/dgmres.f 
dlaputil.$(OBJ_EXT) : ./leq_slap/dlaputil.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./leq_slap/dlaputil.f 
dmset.$(OBJ_EXT) : ./leq_slap/dmset.f \
            solblk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./leq_slap/dmset.f 
dmvops.$(OBJ_EXT) : ./leq_slap/dmvops.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./leq_slap/dmvops.f 
mach.$(OBJ_EXT) : ./leq_slap/mach.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./leq_slap/mach.f 
xersla.$(OBJ_EXT) : ./leq_slap/xersla.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./leq_slap/xersla.f 
