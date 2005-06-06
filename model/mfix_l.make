.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
mfix.exe : \
    ambm.mod \
    bc.mod \
    boundfunijk3.mod \
    boundfunijk.mod \
    check.mod \
    chischeme.mod \
    coeff.mod \
    constant.mod \
    cont.mod \
    corner.mod \
    drag.mod \
    energy.mod \
    fldvar.mod \
    function.mod \
    funits.mod \
    geometry.mod \
    ic.mod \
    indices.mod \
    is.mod \
    leqsol.mod \
    machine.mod \
    matrix.mod \
    mflux.mod \
    output.mod \
    parallel.mod \
    param1.mod \
    param.mod \
    parse.mod \
    pgcor.mod \
    physprop.mod \
    pscor.mod \
    residual.mod \
    run.mod \
    rxns.mod \
    scalars.mod \
    scales.mod \
    tau_g.mod \
    tau_s.mod \
    time_cpu.mod \
    tmp_array1.mod \
    tmp_array.mod \
    toleranc.mod \
    trace.mod \
    turb.mod \
    ur_facs.mod \
    usr.mod \
    visc_g.mod \
    visc_s.mod \
    vshear.mod \
    xsi_array.mod \
    discretelement.mod \
    compar.mod \
    dbg_util.mod \
    debug.mod \
    gridmap.mod \
    mpi.mod \
    mpi_utility.mod \
    parallel_mpi.mod \
    sendrecv3.mod \
    sendrecv.mod \
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
    bc_phi.$(OBJ_EXT) \
    bc_theta.$(OBJ_EXT) \
    b_m_p_star.$(OBJ_EXT) \
    bound_x.$(OBJ_EXT) \
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
    calc_mflux.$(OBJ_EXT) \
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
    cal_d.$(OBJ_EXT) \
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
    check_mass_balance.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    check_plane.$(OBJ_EXT) \
    cn_extrapol.$(OBJ_EXT) \
    compare.$(OBJ_EXT) \
    conv_dif_phi.$(OBJ_EXT) \
    conv_dif_u_g.$(OBJ_EXT) \
    conv_dif_u_s.$(OBJ_EXT) \
    conv_dif_v_g.$(OBJ_EXT) \
    conv_dif_v_s.$(OBJ_EXT) \
    conv_dif_w_g.$(OBJ_EXT) \
    conv_dif_w_s.$(OBJ_EXT) \
    conv_pp_g.$(OBJ_EXT) \
    conv_rop.$(OBJ_EXT) \
    conv_rop_g.$(OBJ_EXT) \
    conv_rop_s.$(OBJ_EXT) \
    conv_source_epp.$(OBJ_EXT) \
    copy_a.$(OBJ_EXT) \
    corner.$(OBJ_EXT) \
    correct_0.$(OBJ_EXT) \
    correct_1.$(OBJ_EXT) \
    dgtsl.$(OBJ_EXT) \
    dif_u_is.$(OBJ_EXT) \
    dif_v_is.$(OBJ_EXT) \
    dif_w_is.$(OBJ_EXT) \
    discretize.$(OBJ_EXT) \
    display_resid.$(OBJ_EXT) \
    drag_gs.$(OBJ_EXT) \
    drag_ss.$(OBJ_EXT) \
    eosg.$(OBJ_EXT) \
    equal.$(OBJ_EXT) \
    error_routine.$(OBJ_EXT) \
    exchange.$(OBJ_EXT) \
    exit.$(OBJ_EXT) \
    flow_to_vel.$(OBJ_EXT) \
    g_0.$(OBJ_EXT) \
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
    in_bin_512.$(OBJ_EXT) \
    in_bin_512i.$(OBJ_EXT) \
    init_ab_m.$(OBJ_EXT) \
    init_fvars.$(OBJ_EXT) \
    init_namelist.$(OBJ_EXT) \
    init_resid.$(OBJ_EXT) \
    iterate.$(OBJ_EXT) \
    k_epsilon_prop.$(OBJ_EXT) \
    leq_bicgs.$(OBJ_EXT) \
    leq_gmres.$(OBJ_EXT) \
    leq_sor.$(OBJ_EXT) \
    line_too_big.$(OBJ_EXT) \
    location_check.$(OBJ_EXT) \
    location.$(OBJ_EXT) \
    machine.$(OBJ_EXT) \
    make_upper_case.$(OBJ_EXT) \
    mark_phase_4_cor.$(OBJ_EXT) \
    mfix.$(OBJ_EXT) \
    mod_bc_i.$(OBJ_EXT) \
    mod_bc_j.$(OBJ_EXT) \
    mod_bc_k.$(OBJ_EXT) \
    open_file.$(OBJ_EXT) \
    open_files.$(OBJ_EXT) \
    out_array_c.$(OBJ_EXT) \
    out_array.$(OBJ_EXT) \
    out_array_kc.$(OBJ_EXT) \
    out_array_k.$(OBJ_EXT) \
    out_bin_512.$(OBJ_EXT) \
    out_bin_512i.$(OBJ_EXT) \
    out_bin_512r.$(OBJ_EXT) \
    out_bin_r.$(OBJ_EXT) \
    parse_line.$(OBJ_EXT) \
    parse_resid_string.$(OBJ_EXT) \
    parse_rxn.$(OBJ_EXT) \
    partial_elim.$(OBJ_EXT) \
    physical_prop.$(OBJ_EXT) \
    read_namelist.$(OBJ_EXT) \
    read_res0.$(OBJ_EXT) \
    read_res1.$(OBJ_EXT) \
    remove_comment.$(OBJ_EXT) \
    reset_new.$(OBJ_EXT) \
    rrates0.$(OBJ_EXT) \
    rrates.$(OBJ_EXT) \
    rrates_init.$(OBJ_EXT) \
    scalar_prop.$(OBJ_EXT) \
    seek_comment.$(OBJ_EXT) \
    seek_end.$(OBJ_EXT) \
    set_bc0.$(OBJ_EXT) \
    set_bc1.$(OBJ_EXT) \
    set_constants.$(OBJ_EXT) \
    set_constprop.$(OBJ_EXT) \
    set_flags.$(OBJ_EXT) \
    set_fluidbed_p.$(OBJ_EXT) \
    set_geometry1.$(OBJ_EXT) \
    set_geometry.$(OBJ_EXT) \
    set_ic.$(OBJ_EXT) \
    set_increments3.$(OBJ_EXT) \
    set_increments.$(OBJ_EXT) \
    set_index1a3.$(OBJ_EXT) \
    set_index1a.$(OBJ_EXT) \
    set_index1.$(OBJ_EXT) \
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
    solve_k_epsilon_eq.$(OBJ_EXT) \
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
    tau_u_g.$(OBJ_EXT) \
    tau_u_s.$(OBJ_EXT) \
    tau_v_g.$(OBJ_EXT) \
    tau_v_s.$(OBJ_EXT) \
    tau_w_g.$(OBJ_EXT) \
    tau_w_s.$(OBJ_EXT) \
    test_lin_eq.$(OBJ_EXT) \
    time_march.$(OBJ_EXT) \
    transfer.$(OBJ_EXT) \
    transport_prop.$(OBJ_EXT) \
    undef_2_0.$(OBJ_EXT) \
    under_relax.$(OBJ_EXT) \
    update_old.$(OBJ_EXT) \
    usr0.$(OBJ_EXT) \
    usr1.$(OBJ_EXT) \
    usr2.$(OBJ_EXT) \
    usr3.$(OBJ_EXT) \
    usr_init_namelist.$(OBJ_EXT) \
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
    write_ab_m.$(OBJ_EXT) \
    write_ab_m_var.$(OBJ_EXT) \
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
    zero_array.$(OBJ_EXT) \
    zero_norm_vel.$(OBJ_EXT) \
    add_part_to_link_list.$(OBJ_EXT) \
    calc_app_coh_force.$(OBJ_EXT) \
    calc_cap_coh_force.$(OBJ_EXT) \
    calc_cohesive_forces.$(OBJ_EXT) \
    calc_esc_coh_force.$(OBJ_EXT) \
    calc_square_well.$(OBJ_EXT) \
    calc_van_der_waals.$(OBJ_EXT) \
    check_link.$(OBJ_EXT) \
    check_sw_wall_interaction.$(OBJ_EXT) \
    check_vdw_wall_interaction.$(OBJ_EXT) \
    initialize_cohesion_parameters.$(OBJ_EXT) \
    initialize_coh_int_search.$(OBJ_EXT) \
    linked_interaction_eval.$(OBJ_EXT) \
    remove_part_from_link_list.$(OBJ_EXT) \
    unlinked_interaction_eval.$(OBJ_EXT) \
    update_search_grids.$(OBJ_EXT) \
    calc_force_des.$(OBJ_EXT) \
    cfassign.$(OBJ_EXT) \
    cffctow.$(OBJ_EXT) \
    cffn.$(OBJ_EXT) \
    cffnwall.$(OBJ_EXT) \
    cfft.$(OBJ_EXT) \
    cfftwall.$(OBJ_EXT) \
    cfincrementaloverlaps.$(OBJ_EXT) \
    cfnewvalues.$(OBJ_EXT) \
    cfnocontact.$(OBJ_EXT) \
    cfnormal.$(OBJ_EXT) \
    cfoutofbox.$(OBJ_EXT) \
    cfperiodicwallneighbourx.$(OBJ_EXT) \
    cfperiodicwallneighboury.$(OBJ_EXT) \
    cfperiodicwallneighbourz.$(OBJ_EXT) \
    cfperiodicwallx.$(OBJ_EXT) \
    cfperiodicwally.$(OBJ_EXT) \
    cfperiodicwallz.$(OBJ_EXT) \
    cfrelvel.$(OBJ_EXT) \
    cfslide.$(OBJ_EXT) \
    cfslidewall.$(OBJ_EXT) \
    cftangent.$(OBJ_EXT) \
    cftotaloverlaps.$(OBJ_EXT) \
    cfupdateold.$(OBJ_EXT) \
    cfvrn.$(OBJ_EXT) \
    cfvrt.$(OBJ_EXT) \
    cfwallcontact.$(OBJ_EXT) \
    cfwallposvel.$(OBJ_EXT) \
    des_calc_d.$(OBJ_EXT) \
    des_granular_temperature.$(OBJ_EXT) \
    des_init_namelist.$(OBJ_EXT) \
    des_inlet_outlet.$(OBJ_EXT) \
    des_time_march.$(OBJ_EXT) \
    drag_fgs.$(OBJ_EXT) \
    gas_drag.$(OBJ_EXT) \
    make_arrays_des.$(OBJ_EXT) \
    neighbour.$(OBJ_EXT) \
    nsquare.$(OBJ_EXT) \
    octree.$(OBJ_EXT) \
    particles_in_cell.$(OBJ_EXT) \
    periodic_wall_calc_force_des.$(OBJ_EXT) \
    pressure_drop.$(OBJ_EXT) \
    print_vel.$(OBJ_EXT) \
    quadtree.$(OBJ_EXT) \
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
    bc_mod.$(OBJ_EXT) \
    bc_phi.$(OBJ_EXT) \
    bc_theta.$(OBJ_EXT) \
    b_m_p_star.$(OBJ_EXT) \
    boundfunijk3_mod.$(OBJ_EXT) \
    boundfunijk_mod.$(OBJ_EXT) \
    bound_x.$(OBJ_EXT) \
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
    calc_mflux.$(OBJ_EXT) \
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
    cal_d.$(OBJ_EXT) \
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
    check_mass_balance.$(OBJ_EXT) \
    check_mod.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    check_plane.$(OBJ_EXT) \
    chischeme_mod.$(OBJ_EXT) \
    cn_extrapol.$(OBJ_EXT) \
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
    conv_rop.$(OBJ_EXT) \
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
    k_epsilon_prop.$(OBJ_EXT) \
    leq_bicgs.$(OBJ_EXT) \
    leq_gmres.$(OBJ_EXT) \
    leqsol_mod.$(OBJ_EXT) \
    leq_sor.$(OBJ_EXT) \
    line_too_big.$(OBJ_EXT) \
    location_check.$(OBJ_EXT) \
    location.$(OBJ_EXT) \
    machine.$(OBJ_EXT) \
    machine_mod.$(OBJ_EXT) \
    make_upper_case.$(OBJ_EXT) \
    mark_phase_4_cor.$(OBJ_EXT) \
    matrix_mod.$(OBJ_EXT) \
    mfix.$(OBJ_EXT) \
    mflux_mod.$(OBJ_EXT) \
    mod_bc_i.$(OBJ_EXT) \
    mod_bc_j.$(OBJ_EXT) \
    mod_bc_k.$(OBJ_EXT) \
    open_file.$(OBJ_EXT) \
    open_files.$(OBJ_EXT) \
    out_array_c.$(OBJ_EXT) \
    out_array.$(OBJ_EXT) \
    out_array_kc.$(OBJ_EXT) \
    out_array_k.$(OBJ_EXT) \
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
    rrates0.$(OBJ_EXT) \
    rrates.$(OBJ_EXT) \
    rrates_init.$(OBJ_EXT) \
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
    set_geometry1.$(OBJ_EXT) \
    set_geometry.$(OBJ_EXT) \
    set_ic.$(OBJ_EXT) \
    set_increments3.$(OBJ_EXT) \
    set_increments.$(OBJ_EXT) \
    set_index1a3.$(OBJ_EXT) \
    set_index1a.$(OBJ_EXT) \
    set_index1.$(OBJ_EXT) \
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
    solve_k_epsilon_eq.$(OBJ_EXT) \
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
    turb_mod.$(OBJ_EXT) \
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
    write_ab_m_var.$(OBJ_EXT) \
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
    add_part_to_link_list.$(OBJ_EXT) \
    calc_app_coh_force.$(OBJ_EXT) \
    calc_cap_coh_force.$(OBJ_EXT) \
    calc_cohesive_forces.$(OBJ_EXT) \
    calc_esc_coh_force.$(OBJ_EXT) \
    calc_square_well.$(OBJ_EXT) \
    calc_van_der_waals.$(OBJ_EXT) \
    check_link.$(OBJ_EXT) \
    check_sw_wall_interaction.$(OBJ_EXT) \
    check_vdw_wall_interaction.$(OBJ_EXT) \
    initialize_cohesion_parameters.$(OBJ_EXT) \
    initialize_coh_int_search.$(OBJ_EXT) \
    linked_interaction_eval.$(OBJ_EXT) \
    remove_part_from_link_list.$(OBJ_EXT) \
    unlinked_interaction_eval.$(OBJ_EXT) \
    update_search_grids.$(OBJ_EXT) \
    calc_force_des.$(OBJ_EXT) \
    cfassign.$(OBJ_EXT) \
    cffctow.$(OBJ_EXT) \
    cffn.$(OBJ_EXT) \
    cffnwall.$(OBJ_EXT) \
    cfft.$(OBJ_EXT) \
    cfftwall.$(OBJ_EXT) \
    cfincrementaloverlaps.$(OBJ_EXT) \
    cfnewvalues.$(OBJ_EXT) \
    cfnocontact.$(OBJ_EXT) \
    cfnormal.$(OBJ_EXT) \
    cfoutofbox.$(OBJ_EXT) \
    cfperiodicwallneighbourx.$(OBJ_EXT) \
    cfperiodicwallneighboury.$(OBJ_EXT) \
    cfperiodicwallneighbourz.$(OBJ_EXT) \
    cfperiodicwallx.$(OBJ_EXT) \
    cfperiodicwally.$(OBJ_EXT) \
    cfperiodicwallz.$(OBJ_EXT) \
    cfrelvel.$(OBJ_EXT) \
    cfslide.$(OBJ_EXT) \
    cfslidewall.$(OBJ_EXT) \
    cftangent.$(OBJ_EXT) \
    cftotaloverlaps.$(OBJ_EXT) \
    cfupdateold.$(OBJ_EXT) \
    cfvrn.$(OBJ_EXT) \
    cfvrt.$(OBJ_EXT) \
    cfwallcontact.$(OBJ_EXT) \
    cfwallposvel.$(OBJ_EXT) \
    des_calc_d.$(OBJ_EXT) \
    des_granular_temperature.$(OBJ_EXT) \
    des_init_namelist.$(OBJ_EXT) \
    des_inlet_outlet.$(OBJ_EXT) \
    des_time_march.$(OBJ_EXT) \
    discretelement_mod.$(OBJ_EXT) \
    drag_fgs.$(OBJ_EXT) \
    gas_drag.$(OBJ_EXT) \
    make_arrays_des.$(OBJ_EXT) \
    neighbour.$(OBJ_EXT) \
    nsquare.$(OBJ_EXT) \
    octree.$(OBJ_EXT) \
    particles_in_cell.$(OBJ_EXT) \
    periodic_wall_calc_force_des.$(OBJ_EXT) \
    pressure_drop.$(OBJ_EXT) \
    print_vel.$(OBJ_EXT) \
    quadtree.$(OBJ_EXT) \
    compar_mod.$(OBJ_EXT) \
    dbg_util_mod.$(OBJ_EXT) \
    debug_mod.$(OBJ_EXT) \
    gridmap_mod.$(OBJ_EXT) \
    mpi_mod.$(OBJ_EXT) \
    mpi_utility_mod.$(OBJ_EXT) \
    parallel_mpi_mod.$(OBJ_EXT) \
    sendrecv3_mod.$(OBJ_EXT) \
    sendrecv_mod.$(OBJ_EXT) \
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
boundfunijk3.mod : boundfunijk3_mod.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk3_mod.f 
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
check.mod : check_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_mod.f 
chischeme.mod : chischeme_mod.f \
            param.mod \
            param1.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) chischeme_mod.f 
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
mflux.mod : mflux_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) mflux_mod.f 
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
scalars.mod : scalars_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scalars_mod.f 
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
turb.mod : turb_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) turb_mod.f 
ur_facs.mod : ur_facs_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ur_facs_mod.f 
usr.mod : usr_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_mod.f 
visc_g.mod : visc_g_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_g_mod.f 
visc_s.mod : visc_s_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_s_mod.f 
vshear.mod : vshear_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) vshear_mod.f 
xsi_array.mod : xsi_array_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) xsi_array_mod.f 
discretelement.mod : ./des/discretelement_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/discretelement_mod.f 
compar.mod : ./dmp_modules/compar_mod.f \
            mpi.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/compar_mod.f 
dbg_util.mod : ./dmp_modules/dbg_util_mod.f \
            compar.mod \
            geometry.mod \
            parallel_mpi.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/dbg_util_mod.f 
debug.mod : ./dmp_modules/debug_mod.f \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/debug_mod.f 
gridmap.mod : ./dmp_modules/gridmap_mod.f \
            mpi_utility.mod \
            parallel_mpi.mod \
            geometry.mod \
            sendrecv.mod \
            compar.mod \
            run.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/gridmap_mod.f 
mpi.mod : ./dmp_modules/mpi_mod.f \
            mpif.h                                                      
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
sendrecv3.mod : ./dmp_modules/sendrecv3_mod.f \
            parallel_mpi.mod \
            debug.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            mpi.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv3_mod.f 
sendrecv.mod : ./dmp_modules/sendrecv_mod.f \
            parallel_mpi.mod \
            debug.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            mpi.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv_mod.f 
adjust_a_u_g.$(OBJ_EXT) : adjust_a_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            usr.mod \
            compar.mod \
            sendrecv.mod \
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
            sendrecv.mod \
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
            sendrecv.mod \
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
            sendrecv.mod \
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
            sendrecv.mod \
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
            sendrecv.mod \
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
            sendrecv.mod \
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
ar             .$(OBJ_EXT) : tant           .f g
            param.mod g
            param1.mod g
            fldvar.mod g
            mflux.mod g
            physprop.mod g
            run.mod g
            parallel.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            sendrecv.mod g
            xsi_array.mod g
            mpi_utility.mod g
            matrix.mod g
            toleranc.mod g
            sendrecv3.mod g
            tmp_array.mod g
            vshear.mod g
            scales.mod g
            constant.mod g
            visc_s.mod g
            output.mod g
            is.mod g
            visc_g.mod g
            pgcor.mod g
            ur_facs.mod g
            funits.mod g
            time_cpu.mod g
            pscor.mod g
            coeff.mod g
            leqsol.mod g
            cont.mod g
            scalars.mod g
            discretelement.mod g
            bc.mod g
            param1.mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            nt.mod g
            .mod g
            m.mod g
            m1.mod g
            ar.mod g
            x.mod g
            prop.mod g
            .mod g
            llel.mod g
            etry.mod g
            ces.mod g
            ar.mod g
            recv.mod g
            array.mod g
            utility.mod g
            ix.mod g
            ranc.mod g
            recv3.mod g
            array.mod g
            ar.mod g
            es.mod g
            tant.mod g
            _s.mod g
            ut.mod g
            .mod g
            _g.mod g
            r.mod g
            acs.mod g
            ts.mod g
            _cpu.mod g
            r.mod g
            f.mod g
            ol.mod g
            .mod g
            ars.mod g
            retelement.mod g
            m1.mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            .mod g
            function.inc                                                 g
            fun_avg1.inc                                                 g
            fun_avg2.inc                                                 g
            function3.inc                                                g
            ep_s1.inc                                                    g
            ep_s2.inc                                                    g
            radtn1.inc                                                   g
            radtn2.inc                                                   g
            b_force1.inc                                                 g
            b_force2.inc                                                
bc_phi.$(OBJ_EXT) : bc_phi.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            geometry.mod g
            output.mod g
            indices.mod g
            bc.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
bc_theta.$(OBJ_EXT) : bc_theta.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            geometry.mod g
            output.mod g
            indices.mod g
            bc.mod g
            compar.mod g
            mpi_utility.mod g
            turb.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
b_m_p_star.$(OBJ_EXT) : b_m_p_star.f g
            param.mod g
            param1.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            run.mod g
            rxns.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
bound_x.$(OBJ_EXT) : bound_x.f g
            param.mod g
            param1.mod 
calc_cell.$(OBJ_EXT) : calc_cell.f g
            param.mod g
            param1.mod 
calc_coeff.$(OBJ_EXT) : calc_coeff.f g
            param.mod g
            param1.mod g
            physprop.mod g
            rxns.mod g
            funits.mod g
            compar.mod 
calc_d.$(OBJ_EXT) : calc_d.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            run.mod g
            scales.mod g
            compar.mod g
            sendrecv.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
calc_dif_g.$(OBJ_EXT) : calc_dif_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            constant.mod g
            compar.mod g
            sendrecv.mod g
            run.mod g
            function.inc                                                
calc_dif_s.$(OBJ_EXT) : calc_dif_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            constant.mod g
            toleranc.mod g
            compar.mod g
            sendrecv.mod g
            run.mod g
            function.inc                                                
calc_drag.$(OBJ_EXT) : calc_drag.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            run.mod g
            drag.mod g
            compar.mod g
            discretelement.mod 
calc_e.$(OBJ_EXT) : calc_e.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            run.mod g
            constant.mod g
            compar.mod g
            sendrecv.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
calc_gama.$(OBJ_EXT) : calc_gama.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            energy.mod g
            rxns.mod g
            indices.mod g
            compar.mod g
            sendrecv.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
calc_grbdry.$(OBJ_EXT) : calc_grbdry.f g
            param.mod g
            param1.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            run.mod g
            turb.mod g
            visc_s.mod g
            geometry.mod g
            indices.mod g
            bc.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
calc_k_cp.$(OBJ_EXT) : calc_k_cp.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            physprop.mod g
            indices.mod g
            pscor.mod g
            geometry.mod g
            constant.mod g
            run.mod g
            visc_s.mod g
            trace.mod g
            compar.mod g
            sendrecv.mod g
            ep_s1.inc                                                    g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            ep_s2.inc                                                   
calc_k_g.$(OBJ_EXT) : calc_k_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            constant.mod g
            compar.mod g
            run.mod g
            sendrecv.mod g
            function.inc                                                
calc_k_s.$(OBJ_EXT) : calc_k_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            constant.mod g
            toleranc.mod g
            compar.mod g
            sendrecv.mod g
            run.mod g
            function.inc                                                
calc_mflux.$(OBJ_EXT) : calc_mflux.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            mflux.mod g
            physprop.mod g
            run.mod g
            parallel.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
calc_mu_g.$(OBJ_EXT) : calc_mu_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            visc_g.mod g
            visc_s.mod g
            indices.mod g
            constant.mod g
            toleranc.mod g
            compar.mod g
            drag.mod g
            run.mod g
            turb.mod g
            sendrecv.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg2.inc                                                
calc_mu_s.$(OBJ_EXT) : calc_mu_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            drag.mod g
            run.mod g
            geometry.mod g
            fldvar.mod g
            visc_g.mod g
            visc_s.mod g
            trace.mod g
            turb.mod g
            indices.mod g
            constant.mod g
            toleranc.mod g
            vshear.mod g
            compar.mod g
            sendrecv.mod g
            s_pr1.inc                                                    g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg2.inc                                                 g
            s_pr2.inc                                                   
calc_mw.$(OBJ_EXT) : calc_mw.f g
            param.mod g
            param1.mod g
            toleranc.mod 
calc_outflow.$(OBJ_EXT) : calc_outflow.f g
            param.mod g
            param1.mod g
            bc.mod g
            fldvar.mod g
            indices.mod g
            physprop.mod g
            geometry.mod g
            compar.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
calc_p_star.$(OBJ_EXT) : calc_p_star.f g
            param.mod g
            param1.mod g
            parallel.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            constant.mod g
            pgcor.mod g
            pscor.mod g
            ur_facs.mod g
            residual.mod g
            compar.mod g
            fldvar.mod g
            toleranc.mod g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
calc_resid.$(OBJ_EXT) : calc_resid.f g
            param.mod g
            param1.mod g
            matrix.mod g
            parallel.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            mpi_utility.mod g
            fldvar.mod g
            run.mod g
            bc.mod g
            constant.mod g
            physprop.mod g
            residual.mod g
            rxns.mod g
            function.inc                                                
calc_s_ddot_s.$(OBJ_EXT) : calc_s_ddot_s.f g
            param.mod g
            param1.mod g
            constant.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                
calc_trd_g.$(OBJ_EXT) : calc_trd_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            geometry.mod g
            fldvar.mod g
            indices.mod g
            compar.mod g
            sendrecv.mod g
            function.inc                                                
calc_trd_s.$(OBJ_EXT) : calc_trd_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            geometry.mod g
            fldvar.mod g
            indices.mod g
            physprop.mod g
            compar.mod g
            sendrecv.mod g
            function.inc                                                
calc_u_friction.$(OBJ_EXT) : calc_u_friction.f g
            param.mod g
            param1.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            run.mod g
            turb.mod g
            visc_s.mod g
            geometry.mod g
            indices.mod g
            bc.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
calc_vol_fr.$(OBJ_EXT) : calc_vol_fr.f g
            param.mod g
            param1.mod g
            parallel.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            constant.mod g
            pgcor.mod g
            pscor.mod g
            compar.mod g
            sendrecv.mod g
            ep_s1.inc                                                    g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS2) calc_vol_fr.f 
calc_xsi.$(OBJ_EXT) : calc_xsi.f g
            param.mod g
            param1.mod g
            run.mod g
            geometry.mod g
            indices.mod g
            vshear.mod g
            chischeme.mod g
            compar.mod g
            sendrecv.mod g
            xsi1.inc                                                     g
            function.inc                                                 g
            xsi2.inc                                                    
cal_d.$(OBJ_EXT) : cal_d.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            rxns.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            tau_s.mod g
            bc.mod g
            vshear.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
check_ab_m.$(OBJ_EXT) : check_ab_m.f g
            param.mod g
            param1.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
check_convergence.$(OBJ_EXT) : check_convergence.f g
            param.mod g
            param1.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            run.mod g
            residual.mod g
            toleranc.mod g
            mpi_utility.mod 
check_data_01.$(OBJ_EXT) : check_data_01.f g
            param.mod g
            param1.mod g
            constant.mod g
            run.mod g
            physprop.mod g
            indices.mod g
            scalars.mod g
            funits.mod 
check_data_02.$(OBJ_EXT) : check_data_02.f g
            param.mod g
            param1.mod g
            output.mod g
            leqsol.mod g
            geometry.mod g
            run.mod g
            rxns.mod 
check_data_03.$(OBJ_EXT) : check_data_03.f g
            param.mod g
            param1.mod g
            geometry.mod g
            bc.mod g
            funits.mod g
            compar.mod g
            mpi_utility.mod 
check_data_04.$(OBJ_EXT) : check_data_04.f g
            param.mod g
            param1.mod g
            run.mod g
            indices.mod g
            physprop.mod g
            constant.mod g
            funits.mod 
check_data_05.$(OBJ_EXT) : check_data_05.f g
            param.mod g
            param1.mod g
            physprop.mod g
            funits.mod g
            run.mod g
            indices.mod 
check_data_06.$(OBJ_EXT) : check_data_06.f g
            param.mod g
            param1.mod g
            geometry.mod g
            ic.mod g
            fldvar.mod g
            physprop.mod g
            run.mod g
            indices.mod g
            funits.mod g
            scalars.mod g
            compar.mod g
            mpi_utility.mod g
            sendrecv.mod g
            function.inc                                                
check_data_07.$(OBJ_EXT) : check_data_07.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            run.mod g
            bc.mod g
            indices.mod g
            funits.mod g
            scalars.mod g
            compar.mod g
            sendrecv.mod g
            function.inc                                                
check_data_08.$(OBJ_EXT) : check_data_08.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            run.mod g
            is.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            function.inc                                                
check_data_09.$(OBJ_EXT) : check_data_09.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            run.mod g
            rxns.mod g
            indices.mod g
            funits.mod g
            compar.mod 
check_data_20.$(OBJ_EXT) : check_data_20.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            fldvar.mod g
            run.mod g
            geometry.mod g
            constant.mod g
            physprop.mod g
            indices.mod g
            funits.mod g
            visc_g.mod g
            rxns.mod g
            scalars.mod g
            compar.mod g
            sendrecv.mod g
            function.inc                                                
check_data_30.$(OBJ_EXT) : check_data_30.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            fldvar.mod g
            rxns.mod g
            visc_s.mod g
            visc_g.mod g
            geometry.mod g
            run.mod g
            constant.mod g
            physprop.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            mpi_utility.mod g
            function.inc                                                
check_mass_balance.$(OBJ_EXT) : check_mass_balance.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            fldvar.mod g
            rxns.mod g
            geometry.mod g
            run.mod g
            bc.mod g
            constant.mod g
            physprop.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            mpi_utility.mod g
            output.mod g
            check.mod g
            parallel.mod g
            matrix.mod g
            function.inc                                                
check_one_axis.$(OBJ_EXT) : check_one_axis.f g
            param.mod g
            param1.mod g
            funits.mod 
check_plane.$(OBJ_EXT) : check_plane.f g
            funits.mod g
            compar.mod 
cn_extrapol.$(OBJ_EXT) : cn_extrapol.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            scalars.mod g
            trace.mod g
            run.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            mpi_utility.mod g
            function.inc                                                
compare.$(OBJ_EXT) : compare.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
conv_dif_phi.$(OBJ_EXT) : conv_dif_phi.f g
            param.mod g
            param1.mod g
            run.mod g
            geometry.mod g
            compar.mod g
            sendrecv.mod g
            xsi_array.mod g
            mpi_utility.mod g
            indices.mod g
            parallel.mod g
            matrix.mod g
            toleranc.mod g
            sendrecv3.mod g
            tmp_array.mod g
            vshear.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            output.mod g
            is.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            function3.inc                                                g
            ep_s1.inc                                                    g
            ep_s2.inc                                                   
conv_dif_u_g.$(OBJ_EXT) : conv_dif_u_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            run.mod g
            visc_g.mod g
            compar.mod g
            toleranc.mod g
            physprop.mod g
            fldvar.mod g
            output.mod g
            mflux.mod g
            vshear.mod g
            xsi_array.mod g
            tmp_array.mod g
            sendrecv.mod g
            sendrecv3.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_dif_u_s.$(OBJ_EXT) : conv_dif_u_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            run.mod g
            physprop.mod g
            visc_s.mod g
            compar.mod g
            toleranc.mod g
            fldvar.mod g
            output.mod g
            mflux.mod g
            xsi_array.mod g
            tmp_array.mod g
            sendrecv.mod g
            sendrecv3.mod g
            vshear.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_dif_v_g.$(OBJ_EXT) : conv_dif_v_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            run.mod g
            visc_g.mod g
            compar.mod g
            toleranc.mod g
            physprop.mod g
            fldvar.mod g
            output.mod g
            mflux.mod g
            xsi_array.mod g
            vshear.mod g
            tmp_array.mod g
            sendrecv.mod g
            sendrecv3.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_dif_v_s.$(OBJ_EXT) : conv_dif_v_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            run.mod g
            physprop.mod g
            visc_s.mod g
            compar.mod g
            toleranc.mod g
            fldvar.mod g
            output.mod g
            mflux.mod g
            xsi_array.mod g
            tmp_array.mod g
            sendrecv.mod g
            sendrecv3.mod g
            vshear.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_dif_w_g.$(OBJ_EXT) : conv_dif_w_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            run.mod g
            visc_g.mod g
            compar.mod g
            toleranc.mod g
            physprop.mod g
            fldvar.mod g
            output.mod g
            mflux.mod g
            xsi_array.mod g
            tmp_array.mod g
            sendrecv.mod g
            sendrecv3.mod g
            vshear.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_dif_w_s.$(OBJ_EXT) : conv_dif_w_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            run.mod g
            physprop.mod g
            visc_s.mod g
            compar.mod g
            toleranc.mod g
            fldvar.mod g
            output.mod g
            mflux.mod g
            xsi_array.mod g
            tmp_array.mod g
            sendrecv.mod g
            sendrecv3.mod g
            vshear.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_pp_g.$(OBJ_EXT) : conv_pp_g.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            run.mod g
            parallel.mod g
            matrix.mod g
            physprop.mod g
            geometry.mod g
            indices.mod g
            pgcor.mod g
            compar.mod g
            mflux.mod g
            function.inc                                                
conv_rop.$(OBJ_EXT) : conv_rop.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            mflux.mod g
            physprop.mod g
            run.mod g
            parallel.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            xsi_array.mod g
            function.inc                                                
conv_rop_g.$(OBJ_EXT) : conv_rop_g.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            run.mod g
            compar.mod g
            parallel.mod g
            matrix.mod g
            physprop.mod g
            geometry.mod g
            indices.mod g
            pgcor.mod g
            xsi_array.mod g
            function.inc                                                
conv_rop_s.$(OBJ_EXT) : conv_rop_s.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            run.mod g
            compar.mod g
            parallel.mod g
            matrix.mod g
            physprop.mod g
            geometry.mod g
            indices.mod g
            pgcor.mod g
            pscor.mod g
            xsi_array.mod g
            function.inc                                                
conv_source_epp.$(OBJ_EXT) : conv_source_epp.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            run.mod g
            geometry.mod g
            compar.mod g
            sendrecv.mod g
            xsi_array.mod g
            parallel.mod g
            matrix.mod g
            constant.mod g
            physprop.mod g
            rxns.mod g
            indices.mod g
            pgcor.mod g
            pscor.mod g
            vshear.mod g
            ep_s1.inc                                                    g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            ep_s2.inc                                                   
copy_a.$(OBJ_EXT) : copy_a.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            physprop.mod g
            function.inc                                                
corner.$(OBJ_EXT) : corner.f g
            param.mod g
            param1.mod g
            geometry.mod g
            physprop.mod g
            indices.mod g
            matrix.mod g
            corner.mod g
            funits.mod g
            compar.mod g
            function.inc                                                
correct_0.$(OBJ_EXT) : correct_0.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            pgcor.mod g
            ur_facs.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            compar.mod g
            function.inc                                                
correct_1.$(OBJ_EXT) : correct_1.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            physprop.mod g
            indices.mod g
            geometry.mod g
            pscor.mod g
            ur_facs.mod g
            constant.mod g
            compar.mod g
            sendrecv.mod g
            ep_s1.inc                                                    g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            ep_s2.inc                                                   
dgtsl.$(OBJ_EXT) : dgtsl.f 
dif_u_is.$(OBJ_EXT) : dif_u_is.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            output.mod g
            indices.mod g
            is.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
dif_v_is.$(OBJ_EXT) : dif_v_is.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            output.mod g
            indices.mod g
            is.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
dif_w_is.$(OBJ_EXT) : dif_w_is.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            output.mod g
            indices.mod g
            is.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
discretize.$(OBJ_EXT) : discretize.f g
            param.mod g
            param1.mod g
            run.mod 
display_resid.$(OBJ_EXT) : display_resid.f g
            param.mod g
            param1.mod g
            physprop.mod g
            residual.mod g
            fldvar.mod g
            compar.mod g
            geometry.mod 
drag_gs.$(OBJ_EXT) : drag_gs.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            run.mod g
            constant.mod g
            compar.mod g
            drag.mod g
            sendrecv.mod g
            discretelement.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
drag_ss.$(OBJ_EXT) : drag_ss.f g
            param.mod g
            param1.mod g
            parallel.mod g
            constant.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            compar.mod g
            sendrecv.mod g
            drag.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                
eosg.$(OBJ_EXT) : eosg.f g
            param.mod g
            param1.mod g
            constant.mod g
            physprop.mod g
            scales.mod g
            sc_p_g1.inc                                                  g
            sc_p_g2.inc                                                 
equal.$(OBJ_EXT) : equal.f g
            param.mod g
            param1.mod g
            indices.mod g
            physprop.mod 
error_routine.$(OBJ_EXT) : error_routine.f g
            funits.mod g
            compar.mod g
            mpi_utility.mod 
exchange.$(OBJ_EXT) : exchange.f g
            param.mod g
            param1.mod g
            compar.mod 
exit.$(OBJ_EXT) : exit.f g
            funits.mod g
            compar.mod g
            mpi_utility.mod 
flow_to_vel.$(OBJ_EXT) : flow_to_vel.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            run.mod g
            bc.mod g
            indices.mod g
            funits.mod g
            compar.mod 
g_0.$(OBJ_EXT) : g_0.f g
            param.mod g
            param1.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg1.inc                                                 g
            fun_avg2.inc                                                
get_bc_area.$(OBJ_EXT) : get_bc_area.f g
            param.mod g
            param1.mod g
            geometry.mod g
            bc.mod g
            compar.mod 
get_data.$(OBJ_EXT) : get_data.f g
            param.mod g
            param1.mod g
            run.mod g
            funits.mod g
            compar.mod g
            gridmap.mod 
get_eq.$(OBJ_EXT) : get_eq.f g
            param.mod g
            param1.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            indices.mod 
get_flow_bc.$(OBJ_EXT) : get_flow_bc.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            bc.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            sendrecv.mod g
            function.inc                                                
get_hloss.$(OBJ_EXT) : get_hloss.f g
            param.mod g
            param1.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            bc.mod g
            indices.mod g
            energy.mod 
get_is.$(OBJ_EXT) : get_is.f g
            param.mod g
            param1.mod g
            geometry.mod g
            is.mod g
            indices.mod g
            funits.mod g
            compar.mod 
get_philoss.$(OBJ_EXT) : get_philoss.f g
            param.mod g
            param1.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            bc.mod g
            indices.mod g
            energy.mod g
            compar.mod g
            function.inc                                                
get_smass.$(OBJ_EXT) : get_smass.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            indices.mod g
            compar.mod g
            mpi_utility.mod g
            function.inc                                                
get_stats.$(OBJ_EXT) : get_stats.f g
            param.mod g
            param1.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            indices.mod g
            funits.mod g
            residual.mod g
            run.mod g
            compar.mod g
            function.inc                                                
get_walls_bc.$(OBJ_EXT) : get_walls_bc.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            bc.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            sendrecv.mod g
            function.inc                                                
in_bin_512.$(OBJ_EXT) : in_bin_512.f g
            machine.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
in_bin_512i.$(OBJ_EXT) : in_bin_512i.f g
            machine.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
init_ab_m.$(OBJ_EXT) : init_ab_m.f g
            param.mod g
            param1.mod g
            matrix.mod g
            parallel.mod g
            compar.mod 
init_fvars.$(OBJ_EXT) : init_fvars.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            geometry.mod g
            physprop.mod g
            indices.mod g
            scalars.mod g
            rxns.mod g
            run.mod g
            compar.mod 
init_namelist.$(OBJ_EXT) : init_namelist.f g
            param.mod g
            param1.mod g
            run.mod g
            output.mod g
            physprop.mod g
            geometry.mod g
            ic.mod g
            bc.mod g
            fldvar.mod g
            constant.mod g
            indices.mod g
            is.mod g
            toleranc.mod g
            scales.mod g
            ur_facs.mod g
            leqsol.mod g
            residual.mod g
            rxns.mod g
            scalars.mod g
            compar.mod g
            parallel.mod g
            namelist.inc                                                
init_resid.$(OBJ_EXT) : init_resid.f g
            param.mod g
            param1.mod g
            physprop.mod g
            residual.mod 
iterate.$(OBJ_EXT) : iterate.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            output.mod g
            indices.mod g
            funits.mod g
            time_cpu.mod g
            pscor.mod g
            coeff.mod g
            leqsol.mod g
            visc_g.mod g
            pgcor.mod g
            cont.mod g
            scalars.mod g
            compar.mod g
            mpi_utility.mod g
            discretelement.mod g
            bc.mod g
            constant.mod 
k_epsilon_prop.$(OBJ_EXT) : k_epsilon_prop.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            drag.mod g
            run.mod g
            output.mod g
            geometry.mod g
            fldvar.mod g
            visc_g.mod g
            visc_s.mod g
            trace.mod g
            indices.mod g
            constant.mod g
            vshear.mod g
            turb.mod g
            toleranc.mod g
            compar.mod g
            tau_g.mod g
            sendrecv.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg2.inc                                                
leq_bicgs.$(OBJ_EXT) : leq_bicgs.f g
            param.mod g
            param1.mod g
            matrix.mod g
            geometry.mod g
            compar.mod g
            indices.mod g
            leqsol.mod g
            funits.mod g
            parallel.mod g
            mpi_utility.mod g
            sendrecv.mod g
            function.inc                                                
leq_gmres.$(OBJ_EXT) : leq_gmres.f g
            param.mod g
            param1.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            debug.mod g
            compar.mod g
            mpi_utility.mod g
            parallel.mod g
            funits.mod g
            gridmap.mod g
            function.inc                                                
leq_sor.$(OBJ_EXT) : leq_sor.f g
            param.mod g
            param1.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
line_too_big.$(OBJ_EXT) : line_too_big.f 
location_check.$(OBJ_EXT) : location_check.f g
            param.mod g
            param1.mod g
            funits.mod g
            geometry.mod 
location.$(OBJ_EXT) : location.f g
            param.mod g
            param1.mod 
machine.$(OBJ_EXT) : machine.f g
            machine.mod g
            param.mod g
            run.mod g
            funits.mod 
make_upper_case.$(OBJ_EXT) : make_upper_case.f 
mark_phase_4_cor.$(OBJ_EXT) : mark_phase_4_cor.f g
            param.mod g
            param1.mod g
            geometry.mod g
            indices.mod g
            fldvar.mod g
            physprop.mod g
            constant.mod g
            compar.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS2) mark_phase_4_cor.f 
mfix.$(OBJ_EXT) : mfix.f g
            param.mod g
            param1.mod g
            run.mod g
            time_cpu.mod g
            funits.mod g
            output.mod g
            compar.mod g
            mpi_utility.mod g
            parallel_mpi.mod g
            function.inc                                                
mod_bc_i.$(OBJ_EXT) : mod_bc_i.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            mpi_utility.mod g
            function.inc                                                
mod_bc_j.$(OBJ_EXT) : mod_bc_j.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            mpi_utility.mod g
            function.inc                                                
mod_bc_k.$(OBJ_EXT) : mod_bc_k.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            mpi_utility.mod g
            function.inc                                                
open_file.$(OBJ_EXT) : open_file.f 
open_files.$(OBJ_EXT) : open_files.f g
            machine.mod g
            funits.mod g
            compar.mod 
out_array_c.$(OBJ_EXT) : out_array_c.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            function.inc                                                
out_array.$(OBJ_EXT) : out_array.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            function.inc                                                
out_array_kc.$(OBJ_EXT) : out_array_kc.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            mpi_utility.mod g
            function.inc                                                
out_array_k.$(OBJ_EXT) : out_array_k.f g
            param.mod g
            param1.mod g
            geometry.mod g
            fldvar.mod g
            physprop.mod g
            indices.mod g
            funits.mod g
            compar.mod g
            function.inc                                                
out_bin_512.$(OBJ_EXT) : out_bin_512.f g
            machine.mod 
out_bin_512i.$(OBJ_EXT) : out_bin_512i.f g
            machine.mod 
out_bin_512r.$(OBJ_EXT) : out_bin_512r.f g
            machine.mod 
out_bin_r.$(OBJ_EXT) : out_bin_r.f g
            param.mod 
parse_line.$(OBJ_EXT) : parse_line.f g
            param.mod g
            param1.mod g
            parse.mod g
            compar.mod 
parse_resid_string.$(OBJ_EXT) : parse_resid_string.f g
            param.mod g
            param1.mod g
            physprop.mod g
            residual.mod g
            funits.mod g
            compar.mod 
parse_rxn.$(OBJ_EXT) : parse_rxn.f g
            param.mod g
            param1.mod g
            parse.mod g
            rxns.mod g
            compar.mod 
partial_elim.$(OBJ_EXT) : partial_elim.f g
            param.mod g
            param1.mod g
            parallel.mod g
            geometry.mod g
            matrix.mod g
            physprop.mod g
            indices.mod g
            compar.mod g
            drag.mod g
            fldvar.mod g
            run.mod g
            function.inc                                                 g
            fun_avg1.inc                                                 g
            fun_avg2.inc                                                
physical_prop.$(OBJ_EXT) : physical_prop.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            physprop.mod g
            geometry.mod g
            indices.mod g
            run.mod g
            toleranc.mod g
            constant.mod g
            compar.mod g
            funits.mod g
            cp_fun1.inc                                                  g
            function.inc                                                 g
            cp_fun2.inc                                                 
read_namelist.$(OBJ_EXT) : read_namelist.f g
            param.mod g
            param1.mod g
            run.mod g
            output.mod g
            physprop.mod g
            geometry.mod g
            ic.mod g
            is.mod g
            bc.mod g
            fldvar.mod g
            constant.mod g
            indices.mod g
            toleranc.mod g
            funits.mod g
            scales.mod g
            ur_facs.mod g
            leqsol.mod g
            residual.mod g
            rxns.mod g
            scalars.mod g
            compar.mod g
            parallel.mod g
            discretelement.mod g
            usrnlst.inc                                                  g
            namelist.inc                                                 g
            des/desnamelist.inc                                         
read_res0.$(OBJ_EXT) : read_res0.f g
            param.mod g
            param1.mod g
            geometry.mod g
            physprop.mod g
            run.mod g
            ic.mod g
            bc.mod g
            is.mod g
            constant.mod g
            funits.mod g
            output.mod g
            scales.mod g
            ur_facs.mod g
            toleranc.mod g
            leqsol.mod g
            scalars.mod g
            rxns.mod g
            compar.mod g
            mpi_utility.mod g
            fldvar.mod 
read_res1.$(OBJ_EXT) : read_res1.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            geometry.mod g
            physprop.mod g
            run.mod g
            rxns.mod g
            scalars.mod g
            funits.mod g
            energy.mod g
            compar.mod g
            mpi_utility.mod g
            sendrecv.mod 
remove_comment.$(OBJ_EXT) : remove_comment.f 
reset_new.$(OBJ_EXT) : reset_new.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            trace.mod g
            run.mod g
            scalars.mod 
rrates0.$(OBJ_EXT) : rrates0.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            rxns.mod g
            energy.mod g
            geometry.mod g
            run.mod g
            indices.mod g
            physprop.mod g
            constant.mod g
            funits.mod g
            compar.mod g
            sendrecv.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
rrates.$(OBJ_EXT) : rrates.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            rxns.mod g
            energy.mod g
            geometry.mod g
            run.mod g
            indices.mod g
            physprop.mod g
            constant.mod g
            funits.mod g
            compar.mod g
            sendrecv.mod g
            function.inc                                                
rrates_init.$(OBJ_EXT) : rrates_init.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            rxns.mod g
            energy.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
scalar_prop.$(OBJ_EXT) : scalar_prop.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            physprop.mod g
            geometry.mod g
            indices.mod g
            run.mod g
            scalars.mod g
            toleranc.mod g
            compar.mod g
            sendrecv.mod g
            function.inc                                                
seek_comment.$(OBJ_EXT) : seek_comment.f 
seek_end.$(OBJ_EXT) : seek_end.f 
set_bc0.$(OBJ_EXT) : set_bc0.f g
            param.mod g
            param1.mod g
            geometry.mod g
            compar.mod g
            mpi_utility.mod g
            physprop.mod g
            bc.mod g
            fldvar.mod g
            indices.mod g
            run.mod g
            funits.mod g
            scales.mod g
            scalars.mod g
            boundfunijk.mod g
            toleranc.mod g
            sc_p_g1.inc                                                  g
            function.inc                                                 g
            sc_p_g2.inc                                                 
set_bc1.$(OBJ_EXT) : set_bc1.f g
            param.mod g
            param1.mod g
            bc.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            run.mod g
            funits.mod g
            compar.mod g
            function.inc                                                
set_constants.$(OBJ_EXT) : set_constants.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            visc_s.mod g
            energy.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            constant.mod g
            run.mod g
            funits.mod g
            drag.mod g
            compar.mod 
set_constprop.$(OBJ_EXT) : set_constprop.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            visc_s.mod g
            visc_g.mod g
            energy.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            constant.mod g
            run.mod g
            funits.mod g
            drag.mod g
            compar.mod g
            function.inc                                                
set_flags.$(OBJ_EXT) : set_flags.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            geometry.mod g
            bc.mod g
            is.mod g
            indices.mod g
            physprop.mod g
            funits.mod g
            compar.mod g
            sendrecv.mod g
            sendrecv3.mod g
            boundfunijk.mod g
            mpi_utility.mod g
            function.inc                                                 g
            function3.inc                                               
set_fluidbed_p.$(OBJ_EXT) : set_fluidbed_p.f g
            param.mod g
            param1.mod g
            physprop.mod g
            geometry.mod g
            bc.mod g
            ic.mod g
            fldvar.mod g
            constant.mod g
            indices.mod g
            funits.mod g
            scales.mod g
            compar.mod g
            mpi_utility.mod g
            sendrecv.mod g
            sc_p_g1.inc                                                  g
            b_force1.inc                                                 g
            function.inc                                                 g
            b_force2.inc                                                 g
            sc_p_g2.inc                                                 
set_geometry1.$(OBJ_EXT) : set_geometry1.f g
            param.mod g
            param1.mod g
            parallel.mod g
            run.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
set_geometry.$(OBJ_EXT) : set_geometry.f g
            param.mod g
            param1.mod g
            run.mod g
            geometry.mod g
            compar.mod 
set_ic.$(OBJ_EXT) : set_ic.f g
            param.mod g
            param1.mod g
            geometry.mod g
            constant.mod g
            physprop.mod g
            ic.mod g
            fldvar.mod g
            visc_g.mod g
            indices.mod g
            scales.mod g
            energy.mod g
            scalars.mod g
            compar.mod g
            run.mod g
            sendrecv.mod g
            sc_p_g1.inc                                                  g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            sc_p_g2.inc                                                 
set_increments3.$(OBJ_EXT) : set_increments3.f g
            param.mod g
            param1.mod g
            indices.mod g
            geometry.mod g
            compar.mod g
            physprop.mod g
            fldvar.mod g
            funits.mod g
            function.inc                                                 g
            function3.inc                                               
set_increments.$(OBJ_EXT) : set_increments.f g
            param.mod g
            param1.mod g
            indices.mod g
            geometry.mod g
            compar.mod g
            physprop.mod g
            fldvar.mod g
            funits.mod g
            function.inc                                                
set_index1a3.$(OBJ_EXT) : set_index1a3.f g
            param.mod g
            param1.mod g
            physprop.mod g
            geometry.mod g
            compar.mod g
            fldvar.mod g
            indices.mod g
            boundfunijk3.mod g
            function.inc                                                
set_index1a.$(OBJ_EXT) : set_index1a.f g
            param.mod g
            param1.mod g
            physprop.mod g
            geometry.mod g
            compar.mod g
            fldvar.mod g
            indices.mod g
            boundfunijk.mod g
            function.inc                                                
set_index1.$(OBJ_EXT) : set_index1.f g
            param.mod g
            param1.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            constant.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
set_l_scale.$(OBJ_EXT) : set_l_scale.f g
            param.mod g
            param1.mod g
            parallel.mod g
            constant.mod g
            visc_g.mod g
            geometry.mod g
            indices.mod g
            compar.mod 
set_max2.$(OBJ_EXT) : set_max2.f g
            param.mod g
            param1.mod g
            geometry.mod g
            compar.mod 
set_mw_mix_g.$(OBJ_EXT) : set_mw_mix_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            constant.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
set_outflow.$(OBJ_EXT) : set_outflow.f g
            param.mod g
            param1.mod g
            bc.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            scalars.mod g
            run.mod g
            compar.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
set_ro_g.$(OBJ_EXT) : set_ro_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            constant.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
set_wall_bc.$(OBJ_EXT) : set_wall_bc.f g
            param.mod g
            param1.mod g
            bc.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            run.mod g
            funits.mod g
            compar.mod g
            function.inc                                                
shift_dxyz.$(OBJ_EXT) : shift_dxyz.f g
            param.mod g
            param1.mod g
            geometry.mod 
solve_continuity.$(OBJ_EXT) : solve_continuity.f g
            param.mod g
            param1.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            indices.mod g
            residual.mod g
            cont.mod g
            leqsol.mod g
            ambm.mod 
solve_energy_eq.$(OBJ_EXT) : solve_energy_eq.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            output.mod g
            indices.mod g
            drag.mod g
            residual.mod g
            ur_facs.mod g
            pgcor.mod g
            pscor.mod g
            leqsol.mod g
            bc.mod g
            energy.mod g
            rxns.mod g
            ambm.mod g
            tmp_array.mod g
            tmp_array1.mod g
            compar.mod g
            discretelement.mod g
            mflux.mod g
            radtn1.inc                                                   g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                    g
            radtn2.inc                                                  
solve_epp.$(OBJ_EXT) : solve_epp.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            geometry.mod g
            pscor.mod g
            residual.mod g
            leqsol.mod g
            physprop.mod g
            ambm.mod 
solve_granular_energy.$(OBJ_EXT) : solve_granular_energy.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            constant.mod g
            output.mod g
            indices.mod g
            drag.mod g
            residual.mod g
            ur_facs.mod g
            pgcor.mod g
            pscor.mod g
            leqsol.mod g
            bc.mod g
            energy.mod g
            rxns.mod g
            ambm.mod g
            tmp_array.mod g
            compar.mod g
            mflux.mod g
            radtn1.inc                                                   g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                    g
            radtn2.inc                                                  
solve_k_epsilon_eq.$(OBJ_EXT) : solve_k_epsilon_eq.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            output.mod g
            indices.mod g
            drag.mod g
            residual.mod g
            ur_facs.mod g
            pgcor.mod g
            pscor.mod g
            leqsol.mod g
            bc.mod g
            energy.mod g
            rxns.mod g
            turb.mod g
            usr.mod g
            ambm.mod g
            tmp_array.mod g
            compar.mod g
            mflux.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg1.inc                                                 g
            fun_avg2.inc                                                
solve_lin_eq.$(OBJ_EXT) : solve_lin_eq.f g
            param.mod g
            param1.mod g
            geometry.mod g
            compar.mod 
solve_pp_g.$(OBJ_EXT) : solve_pp_g.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            physprop.mod g
            geometry.mod g
            pgcor.mod g
            residual.mod g
            leqsol.mod g
            run.mod g
            ambm.mod 
solve_scalar_eq.$(OBJ_EXT) : solve_scalar_eq.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            output.mod g
            indices.mod g
            drag.mod g
            residual.mod g
            ur_facs.mod g
            pgcor.mod g
            pscor.mod g
            leqsol.mod g
            bc.mod g
            energy.mod g
            rxns.mod g
            scalars.mod g
            ambm.mod g
            tmp_array.mod g
            compar.mod g
            mflux.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
solve_species_eq.$(OBJ_EXT) : solve_species_eq.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            output.mod g
            indices.mod g
            drag.mod g
            residual.mod g
            ur_facs.mod g
            pgcor.mod g
            pscor.mod g
            leqsol.mod g
            bc.mod g
            energy.mod g
            rxns.mod g
            ambm.mod g
            matrix.mod g
            chischeme.mod g
            tmp_array.mod g
            compar.mod g
            mpi_utility.mod g
            sendrecv.mod g
            mflux.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
solve_vel_star.$(OBJ_EXT) : solve_vel_star.f g
            param.mod g
            param1.mod g
            toleranc.mod g
            run.mod g
            physprop.mod g
            geometry.mod g
            fldvar.mod g
            output.mod g
            indices.mod g
            drag.mod g
            residual.mod g
            ur_facs.mod g
            pgcor.mod g
            pscor.mod g
            leqsol.mod g
            ambm.mod g
            tmp_array1.mod g
            tmp_array.mod g
            compar.mod g
            discretelement.mod 
source_granular_energy.$(OBJ_EXT) : source_granular_energy.f g
            param.mod g
            param1.mod g
            parallel.mod g
            physprop.mod g
            run.mod g
            drag.mod g
            geometry.mod g
            fldvar.mod g
            visc_g.mod g
            visc_s.mod g
            trace.mod g
            turb.mod g
            indices.mod g
            constant.mod g
            toleranc.mod g
            compar.mod g
            s_pr1.inc                                                    g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg2.inc                                                 g
            s_pr2.inc                                                   
source_phi.$(OBJ_EXT) : source_phi.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            tau_s.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
source_pp_g.$(OBJ_EXT) : source_pp_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            physprop.mod g
            fldvar.mod g
            rxns.mod g
            run.mod g
            geometry.mod g
            indices.mod g
            pgcor.mod g
            bc.mod g
            vshear.mod g
            xsi_array.mod g
            compar.mod g
            ur_facs.mod g
            function.inc                                                
source_rop_g.$(OBJ_EXT) : source_rop_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            fldvar.mod g
            rxns.mod g
            run.mod g
            geometry.mod g
            indices.mod g
            pgcor.mod g
            compar.mod g
            function.inc                                                
source_rop_s.$(OBJ_EXT) : source_rop_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            fldvar.mod g
            rxns.mod g
            run.mod g
            geometry.mod g
            indices.mod g
            pgcor.mod g
            pscor.mod g
            compar.mod g
            function.inc                                                
source_u_g.$(OBJ_EXT) : source_u_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_g.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            tau_g.mod g
            bc.mod g
            compar.mod g
            sendrecv.mod g
            output.mod g
            turb.mod g
            mpi_utility.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
source_u_s.$(OBJ_EXT) : source_u_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            tau_s.mod g
            bc.mod g
            compar.mod g
            sendrecv.mod g
            output.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
source_v_g.$(OBJ_EXT) : source_v_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_g.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            tau_g.mod g
            bc.mod g
            vshear.mod g
            compar.mod g
            sendrecv.mod g
            output.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
source_v_s.$(OBJ_EXT) : source_v_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            tau_s.mod g
            bc.mod g
            vshear.mod g
            compar.mod g
            sendrecv.mod g
            output.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
source_w_g.$(OBJ_EXT) : source_w_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_g.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            tau_g.mod g
            bc.mod g
            compar.mod g
            sendrecv.mod g
            output.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
source_w_s.$(OBJ_EXT) : source_w_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            tau_s.mod g
            bc.mod g
            compar.mod g
            sendrecv.mod g
            output.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
tau_u_g.$(OBJ_EXT) : tau_u_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_g.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            compar.mod g
            sendrecv.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
tau_u_s.$(OBJ_EXT) : tau_u_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            vshear.mod g
            sendrecv.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
tau_v_g.$(OBJ_EXT) : tau_v_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_g.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            sendrecv.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
tau_v_s.$(OBJ_EXT) : tau_v_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            sendrecv.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
tau_w_g.$(OBJ_EXT) : tau_w_g.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_g.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            sendrecv.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
tau_w_s.$(OBJ_EXT) : tau_w_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_s.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            sendrecv.mod g
            compar.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
test_lin_eq.$(OBJ_EXT) : test_lin_eq.f g
            param.mod g
            param1.mod g
            matrix.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            function.inc                                                
time_march.$(OBJ_EXT) : time_march.f g
            param.mod g
            param1.mod g
            run.mod g
            output.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            pgcor.mod g
            pscor.mod g
            cont.mod g
            coeff.mod g
            tau_g.mod g
            tau_s.mod g
            visc_g.mod g
            visc_s.mod g
            funits.mod g
            vshear.mod g
            scalars.mod g
            drag.mod g
            rxns.mod g
            compar.mod g
            time_cpu.mod g
            discretelement.mod 
transfer.$(OBJ_EXT) : transfer.f g
            param.mod g
            param1.mod g
            geometry.mod g
            indices.mod 
transport_prop.$(OBJ_EXT) : transport_prop.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            physprop.mod g
            geometry.mod g
            indices.mod g
            run.mod g
            toleranc.mod g
            compar.mod 
undef_2_0.$(OBJ_EXT) : undef_2_0.f g
            param.mod g
            param1.mod g
            geometry.mod g
            compar.mod 
under_relax.$(OBJ_EXT) : under_relax.f g
            param.mod g
            param1.mod g
            geometry.mod g
            indices.mod g
            compar.mod g
            sendrecv.mod g
            function.inc                                                
update_old.$(OBJ_EXT) : update_old.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            run.mod g
            trace.mod g
            visc_s.mod g
            scalars.mod 
usr0.$(OBJ_EXT) : usr0.f g
            usr.mod 
usr1.$(OBJ_EXT) : usr1.f g
            usr.mod 
usr2.$(OBJ_EXT) : usr2.f g
            usr.mod 
usr3.$(OBJ_EXT) : usr3.f g
            usr.mod 
usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f g
            usrnlst.inc                                                 
usr_write_out0.$(OBJ_EXT) : usr_write_out0.f 
usr_write_out1.$(OBJ_EXT) : usr_write_out1.f 
vavg_u_g.$(OBJ_EXT) : vavg_u_g.f g
            param.mod g
            param1.mod g
            run.mod g
            parallel.mod g
            fldvar.mod g
            bc.mod g
            geometry.mod g
            physprop.mod g
            indices.mod g
            compar.mod g
            mpi_utility.mod g
            function.inc                                                
vavg_u_s.$(OBJ_EXT) : vavg_u_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            bc.mod g
            geometry.mod g
            physprop.mod g
            indices.mod g
            compar.mod g
            mpi_utility.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
vavg_v_g.$(OBJ_EXT) : vavg_v_g.f g
            param.mod g
            param1.mod g
            run.mod g
            parallel.mod g
            fldvar.mod g
            bc.mod g
            geometry.mod g
            physprop.mod g
            indices.mod g
            compar.mod g
            mpi_utility.mod g
            function.inc                                                
vavg_v_s.$(OBJ_EXT) : vavg_v_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            bc.mod g
            geometry.mod g
            physprop.mod g
            indices.mod g
            compar.mod g
            mpi_utility.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
vavg_w_g.$(OBJ_EXT) : vavg_w_g.f g
            param.mod g
            param1.mod g
            run.mod g
            parallel.mod g
            fldvar.mod g
            bc.mod g
            geometry.mod g
            physprop.mod g
            indices.mod g
            compar.mod g
            mpi_utility.mod g
            function.inc                                                
vavg_w_s.$(OBJ_EXT) : vavg_w_s.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            bc.mod g
            geometry.mod g
            physprop.mod g
            indices.mod g
            compar.mod g
            mpi_utility.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
vf_gs_x.$(OBJ_EXT) : vf_gs_x.f g
            param.mod g
            param1.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            compar.mod g
            drag.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                
vf_gs_y.$(OBJ_EXT) : vf_gs_y.f g
            param.mod g
            param1.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            compar.mod g
            drag.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                
vf_gs_z.$(OBJ_EXT) : vf_gs_z.f g
            param.mod g
            param1.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            compar.mod g
            drag.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                
write_ab_m.$(OBJ_EXT) : write_ab_m.f g
            param.mod g
            param1.mod g
            matrix.mod g
            compar.mod g
            mpi_utility.mod g
            indices.mod g
            function.inc                                                
write_ab_m_var.$(OBJ_EXT) : write_ab_m_var.f g
            param.mod g
            param1.mod g
            matrix.mod g
            geometry.mod g
            compar.mod g
            mpi_utility.mod g
            indices.mod g
            function.inc                                                
write_error.$(OBJ_EXT) : write_error.f g
            param.mod g
            param1.mod g
            funits.mod 
write_header.$(OBJ_EXT) : write_header.f g
            param.mod g
            param1.mod g
            run.mod g
            output.mod g
            funits.mod g
            compar.mod 
write_out0.$(OBJ_EXT) : write_out0.f g
            param.mod g
            param1.mod g
            run.mod g
            output.mod g
            physprop.mod g
            geometry.mod g
            ic.mod g
            bc.mod g
            is.mod g
            fldvar.mod g
            constant.mod g
            indices.mod g
            funits.mod g
            toleranc.mod g
            scales.mod g
            scalars.mod g
            ur_facs.mod g
            leqsol.mod g
            compar.mod g
            mpi_utility.mod g
            sendrecv.mod g
            function.inc                                                
write_out1.$(OBJ_EXT) : write_out1.f g
            param.mod g
            param1.mod g
            physprop.mod g
            fldvar.mod g
            run.mod g
            scalars.mod g
            funits.mod g
            rxns.mod g
            compar.mod g
            mpi_utility.mod 
write_out3.$(OBJ_EXT) : write_out3.f g
            funits.mod g
            compar.mod 
write_res0.$(OBJ_EXT) : write_res0.f g
            param.mod g
            param1.mod g
            geometry.mod g
            physprop.mod g
            run.mod g
            ic.mod g
            is.mod g
            bc.mod g
            constant.mod g
            funits.mod g
            output.mod g
            scales.mod g
            scalars.mod g
            rxns.mod g
            ur_facs.mod g
            leqsol.mod g
            toleranc.mod g
            compar.mod g
            mpi_utility.mod g
            sendrecv.mod 
write_res1.$(OBJ_EXT) : write_res1.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            geometry.mod g
            physprop.mod g
            run.mod g
            scalars.mod g
            rxns.mod g
            funits.mod g
            output.mod g
            energy.mod g
            compar.mod g
            mpi_utility.mod g
            sendrecv.mod 
write_spx0.$(OBJ_EXT) : write_spx0.f g
            param.mod g
            param1.mod g
            run.mod g
            funits.mod g
            compar.mod g
            mpi_utility.mod 
write_spx1.$(OBJ_EXT) : write_spx1.f g
            param.mod g
            param1.mod g
            fldvar.mod g
            geometry.mod g
            physprop.mod g
            run.mod g
            funits.mod g
            scalars.mod g
            output.mod g
            rxns.mod g
            compar.mod g
            mpi_utility.mod g
            sendrecv.mod 
write_table.$(OBJ_EXT) : write_table.f g
            param.mod g
            param1.mod g
            funits.mod 
write_usr0.$(OBJ_EXT) : write_usr0.f 
write_usr1.$(OBJ_EXT) : write_usr1.f 
xerbla.$(OBJ_EXT) : xerbla.f g
            compar.mod 
zero_array.$(OBJ_EXT) : zero_array.f g
            param.mod g
            param1.mod 
zero_norm_vel.$(OBJ_EXT) : zero_norm_vel.f g
            param.mod g
            param1.mod g
            parallel.mod g
            geometry.mod g
            physprop.mod g
            fldvar.mod g
            indices.mod g
            is.mod g
            compar.mod g
            function.inc                                                
add_part_to_link_list.$(OBJ_EXT) : ./cohesion/add_part_to_link_list.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/add_part_to_link_list.f 
calc_app_coh_force.$(OBJ_EXT) : ./cohesion/calc_app_coh_force.f g
            discretelement.mod g
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_app_coh_force.f 
calc_cap_coh_force.$(OBJ_EXT) : ./cohesion/calc_cap_coh_force.f g
            discretelement.mod g
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_cap_coh_force.f 
calc_cohesive_forces.$(OBJ_EXT) : ./cohesion/calc_cohesive_forces.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_cohesive_forces.f 
calc_esc_coh_force.$(OBJ_EXT) : ./cohesion/calc_esc_coh_force.f g
            discretelement.mod g
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_esc_coh_force.f 
calc_square_well.$(OBJ_EXT) : ./cohesion/calc_square_well.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_square_well.f 
calc_van_der_waals.$(OBJ_EXT) : ./cohesion/calc_van_der_waals.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_van_der_waals.f 
check_link.$(OBJ_EXT) : ./cohesion/check_link.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_link.f 
check_sw_wall_interaction.$(OBJ_EXT) : ./cohesion/check_sw_wall_interaction.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_sw_wall_interaction.f 
check_vdw_wall_interaction.$(OBJ_EXT) : ./cohesion/check_vdw_wall_interaction.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_vdw_wall_interaction.f 
initialize_cohesion_parameters.$(OBJ_EXT) : ./cohesion/initialize_cohesion_parameters.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/initialize_cohesion_parameters.f 
initialize_coh_int_search.$(OBJ_EXT) : ./cohesion/initialize_coh_int_search.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/initialize_coh_int_search.f 
linked_interaction_eval.$(OBJ_EXT) : ./cohesion/linked_interaction_eval.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/linked_interaction_eval.f 
remove_part_from_link_list.$(OBJ_EXT) : ./cohesion/remove_part_from_link_list.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/remove_part_from_link_list.f 
unlinked_interaction_eval.$(OBJ_EXT) : ./cohesion/unlinked_interaction_eval.f g
            discretelement.mod g
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/unlinked_interaction_eval.f 
update_search_grids.$(OBJ_EXT) : ./cohesion/update_search_grids.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/update_search_grids.f 
calc_force_des.$(OBJ_EXT) : ./des/calc_force_des.f g
            discretelement.mod g
            geometry.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_force_des.f 
cfassign.$(OBJ_EXT) : ./des/cfassign.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfassign.f 
cffctow.$(OBJ_EXT) : ./des/cffctow.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffctow.f 
cffn.$(OBJ_EXT) : ./des/cffn.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffn.f 
cffnwall.$(OBJ_EXT) : ./des/cffnwall.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffnwall.f 
cfft.$(OBJ_EXT) : ./des/cfft.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfft.f 
cfftwall.$(OBJ_EXT) : ./des/cfftwall.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfftwall.f 
cfincrementaloverlaps.$(OBJ_EXT) : ./des/cfincrementaloverlaps.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfincrementaloverlaps.f 
cfnewvalues.$(OBJ_EXT) : ./des/cfnewvalues.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            compar.mod g
            sendrecv.mod g
            physprop.mod g
            fldvar.mod g
            visc_g.mod g
            rxns.mod g
            run.mod g
            geometry.mod g
            indices.mod g
            drag.mod g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfnewvalues.f 
cfnocontact.$(OBJ_EXT) : ./des/cfnocontact.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfnocontact.f 
cfnormal.$(OBJ_EXT) : ./des/cfnormal.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfnormal.f 
cfoutofbox.$(OBJ_EXT) : ./des/cfoutofbox.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfoutofbox.f 
cfperiodicwallneighbourx.$(OBJ_EXT) : ./des/cfperiodicwallneighbourx.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwallneighbourx.f 
cfperiodicwallneighboury.$(OBJ_EXT) : ./des/cfperiodicwallneighboury.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwallneighboury.f 
cfperiodicwallneighbourz.$(OBJ_EXT) : ./des/cfperiodicwallneighbourz.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwallneighbourz.f 
cfperiodicwallx.$(OBJ_EXT) : ./des/cfperiodicwallx.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwallx.f 
cfperiodicwally.$(OBJ_EXT) : ./des/cfperiodicwally.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwally.f 
cfperiodicwallz.$(OBJ_EXT) : ./des/cfperiodicwallz.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwallz.f 
cfrelvel.$(OBJ_EXT) : ./des/cfrelvel.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfrelvel.f 
cfslide.$(OBJ_EXT) : ./des/cfslide.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslide.f 
cfslidewall.$(OBJ_EXT) : ./des/cfslidewall.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslidewall.f 
cftangent.$(OBJ_EXT) : ./des/cftangent.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cftangent.f 
cftotaloverlaps.$(OBJ_EXT) : ./des/cftotaloverlaps.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cftotaloverlaps.f 
cfupdateold.$(OBJ_EXT) : ./des/cfupdateold.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfupdateold.f 
cfvrn.$(OBJ_EXT) : ./des/cfvrn.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfvrn.f 
cfvrt.$(OBJ_EXT) : ./des/cfvrt.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfvrt.f 
cfwallcontact.$(OBJ_EXT) : ./des/cfwallcontact.f g
            discretelement.mod g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            run.mod g
            geometry.mod g
            matrix.mod g
            indices.mod g
            physprop.mod g
            drag.mod g
            constant.mod g
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallcontact.f 
cfwallposvel.$(OBJ_EXT) : ./des/cfwallposvel.f g
            discretelement.mod g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            run.mod g
            geometry.mod g
            matrix.mod g
            indices.mod g
            physprop.mod g
            drag.mod g
            constant.mod g
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallposvel.f 
des_calc_d.$(OBJ_EXT) : ./des/des_calc_d.f g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            geometry.mod g
            indices.mod g
            physprop.mod g
            run.mod g
            scales.mod g
            compar.mod g
            sendrecv.mod g
            discretelement.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_calc_d.f 
des_granular_temperature.$(OBJ_EXT) : ./des/des_granular_temperature.f g
            discretelement.mod g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            run.mod g
            geometry.mod g
            matrix.mod g
            indices.mod g
            physprop.mod g
            drag.mod g
            constant.mod g
            compar.mod g
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_granular_temperature.f 
des_init_namelist.$(OBJ_EXT) : ./des/des_init_namelist.f g
            discretelement.mod g
            des/desnamelist.inc                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_namelist.f 
des_inlet_outlet.$(OBJ_EXT) : ./des/des_inlet_outlet.f g
            discretelement.mod g
            geometry.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_inlet_outlet.f 
des_time_march.$(OBJ_EXT) : ./des/des_time_march.f g
            param.mod g
            param1.mod g
            run.mod g
            output.mod g
            physprop.mod g
            fldvar.mod g
            geometry.mod g
            pgcor.mod g
            pscor.mod g
            cont.mod g
            coeff.mod g
            tau_g.mod g
            tau_s.mod g
            visc_g.mod g
            visc_s.mod g
            funits.mod g
            vshear.mod g
            scalars.mod g
            drag.mod g
            rxns.mod g
            compar.mod g
            time_cpu.mod g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_time_march.f 
drag_fgs.$(OBJ_EXT) : ./des/drag_fgs.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_g.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            tau_g.mod g
            bc.mod g
            compar.mod g
            sendrecv.mod g
            discretelement.mod g
            drag.mod g
            function.inc                                                 g
            ep_s1.inc                                                    g
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/drag_fgs.f 
gas_drag.$(OBJ_EXT) : ./des/gas_drag.f g
            param.mod g
            param1.mod g
            parallel.mod g
            matrix.mod g
            scales.mod g
            constant.mod g
            physprop.mod g
            fldvar.mod g
            visc_g.mod g
            rxns.mod g
            run.mod g
            toleranc.mod g
            geometry.mod g
            indices.mod g
            is.mod g
            tau_g.mod g
            bc.mod g
            compar.mod g
            sendrecv.mod g
            discretelement.mod g
            drag.mod g
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/gas_drag.f 
make_arrays_des.$(OBJ_EXT) : ./des/make_arrays_des.f g
            funits.mod g
            compar.mod g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/make_arrays_des.f 
neighbour.$(OBJ_EXT) : ./des/neighbour.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/neighbour.f 
nsquare.$(OBJ_EXT) : ./des/nsquare.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/nsquare.f 
octree.$(OBJ_EXT) : ./des/octree.f g
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/octree.f 
particles_in_cell.$(OBJ_EXT) : ./des/particles_in_cell.f g
            discretelement.mod g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            run.mod g
            geometry.mod g
            matrix.mod g
            indices.mod g
            physprop.mod g
            drag.mod g
            constant.mod g
            compar.mod g
            sendrecv.mod g
            function.inc                                                 g
            ep_s1.inc                                                    g
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/particles_in_cell.f 
periodic_wall_calc_force_des.$(OBJ_EXT) : ./des/periodic_wall_calc_force_des.f g
            discretelement.mod g
            geometry.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/periodic_wall_calc_force_des.f 
pressure_drop.$(OBJ_EXT) : ./des/pressure_drop.f g
            discretelement.mod g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            run.mod g
            geometry.mod g
            matrix.mod g
            indices.mod g
            physprop.mod g
            drag.mod g
            constant.mod g
            compar.mod g
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/pressure_drop.f 
print_vel.$(OBJ_EXT) : ./des/print_vel.f g
            discretelement.mod g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            run.mod g
            geometry.mod g
            matrix.mod g
            indices.mod g
            physprop.mod g
            drag.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/print_vel.f 
quadtree.$(OBJ_EXT) : ./des/quadtree.f g
            discretelement.mod g
            param.mod g
            param1.mod g
            parallel.mod g
            fldvar.mod g
            run.mod g
            geometry.mod g
            matrix.mod g
            indices.mod g
            physprop.mod g
            drag.mod g
            constant.mod g
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/quadtree.f 
