.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
mfix.exe : \
    AMBM.mod \
    BC.mod \
    BOUNDFUNIJK3.mod \
    BOUNDFUNIJK.mod \
    CHECK.mod \
    CHISCHEME.mod \
    COEFF.mod \
    CONSTANT.mod \
    CONT.mod \
    CORNER.mod \
    DRAG.mod \
    ENERGY.mod \
    FLDVAR.mod \
    FUNCTION.mod \
    FUNITS.mod \
    GEOMETRY.mod \
    IC.mod \
    INDICES.mod \
    IS.mod \
    LEQSOL.mod \
    MACHINE.mod \
    MATRIX.mod \
    MFLUX.mod \
    OUTPUT.mod \
    PARALLEL.mod \
    PARAM1.mod \
    PARAM.mod \
    PARSE.mod \
    PGCOR.mod \
    PHYSPROP.mod \
    PSCOR.mod \
    RESIDUAL.mod \
    RUN.mod \
    RXNS.mod \
    SCALARS.mod \
    SCALES.mod \
    TAU_G.mod \
    TAU_S.mod \
    TIME_CPU.mod \
    TMP_ARRAY1.mod \
    TMP_ARRAY.mod \
    TOLERANC.mod \
    TRACE.mod \
    TURB.mod \
    UR_FACS.mod \
    USR.mod \
    VISC_G.mod \
    VISC_S.mod \
    VSHEAR.mod \
    XSI_ARRAY.mod \
    DISCRETELEMENT.mod \
    COMPAR.mod \
    DBG_UTIL.mod \
    DEBUG.mod \
    GRIDMAP.mod \
    MPI.mod \
    MPI_UTILITY.mod \
    PARALLEL_MPI.mod \
    SENDRECV3.mod \
    SENDRECV.mod \
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
BOUNDFUNIJK3.mod : boundfunijk3_mod.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk3_mod.f 
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
CHECK.mod : check_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_mod.f 
CHISCHEME.mod : chischeme_mod.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) chischeme_mod.f 
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
MFLUX.mod : mflux_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) mflux_mod.f 
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
TURB.mod : turb_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) turb_mod.f 
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
DISCRETELEMENT.mod : ./des/discretelement_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/discretelement_mod.f 
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
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/debug_mod.f 
GRIDMAP.mod : ./dmp_modules/gridmap_mod.f \
            MPI_UTILITY.mod \
            PARALLEL_MPI.mod \
            GEOMETRY.mod \
            SENDRECV.mod \
            COMPAR.mod \
            RUN.mod \
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
SENDRECV3.mod : ./dmp_modules/sendrecv3_mod.f \
            PARALLEL_MPI.mod \
            DEBUG.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            MPI.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv3_mod.f 
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
ar             .$(OBJ_EXT) : tant           .f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            MFLUX.mod g
            PHYSPROP.mod g
            RUN.mod g
            PARALLEL.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            SENDRECV.mod g
            XSI_ARRAY.mod g
            MPI_UTILITY.mod g
            MATRIX.mod g
            TOLERANC.mod g
            SENDRECV3.mod g
            TMP_ARRAY.mod g
            VSHEAR.mod g
            SCALES.mod g
            CONSTANT.mod g
            VISC_S.mod g
            OUTPUT.mod g
            IS.mod g
            VISC_G.mod g
            PGCOR.mod g
            UR_FACS.mod g
            FUNITS.mod g
            TIME_CPU.mod g
            PSCOR.mod g
            COEFF.mod g
            LEQSOL.mod g
            CONT.mod g
            SCALARS.mod g
            DISCRETELEMENT.mod g
            BC.mod g
            PARAM1.mod g
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
            NT.mod g
            .mod g
            M.mod g
            M1.mod g
            AR.mod g
            X.mod g
            PROP.mod g
            .mod g
            LLEL.mod g
            ETRY.mod g
            CES.mod g
            AR.mod g
            RECV.mod g
            ARRAY.mod g
            UTILITY.mod g
            IX.mod g
            RANC.mod g
            RECV3.mod g
            ARRAY.mod g
            AR.mod g
            ES.mod g
            TANT.mod g
            _S.mod g
            UT.mod g
            .mod g
            _G.mod g
            R.mod g
            ACS.mod g
            TS.mod g
            _CPU.mod g
            R.mod g
            F.mod g
            OL.mod g
            .mod g
            ARS.mod g
            RETELEMENT.mod g
            M1.mod g
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
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            GEOMETRY.mod g
            OUTPUT.mod g
            INDICES.mod g
            BC.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
bc_theta.$(OBJ_EXT) : bc_theta.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            GEOMETRY.mod g
            OUTPUT.mod g
            INDICES.mod g
            BC.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            TURB.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
b_m_p_star.$(OBJ_EXT) : b_m_p_star.f g
            PARAM.mod g
            PARAM1.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            RUN.mod g
            RXNS.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
bound_x.$(OBJ_EXT) : bound_x.f g
            PARAM.mod g
            PARAM1.mod 
calc_cell.$(OBJ_EXT) : calc_cell.f g
            PARAM.mod g
            PARAM1.mod 
calc_coeff.$(OBJ_EXT) : calc_coeff.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            RXNS.mod g
            FUNITS.mod g
            COMPAR.mod 
calc_d.$(OBJ_EXT) : calc_d.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            RUN.mod g
            SCALES.mod g
            COMPAR.mod g
            SENDRECV.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
calc_dif_g.$(OBJ_EXT) : calc_dif_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            CONSTANT.mod g
            COMPAR.mod g
            SENDRECV.mod g
            RUN.mod g
            function.inc                                                
calc_dif_s.$(OBJ_EXT) : calc_dif_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            CONSTANT.mod g
            TOLERANC.mod g
            COMPAR.mod g
            SENDRECV.mod g
            RUN.mod g
            function.inc                                                
calc_drag.$(OBJ_EXT) : calc_drag.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            RUN.mod g
            DRAG.mod g
            COMPAR.mod g
            DISCRETELEMENT.mod 
calc_e.$(OBJ_EXT) : calc_e.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            RUN.mod g
            CONSTANT.mod g
            COMPAR.mod g
            SENDRECV.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
calc_gama.$(OBJ_EXT) : calc_gama.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            ENERGY.mod g
            RXNS.mod g
            INDICES.mod g
            COMPAR.mod g
            SENDRECV.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
calc_grbdry.$(OBJ_EXT) : calc_grbdry.f g
            PARAM.mod g
            PARAM1.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            RUN.mod g
            TURB.mod g
            VISC_S.mod g
            GEOMETRY.mod g
            INDICES.mod g
            BC.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
calc_k_cp.$(OBJ_EXT) : calc_k_cp.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            INDICES.mod g
            PSCOR.mod g
            GEOMETRY.mod g
            CONSTANT.mod g
            RUN.mod g
            VISC_S.mod g
            TRACE.mod g
            COMPAR.mod g
            SENDRECV.mod g
            ep_s1.inc                                                    g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            ep_s2.inc                                                   
calc_k_g.$(OBJ_EXT) : calc_k_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            CONSTANT.mod g
            COMPAR.mod g
            RUN.mod g
            SENDRECV.mod g
            function.inc                                                
calc_k_s.$(OBJ_EXT) : calc_k_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            CONSTANT.mod g
            TOLERANC.mod g
            COMPAR.mod g
            SENDRECV.mod g
            RUN.mod g
            function.inc                                                
calc_mflux.$(OBJ_EXT) : calc_mflux.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            MFLUX.mod g
            PHYSPROP.mod g
            RUN.mod g
            PARALLEL.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
calc_mu_g.$(OBJ_EXT) : calc_mu_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            VISC_G.mod g
            VISC_S.mod g
            INDICES.mod g
            CONSTANT.mod g
            TOLERANC.mod g
            COMPAR.mod g
            DRAG.mod g
            RUN.mod g
            TURB.mod g
            SENDRECV.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg2.inc                                                
calc_mu_s.$(OBJ_EXT) : calc_mu_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            DRAG.mod g
            RUN.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            VISC_G.mod g
            VISC_S.mod g
            TRACE.mod g
            TURB.mod g
            INDICES.mod g
            CONSTANT.mod g
            TOLERANC.mod g
            VSHEAR.mod g
            COMPAR.mod g
            SENDRECV.mod g
            s_pr1.inc                                                    g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg2.inc                                                 g
            s_pr2.inc                                                   
calc_mw.$(OBJ_EXT) : calc_mw.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod 
calc_outflow.$(OBJ_EXT) : calc_outflow.f g
            PARAM.mod g
            PARAM1.mod g
            BC.mod g
            FLDVAR.mod g
            INDICES.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
calc_p_star.$(OBJ_EXT) : calc_p_star.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            CONSTANT.mod g
            PGCOR.mod g
            PSCOR.mod g
            UR_FACS.mod g
            RESIDUAL.mod g
            COMPAR.mod g
            FLDVAR.mod g
            TOLERANC.mod g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
calc_resid.$(OBJ_EXT) : calc_resid.f g
            PARAM.mod g
            PARAM1.mod g
            MATRIX.mod g
            PARALLEL.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            FLDVAR.mod g
            RUN.mod g
            BC.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            RESIDUAL.mod g
            RXNS.mod g
            function.inc                                                
calc_s_ddot_s.$(OBJ_EXT) : calc_s_ddot_s.f g
            PARAM.mod g
            PARAM1.mod g
            CONSTANT.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                
calc_trd_g.$(OBJ_EXT) : calc_trd_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            INDICES.mod g
            COMPAR.mod g
            SENDRECV.mod g
            function.inc                                                
calc_trd_s.$(OBJ_EXT) : calc_trd_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            INDICES.mod g
            PHYSPROP.mod g
            COMPAR.mod g
            SENDRECV.mod g
            function.inc                                                
calc_u_friction.$(OBJ_EXT) : calc_u_friction.f g
            PARAM.mod g
            PARAM1.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            RUN.mod g
            TURB.mod g
            VISC_S.mod g
            GEOMETRY.mod g
            INDICES.mod g
            BC.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
calc_vol_fr.$(OBJ_EXT) : calc_vol_fr.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            CONSTANT.mod g
            PGCOR.mod g
            PSCOR.mod g
            COMPAR.mod g
            SENDRECV.mod g
            ep_s1.inc                                                    g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS2) calc_vol_fr.f 
calc_xsi.$(OBJ_EXT) : calc_xsi.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            GEOMETRY.mod g
            INDICES.mod g
            VSHEAR.mod g
            CHISCHEME.mod g
            COMPAR.mod g
            SENDRECV.mod g
            xsi1.inc                                                     g
            function.inc                                                 g
            xsi2.inc                                                    
cal_d.$(OBJ_EXT) : cal_d.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            RXNS.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            TAU_S.mod g
            BC.mod g
            VSHEAR.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
check_ab_m.$(OBJ_EXT) : check_ab_m.f g
            PARAM.mod g
            PARAM1.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
check_convergence.$(OBJ_EXT) : check_convergence.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            RUN.mod g
            RESIDUAL.mod g
            TOLERANC.mod g
            MPI_UTILITY.mod 
check_data_01.$(OBJ_EXT) : check_data_01.f g
            PARAM.mod g
            PARAM1.mod g
            CONSTANT.mod g
            RUN.mod g
            PHYSPROP.mod g
            INDICES.mod g
            SCALARS.mod g
            FUNITS.mod 
check_data_02.$(OBJ_EXT) : check_data_02.f g
            PARAM.mod g
            PARAM1.mod g
            OUTPUT.mod g
            LEQSOL.mod g
            GEOMETRY.mod g
            RUN.mod g
            RXNS.mod 
check_data_03.$(OBJ_EXT) : check_data_03.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            BC.mod g
            FUNITS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod 
check_data_04.$(OBJ_EXT) : check_data_04.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            INDICES.mod g
            PHYSPROP.mod g
            CONSTANT.mod g
            FUNITS.mod 
check_data_05.$(OBJ_EXT) : check_data_05.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            FUNITS.mod g
            RUN.mod g
            INDICES.mod 
check_data_06.$(OBJ_EXT) : check_data_06.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            IC.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            RUN.mod g
            INDICES.mod g
            FUNITS.mod g
            SCALARS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            SENDRECV.mod g
            function.inc                                                
check_data_07.$(OBJ_EXT) : check_data_07.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            RUN.mod g
            BC.mod g
            INDICES.mod g
            FUNITS.mod g
            SCALARS.mod g
            COMPAR.mod g
            SENDRECV.mod g
            function.inc                                                
check_data_08.$(OBJ_EXT) : check_data_08.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            RUN.mod g
            IS.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            function.inc                                                
check_data_09.$(OBJ_EXT) : check_data_09.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            RUN.mod g
            RXNS.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod 
check_data_20.$(OBJ_EXT) : check_data_20.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            FLDVAR.mod g
            RUN.mod g
            GEOMETRY.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            INDICES.mod g
            FUNITS.mod g
            VISC_G.mod g
            RXNS.mod g
            SCALARS.mod g
            COMPAR.mod g
            SENDRECV.mod g
            function.inc                                                
check_data_30.$(OBJ_EXT) : check_data_30.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            FLDVAR.mod g
            RXNS.mod g
            VISC_S.mod g
            VISC_G.mod g
            GEOMETRY.mod g
            RUN.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            function.inc                                                
check_mass_balance.$(OBJ_EXT) : check_mass_balance.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            FLDVAR.mod g
            RXNS.mod g
            GEOMETRY.mod g
            RUN.mod g
            BC.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            OUTPUT.mod g
            CHECK.mod g
            PARALLEL.mod g
            MATRIX.mod g
            function.inc                                                
check_one_axis.$(OBJ_EXT) : check_one_axis.f g
            PARAM.mod g
            PARAM1.mod g
            FUNITS.mod 
check_plane.$(OBJ_EXT) : check_plane.f g
            FUNITS.mod g
            COMPAR.mod 
cn_extrapol.$(OBJ_EXT) : cn_extrapol.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            SCALARS.mod g
            TRACE.mod g
            RUN.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            function.inc                                                
compare.$(OBJ_EXT) : compare.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
conv_dif_phi.$(OBJ_EXT) : conv_dif_phi.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            GEOMETRY.mod g
            COMPAR.mod g
            SENDRECV.mod g
            XSI_ARRAY.mod g
            MPI_UTILITY.mod g
            INDICES.mod g
            PARALLEL.mod g
            MATRIX.mod g
            TOLERANC.mod g
            SENDRECV3.mod g
            TMP_ARRAY.mod g
            VSHEAR.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            OUTPUT.mod g
            IS.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            function3.inc                                                g
            ep_s1.inc                                                    g
            ep_s2.inc                                                   
conv_dif_u_g.$(OBJ_EXT) : conv_dif_u_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            RUN.mod g
            VISC_G.mod g
            COMPAR.mod g
            TOLERANC.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            MFLUX.mod g
            VSHEAR.mod g
            XSI_ARRAY.mod g
            TMP_ARRAY.mod g
            SENDRECV.mod g
            SENDRECV3.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_dif_u_s.$(OBJ_EXT) : conv_dif_u_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            RUN.mod g
            PHYSPROP.mod g
            VISC_S.mod g
            COMPAR.mod g
            TOLERANC.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            MFLUX.mod g
            XSI_ARRAY.mod g
            TMP_ARRAY.mod g
            SENDRECV.mod g
            SENDRECV3.mod g
            VSHEAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_dif_v_g.$(OBJ_EXT) : conv_dif_v_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            RUN.mod g
            VISC_G.mod g
            COMPAR.mod g
            TOLERANC.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            MFLUX.mod g
            XSI_ARRAY.mod g
            VSHEAR.mod g
            TMP_ARRAY.mod g
            SENDRECV.mod g
            SENDRECV3.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_dif_v_s.$(OBJ_EXT) : conv_dif_v_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            RUN.mod g
            PHYSPROP.mod g
            VISC_S.mod g
            COMPAR.mod g
            TOLERANC.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            MFLUX.mod g
            XSI_ARRAY.mod g
            TMP_ARRAY.mod g
            SENDRECV.mod g
            SENDRECV3.mod g
            VSHEAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_dif_w_g.$(OBJ_EXT) : conv_dif_w_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            RUN.mod g
            VISC_G.mod g
            COMPAR.mod g
            TOLERANC.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            MFLUX.mod g
            XSI_ARRAY.mod g
            TMP_ARRAY.mod g
            SENDRECV.mod g
            SENDRECV3.mod g
            VSHEAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_dif_w_s.$(OBJ_EXT) : conv_dif_w_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            RUN.mod g
            PHYSPROP.mod g
            VISC_S.mod g
            COMPAR.mod g
            TOLERANC.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            MFLUX.mod g
            XSI_ARRAY.mod g
            TMP_ARRAY.mod g
            SENDRECV.mod g
            SENDRECV3.mod g
            VSHEAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            function3.inc                                               
conv_pp_g.$(OBJ_EXT) : conv_pp_g.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            RUN.mod g
            PARALLEL.mod g
            MATRIX.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PGCOR.mod g
            COMPAR.mod g
            MFLUX.mod g
            function.inc                                                
conv_rop.$(OBJ_EXT) : conv_rop.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            MFLUX.mod g
            PHYSPROP.mod g
            RUN.mod g
            PARALLEL.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            XSI_ARRAY.mod g
            function.inc                                                
conv_rop_g.$(OBJ_EXT) : conv_rop_g.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            RUN.mod g
            COMPAR.mod g
            PARALLEL.mod g
            MATRIX.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PGCOR.mod g
            XSI_ARRAY.mod g
            function.inc                                                
conv_rop_s.$(OBJ_EXT) : conv_rop_s.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            RUN.mod g
            COMPAR.mod g
            PARALLEL.mod g
            MATRIX.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PGCOR.mod g
            PSCOR.mod g
            XSI_ARRAY.mod g
            function.inc                                                
conv_source_epp.$(OBJ_EXT) : conv_source_epp.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            RUN.mod g
            GEOMETRY.mod g
            COMPAR.mod g
            SENDRECV.mod g
            XSI_ARRAY.mod g
            PARALLEL.mod g
            MATRIX.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            RXNS.mod g
            INDICES.mod g
            PGCOR.mod g
            PSCOR.mod g
            VSHEAR.mod g
            ep_s1.inc                                                    g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            ep_s2.inc                                                   
copy_a.$(OBJ_EXT) : copy_a.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            PHYSPROP.mod g
            function.inc                                                
corner.$(OBJ_EXT) : corner.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            INDICES.mod g
            MATRIX.mod g
            CORNER.mod g
            FUNITS.mod g
            COMPAR.mod g
            function.inc                                                
correct_0.$(OBJ_EXT) : correct_0.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            PGCOR.mod g
            UR_FACS.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            COMPAR.mod g
            function.inc                                                
correct_1.$(OBJ_EXT) : correct_1.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            INDICES.mod g
            GEOMETRY.mod g
            PSCOR.mod g
            UR_FACS.mod g
            CONSTANT.mod g
            COMPAR.mod g
            SENDRECV.mod g
            ep_s1.inc                                                    g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            ep_s2.inc                                                   
dgtsl.$(OBJ_EXT) : dgtsl.f 
dif_u_is.$(OBJ_EXT) : dif_u_is.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            OUTPUT.mod g
            INDICES.mod g
            IS.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
dif_v_is.$(OBJ_EXT) : dif_v_is.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            OUTPUT.mod g
            INDICES.mod g
            IS.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
dif_w_is.$(OBJ_EXT) : dif_w_is.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            OUTPUT.mod g
            INDICES.mod g
            IS.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
discretize.$(OBJ_EXT) : discretize.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod 
display_resid.$(OBJ_EXT) : display_resid.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            RESIDUAL.mod g
            FLDVAR.mod g
            COMPAR.mod g
            GEOMETRY.mod 
drag_gs.$(OBJ_EXT) : drag_gs.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            RUN.mod g
            CONSTANT.mod g
            COMPAR.mod g
            DRAG.mod g
            SENDRECV.mod g
            DISCRETELEMENT.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
drag_ss.$(OBJ_EXT) : drag_ss.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            CONSTANT.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            COMPAR.mod g
            SENDRECV.mod g
            DRAG.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                
eosg.$(OBJ_EXT) : eosg.f g
            PARAM.mod g
            PARAM1.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            SCALES.mod g
            sc_p_g1.inc                                                  g
            sc_p_g2.inc                                                 
equal.$(OBJ_EXT) : equal.f g
            PARAM.mod g
            PARAM1.mod g
            INDICES.mod g
            PHYSPROP.mod 
error_routine.$(OBJ_EXT) : error_routine.f g
            FUNITS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod 
exchange.$(OBJ_EXT) : exchange.f g
            PARAM.mod g
            PARAM1.mod g
            COMPAR.mod 
exit.$(OBJ_EXT) : exit.f g
            FUNITS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod 
flow_to_vel.$(OBJ_EXT) : flow_to_vel.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            RUN.mod g
            BC.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod 
g_0.$(OBJ_EXT) : g_0.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg1.inc                                                 g
            fun_avg2.inc                                                
get_bc_area.$(OBJ_EXT) : get_bc_area.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            BC.mod g
            COMPAR.mod 
get_data.$(OBJ_EXT) : get_data.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            FUNITS.mod g
            COMPAR.mod g
            GRIDMAP.mod 
get_eq.$(OBJ_EXT) : get_eq.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            INDICES.mod 
get_flow_bc.$(OBJ_EXT) : get_flow_bc.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            BC.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            SENDRECV.mod g
            function.inc                                                
get_hloss.$(OBJ_EXT) : get_hloss.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            BC.mod g
            INDICES.mod g
            ENERGY.mod 
get_is.$(OBJ_EXT) : get_is.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            IS.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod 
get_philoss.$(OBJ_EXT) : get_philoss.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            BC.mod g
            INDICES.mod g
            ENERGY.mod g
            COMPAR.mod g
            function.inc                                                
get_smass.$(OBJ_EXT) : get_smass.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            INDICES.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            function.inc                                                
get_stats.$(OBJ_EXT) : get_stats.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            INDICES.mod g
            FUNITS.mod g
            RESIDUAL.mod g
            RUN.mod g
            COMPAR.mod g
            function.inc                                                
get_walls_bc.$(OBJ_EXT) : get_walls_bc.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            BC.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            SENDRECV.mod g
            function.inc                                                
in_bin_512.$(OBJ_EXT) : in_bin_512.f g
            MACHINE.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
in_bin_512i.$(OBJ_EXT) : in_bin_512i.f g
            MACHINE.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
init_ab_m.$(OBJ_EXT) : init_ab_m.f g
            PARAM.mod g
            PARAM1.mod g
            MATRIX.mod g
            PARALLEL.mod g
            COMPAR.mod 
init_fvars.$(OBJ_EXT) : init_fvars.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            INDICES.mod g
            SCALARS.mod g
            RXNS.mod g
            RUN.mod g
            COMPAR.mod 
init_namelist.$(OBJ_EXT) : init_namelist.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            OUTPUT.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            IC.mod g
            BC.mod g
            FLDVAR.mod g
            CONSTANT.mod g
            INDICES.mod g
            IS.mod g
            TOLERANC.mod g
            SCALES.mod g
            UR_FACS.mod g
            LEQSOL.mod g
            RESIDUAL.mod g
            RXNS.mod g
            SCALARS.mod g
            COMPAR.mod g
            PARALLEL.mod g
            namelist.inc                                                
init_resid.$(OBJ_EXT) : init_resid.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            RESIDUAL.mod 
iterate.$(OBJ_EXT) : iterate.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            INDICES.mod g
            FUNITS.mod g
            TIME_CPU.mod g
            PSCOR.mod g
            COEFF.mod g
            LEQSOL.mod g
            VISC_G.mod g
            PGCOR.mod g
            CONT.mod g
            SCALARS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            DISCRETELEMENT.mod g
            BC.mod g
            CONSTANT.mod 
k_epsilon_prop.$(OBJ_EXT) : k_epsilon_prop.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            DRAG.mod g
            RUN.mod g
            OUTPUT.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            VISC_G.mod g
            VISC_S.mod g
            TRACE.mod g
            INDICES.mod g
            CONSTANT.mod g
            VSHEAR.mod g
            TURB.mod g
            TOLERANC.mod g
            COMPAR.mod g
            TAU_G.mod g
            SENDRECV.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg2.inc                                                
leq_bicgs.$(OBJ_EXT) : leq_bicgs.f g
            PARAM.mod g
            PARAM1.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            COMPAR.mod g
            INDICES.mod g
            LEQSOL.mod g
            FUNITS.mod g
            PARALLEL.mod g
            MPI_UTILITY.mod g
            SENDRECV.mod g
            function.inc                                                
leq_gmres.$(OBJ_EXT) : leq_gmres.f g
            PARAM.mod g
            PARAM1.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            DEBUG.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            PARALLEL.mod g
            FUNITS.mod g
            GRIDMAP.mod g
            function.inc                                                
leq_sor.$(OBJ_EXT) : leq_sor.f g
            PARAM.mod g
            PARAM1.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
line_too_big.$(OBJ_EXT) : line_too_big.f 
location_check.$(OBJ_EXT) : location_check.f g
            PARAM.mod g
            PARAM1.mod g
            FUNITS.mod g
            GEOMETRY.mod 
location.$(OBJ_EXT) : location.f g
            PARAM.mod g
            PARAM1.mod 
machine.$(OBJ_EXT) : machine.f g
            MACHINE.mod g
            PARAM.mod g
            RUN.mod g
            FUNITS.mod 
make_upper_case.$(OBJ_EXT) : make_upper_case.f 
mark_phase_4_cor.$(OBJ_EXT) : mark_phase_4_cor.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            INDICES.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            CONSTANT.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS2) mark_phase_4_cor.f 
mfix.$(OBJ_EXT) : mfix.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            TIME_CPU.mod g
            FUNITS.mod g
            OUTPUT.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            PARALLEL_MPI.mod g
            function.inc                                                
mod_bc_i.$(OBJ_EXT) : mod_bc_i.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            function.inc                                                
mod_bc_j.$(OBJ_EXT) : mod_bc_j.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            function.inc                                                
mod_bc_k.$(OBJ_EXT) : mod_bc_k.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            function.inc                                                
open_file.$(OBJ_EXT) : open_file.f 
open_files.$(OBJ_EXT) : open_files.f g
            MACHINE.mod g
            FUNITS.mod g
            COMPAR.mod 
out_array_c.$(OBJ_EXT) : out_array_c.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            function.inc                                                
out_array.$(OBJ_EXT) : out_array.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            function.inc                                                
out_array_kc.$(OBJ_EXT) : out_array_kc.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            function.inc                                                
out_array_k.$(OBJ_EXT) : out_array_k.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            INDICES.mod g
            FUNITS.mod g
            COMPAR.mod g
            function.inc                                                
out_bin_512.$(OBJ_EXT) : out_bin_512.f g
            MACHINE.mod 
out_bin_512i.$(OBJ_EXT) : out_bin_512i.f g
            MACHINE.mod 
out_bin_512r.$(OBJ_EXT) : out_bin_512r.f g
            MACHINE.mod 
out_bin_r.$(OBJ_EXT) : out_bin_r.f g
            PARAM.mod 
parse_line.$(OBJ_EXT) : parse_line.f g
            PARAM.mod g
            PARAM1.mod g
            PARSE.mod g
            COMPAR.mod 
parse_resid_string.$(OBJ_EXT) : parse_resid_string.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            RESIDUAL.mod g
            FUNITS.mod g
            COMPAR.mod 
parse_rxn.$(OBJ_EXT) : parse_rxn.f g
            PARAM.mod g
            PARAM1.mod g
            PARSE.mod g
            RXNS.mod g
            COMPAR.mod 
partial_elim.$(OBJ_EXT) : partial_elim.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            GEOMETRY.mod g
            MATRIX.mod g
            PHYSPROP.mod g
            INDICES.mod g
            COMPAR.mod g
            DRAG.mod g
            FLDVAR.mod g
            RUN.mod g
            function.inc                                                 g
            fun_avg1.inc                                                 g
            fun_avg2.inc                                                
physical_prop.$(OBJ_EXT) : physical_prop.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            INDICES.mod g
            RUN.mod g
            TOLERANC.mod g
            CONSTANT.mod g
            COMPAR.mod g
            FUNITS.mod g
            cp_fun1.inc                                                  g
            function.inc                                                 g
            cp_fun2.inc                                                 
read_namelist.$(OBJ_EXT) : read_namelist.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            OUTPUT.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            IC.mod g
            IS.mod g
            BC.mod g
            FLDVAR.mod g
            CONSTANT.mod g
            INDICES.mod g
            TOLERANC.mod g
            FUNITS.mod g
            SCALES.mod g
            UR_FACS.mod g
            LEQSOL.mod g
            RESIDUAL.mod g
            RXNS.mod g
            SCALARS.mod g
            COMPAR.mod g
            PARALLEL.mod g
            DISCRETELEMENT.mod g
            usrnlst.inc                                                  g
            namelist.inc                                                 g
            des/desnamelist.inc                                         
read_res0.$(OBJ_EXT) : read_res0.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            RUN.mod g
            IC.mod g
            BC.mod g
            IS.mod g
            CONSTANT.mod g
            FUNITS.mod g
            OUTPUT.mod g
            SCALES.mod g
            UR_FACS.mod g
            TOLERANC.mod g
            LEQSOL.mod g
            SCALARS.mod g
            RXNS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            FLDVAR.mod 
read_res1.$(OBJ_EXT) : read_res1.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            RUN.mod g
            RXNS.mod g
            SCALARS.mod g
            FUNITS.mod g
            ENERGY.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            SENDRECV.mod 
remove_comment.$(OBJ_EXT) : remove_comment.f 
reset_new.$(OBJ_EXT) : reset_new.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            TRACE.mod g
            RUN.mod g
            SCALARS.mod 
rrates0.$(OBJ_EXT) : rrates0.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            RXNS.mod g
            ENERGY.mod g
            GEOMETRY.mod g
            RUN.mod g
            INDICES.mod g
            PHYSPROP.mod g
            CONSTANT.mod g
            FUNITS.mod g
            COMPAR.mod g
            SENDRECV.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
rrates.$(OBJ_EXT) : rrates.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            RXNS.mod g
            ENERGY.mod g
            GEOMETRY.mod g
            RUN.mod g
            INDICES.mod g
            PHYSPROP.mod g
            CONSTANT.mod g
            FUNITS.mod g
            COMPAR.mod g
            SENDRECV.mod g
            function.inc                                                
rrates_init.$(OBJ_EXT) : rrates_init.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            RXNS.mod g
            ENERGY.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
scalar_prop.$(OBJ_EXT) : scalar_prop.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            INDICES.mod g
            RUN.mod g
            SCALARS.mod g
            TOLERANC.mod g
            COMPAR.mod g
            SENDRECV.mod g
            function.inc                                                
seek_comment.$(OBJ_EXT) : seek_comment.f 
seek_end.$(OBJ_EXT) : seek_end.f 
set_bc0.$(OBJ_EXT) : set_bc0.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            PHYSPROP.mod g
            BC.mod g
            FLDVAR.mod g
            INDICES.mod g
            RUN.mod g
            FUNITS.mod g
            SCALES.mod g
            SCALARS.mod g
            BOUNDFUNIJK.mod g
            TOLERANC.mod g
            sc_p_g1.inc                                                  g
            function.inc                                                 g
            sc_p_g2.inc                                                 
set_bc1.$(OBJ_EXT) : set_bc1.f g
            PARAM.mod g
            PARAM1.mod g
            BC.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            RUN.mod g
            FUNITS.mod g
            COMPAR.mod g
            function.inc                                                
set_constants.$(OBJ_EXT) : set_constants.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            VISC_S.mod g
            ENERGY.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            CONSTANT.mod g
            RUN.mod g
            FUNITS.mod g
            DRAG.mod g
            COMPAR.mod 
set_constprop.$(OBJ_EXT) : set_constprop.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            VISC_S.mod g
            VISC_G.mod g
            ENERGY.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            CONSTANT.mod g
            RUN.mod g
            FUNITS.mod g
            DRAG.mod g
            COMPAR.mod g
            function.inc                                                
set_flags.$(OBJ_EXT) : set_flags.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            BC.mod g
            IS.mod g
            INDICES.mod g
            PHYSPROP.mod g
            FUNITS.mod g
            COMPAR.mod g
            SENDRECV.mod g
            SENDRECV3.mod g
            BOUNDFUNIJK.mod g
            MPI_UTILITY.mod g
            function.inc                                                 g
            function3.inc                                               
set_fluidbed_p.$(OBJ_EXT) : set_fluidbed_p.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            BC.mod g
            IC.mod g
            FLDVAR.mod g
            CONSTANT.mod g
            INDICES.mod g
            FUNITS.mod g
            SCALES.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            SENDRECV.mod g
            sc_p_g1.inc                                                  g
            b_force1.inc                                                 g
            function.inc                                                 g
            b_force2.inc                                                 g
            sc_p_g2.inc                                                 
set_geometry1.$(OBJ_EXT) : set_geometry1.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            RUN.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
set_geometry.$(OBJ_EXT) : set_geometry.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            GEOMETRY.mod g
            COMPAR.mod 
set_ic.$(OBJ_EXT) : set_ic.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            IC.mod g
            FLDVAR.mod g
            VISC_G.mod g
            INDICES.mod g
            SCALES.mod g
            ENERGY.mod g
            SCALARS.mod g
            COMPAR.mod g
            RUN.mod g
            SENDRECV.mod g
            sc_p_g1.inc                                                  g
            s_pr1.inc                                                    g
            function.inc                                                 g
            s_pr2.inc                                                    g
            sc_p_g2.inc                                                 
set_increments3.$(OBJ_EXT) : set_increments3.f g
            PARAM.mod g
            PARAM1.mod g
            INDICES.mod g
            GEOMETRY.mod g
            COMPAR.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            FUNITS.mod g
            function.inc                                                 g
            function3.inc                                               
set_increments.$(OBJ_EXT) : set_increments.f g
            PARAM.mod g
            PARAM1.mod g
            INDICES.mod g
            GEOMETRY.mod g
            COMPAR.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            FUNITS.mod g
            function.inc                                                
set_index1a3.$(OBJ_EXT) : set_index1a3.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            COMPAR.mod g
            FLDVAR.mod g
            INDICES.mod g
            BOUNDFUNIJK3.mod g
            function.inc                                                
set_index1a.$(OBJ_EXT) : set_index1a.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            COMPAR.mod g
            FLDVAR.mod g
            INDICES.mod g
            BOUNDFUNIJK.mod g
            function.inc                                                
set_index1.$(OBJ_EXT) : set_index1.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            CONSTANT.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
set_l_scale.$(OBJ_EXT) : set_l_scale.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            CONSTANT.mod g
            VISC_G.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod 
set_max2.$(OBJ_EXT) : set_max2.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            COMPAR.mod 
set_mw_mix_g.$(OBJ_EXT) : set_mw_mix_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            CONSTANT.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
set_outflow.$(OBJ_EXT) : set_outflow.f g
            PARAM.mod g
            PARAM1.mod g
            BC.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            SCALARS.mod g
            RUN.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
set_ro_g.$(OBJ_EXT) : set_ro_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            CONSTANT.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
set_wall_bc.$(OBJ_EXT) : set_wall_bc.f g
            PARAM.mod g
            PARAM1.mod g
            BC.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            RUN.mod g
            FUNITS.mod g
            COMPAR.mod g
            function.inc                                                
shift_dxyz.$(OBJ_EXT) : shift_dxyz.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod 
solve_continuity.$(OBJ_EXT) : solve_continuity.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            INDICES.mod g
            RESIDUAL.mod g
            CONT.mod g
            LEQSOL.mod g
            AMBM.mod 
solve_energy_eq.$(OBJ_EXT) : solve_energy_eq.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            INDICES.mod g
            DRAG.mod g
            RESIDUAL.mod g
            UR_FACS.mod g
            PGCOR.mod g
            PSCOR.mod g
            LEQSOL.mod g
            BC.mod g
            ENERGY.mod g
            RXNS.mod g
            AMBM.mod g
            TMP_ARRAY.mod g
            TMP_ARRAY1.mod g
            COMPAR.mod g
            DISCRETELEMENT.mod g
            MFLUX.mod g
            radtn1.inc                                                   g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                    g
            radtn2.inc                                                  
solve_epp.$(OBJ_EXT) : solve_epp.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            PSCOR.mod g
            RESIDUAL.mod g
            LEQSOL.mod g
            PHYSPROP.mod g
            AMBM.mod 
solve_granular_energy.$(OBJ_EXT) : solve_granular_energy.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            CONSTANT.mod g
            OUTPUT.mod g
            INDICES.mod g
            DRAG.mod g
            RESIDUAL.mod g
            UR_FACS.mod g
            PGCOR.mod g
            PSCOR.mod g
            LEQSOL.mod g
            BC.mod g
            ENERGY.mod g
            RXNS.mod g
            AMBM.mod g
            TMP_ARRAY.mod g
            COMPAR.mod g
            MFLUX.mod g
            radtn1.inc                                                   g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                    g
            radtn2.inc                                                  
solve_k_epsilon_eq.$(OBJ_EXT) : solve_k_epsilon_eq.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            INDICES.mod g
            DRAG.mod g
            RESIDUAL.mod g
            UR_FACS.mod g
            PGCOR.mod g
            PSCOR.mod g
            LEQSOL.mod g
            BC.mod g
            ENERGY.mod g
            RXNS.mod g
            TURB.mod g
            USR.mod g
            AMBM.mod g
            TMP_ARRAY.mod g
            COMPAR.mod g
            MFLUX.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg1.inc                                                 g
            fun_avg2.inc                                                
solve_lin_eq.$(OBJ_EXT) : solve_lin_eq.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            COMPAR.mod 
solve_pp_g.$(OBJ_EXT) : solve_pp_g.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            PGCOR.mod g
            RESIDUAL.mod g
            LEQSOL.mod g
            RUN.mod g
            AMBM.mod 
solve_scalar_eq.$(OBJ_EXT) : solve_scalar_eq.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            INDICES.mod g
            DRAG.mod g
            RESIDUAL.mod g
            UR_FACS.mod g
            PGCOR.mod g
            PSCOR.mod g
            LEQSOL.mod g
            BC.mod g
            ENERGY.mod g
            RXNS.mod g
            SCALARS.mod g
            AMBM.mod g
            TMP_ARRAY.mod g
            COMPAR.mod g
            MFLUX.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
solve_species_eq.$(OBJ_EXT) : solve_species_eq.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            INDICES.mod g
            DRAG.mod g
            RESIDUAL.mod g
            UR_FACS.mod g
            PGCOR.mod g
            PSCOR.mod g
            LEQSOL.mod g
            BC.mod g
            ENERGY.mod g
            RXNS.mod g
            AMBM.mod g
            MATRIX.mod g
            CHISCHEME.mod g
            TMP_ARRAY.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            SENDRECV.mod g
            MFLUX.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
solve_vel_star.$(OBJ_EXT) : solve_vel_star.f g
            PARAM.mod g
            PARAM1.mod g
            TOLERANC.mod g
            RUN.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            OUTPUT.mod g
            INDICES.mod g
            DRAG.mod g
            RESIDUAL.mod g
            UR_FACS.mod g
            PGCOR.mod g
            PSCOR.mod g
            LEQSOL.mod g
            AMBM.mod g
            TMP_ARRAY1.mod g
            TMP_ARRAY.mod g
            COMPAR.mod g
            DISCRETELEMENT.mod 
source_granular_energy.$(OBJ_EXT) : source_granular_energy.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            PHYSPROP.mod g
            RUN.mod g
            DRAG.mod g
            GEOMETRY.mod g
            FLDVAR.mod g
            VISC_G.mod g
            VISC_S.mod g
            TRACE.mod g
            TURB.mod g
            INDICES.mod g
            CONSTANT.mod g
            TOLERANC.mod g
            COMPAR.mod g
            s_pr1.inc                                                    g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            ep_s2.inc                                                    g
            fun_avg2.inc                                                 g
            s_pr2.inc                                                   
source_phi.$(OBJ_EXT) : source_phi.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            TAU_S.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
source_pp_g.$(OBJ_EXT) : source_pp_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            RXNS.mod g
            RUN.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PGCOR.mod g
            BC.mod g
            VSHEAR.mod g
            XSI_ARRAY.mod g
            COMPAR.mod g
            UR_FACS.mod g
            function.inc                                                
source_rop_g.$(OBJ_EXT) : source_rop_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            FLDVAR.mod g
            RXNS.mod g
            RUN.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PGCOR.mod g
            COMPAR.mod g
            function.inc                                                
source_rop_s.$(OBJ_EXT) : source_rop_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            FLDVAR.mod g
            RXNS.mod g
            RUN.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PGCOR.mod g
            PSCOR.mod g
            COMPAR.mod g
            function.inc                                                
source_u_g.$(OBJ_EXT) : source_u_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_G.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            TAU_G.mod g
            BC.mod g
            COMPAR.mod g
            SENDRECV.mod g
            OUTPUT.mod g
            TURB.mod g
            MPI_UTILITY.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
source_u_s.$(OBJ_EXT) : source_u_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            TAU_S.mod g
            BC.mod g
            COMPAR.mod g
            SENDRECV.mod g
            OUTPUT.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
source_v_g.$(OBJ_EXT) : source_v_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_G.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            TAU_G.mod g
            BC.mod g
            VSHEAR.mod g
            COMPAR.mod g
            SENDRECV.mod g
            OUTPUT.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
source_v_s.$(OBJ_EXT) : source_v_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            TAU_S.mod g
            BC.mod g
            VSHEAR.mod g
            COMPAR.mod g
            SENDRECV.mod g
            OUTPUT.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
source_w_g.$(OBJ_EXT) : source_w_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_G.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            TAU_G.mod g
            BC.mod g
            COMPAR.mod g
            SENDRECV.mod g
            OUTPUT.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
source_w_s.$(OBJ_EXT) : source_w_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            TAU_S.mod g
            BC.mod g
            COMPAR.mod g
            SENDRECV.mod g
            OUTPUT.mod g
            b_force1.inc                                                 g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                    g
            b_force2.inc                                                
tau_u_g.$(OBJ_EXT) : tau_u_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_G.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            COMPAR.mod g
            SENDRECV.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
tau_u_s.$(OBJ_EXT) : tau_u_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            VSHEAR.mod g
            SENDRECV.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
tau_v_g.$(OBJ_EXT) : tau_v_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_G.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            SENDRECV.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
tau_v_s.$(OBJ_EXT) : tau_v_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            SENDRECV.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
tau_w_g.$(OBJ_EXT) : tau_w_g.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_G.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            SENDRECV.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
tau_w_s.$(OBJ_EXT) : tau_w_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_S.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            SENDRECV.mod g
            COMPAR.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
test_lin_eq.$(OBJ_EXT) : test_lin_eq.f g
            PARAM.mod g
            PARAM1.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            function.inc                                                
time_march.$(OBJ_EXT) : time_march.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            OUTPUT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            PGCOR.mod g
            PSCOR.mod g
            CONT.mod g
            COEFF.mod g
            TAU_G.mod g
            TAU_S.mod g
            VISC_G.mod g
            VISC_S.mod g
            FUNITS.mod g
            VSHEAR.mod g
            SCALARS.mod g
            DRAG.mod g
            RXNS.mod g
            COMPAR.mod g
            TIME_CPU.mod g
            DISCRETELEMENT.mod 
transfer.$(OBJ_EXT) : transfer.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            INDICES.mod 
transport_prop.$(OBJ_EXT) : transport_prop.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            INDICES.mod g
            RUN.mod g
            TOLERANC.mod g
            COMPAR.mod 
undef_2_0.$(OBJ_EXT) : undef_2_0.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            COMPAR.mod 
under_relax.$(OBJ_EXT) : under_relax.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            INDICES.mod g
            COMPAR.mod g
            SENDRECV.mod g
            function.inc                                                
update_old.$(OBJ_EXT) : update_old.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            RUN.mod g
            TRACE.mod g
            VISC_S.mod g
            SCALARS.mod 
usr0.$(OBJ_EXT) : usr0.f g
            USR.mod 
usr1.$(OBJ_EXT) : usr1.f g
            USR.mod 
usr2.$(OBJ_EXT) : usr2.f g
            USR.mod 
usr3.$(OBJ_EXT) : usr3.f g
            USR.mod 
usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f g
            usrnlst.inc                                                 
usr_write_out0.$(OBJ_EXT) : usr_write_out0.f 
usr_write_out1.$(OBJ_EXT) : usr_write_out1.f 
vavg_u_g.$(OBJ_EXT) : vavg_u_g.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            BC.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            INDICES.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            function.inc                                                
vavg_u_s.$(OBJ_EXT) : vavg_u_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            BC.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            INDICES.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
vavg_v_g.$(OBJ_EXT) : vavg_v_g.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            BC.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            INDICES.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            function.inc                                                
vavg_v_s.$(OBJ_EXT) : vavg_v_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            BC.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            INDICES.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
vavg_w_g.$(OBJ_EXT) : vavg_w_g.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            BC.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            INDICES.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            function.inc                                                
vavg_w_s.$(OBJ_EXT) : vavg_w_s.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            BC.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            INDICES.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            ep_s1.inc                                                    g
            function.inc                                                 g
            ep_s2.inc                                                   
vf_gs_x.$(OBJ_EXT) : vf_gs_x.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            COMPAR.mod g
            DRAG.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                
vf_gs_y.$(OBJ_EXT) : vf_gs_y.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            COMPAR.mod g
            DRAG.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                
vf_gs_z.$(OBJ_EXT) : vf_gs_z.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            COMPAR.mod g
            DRAG.mod g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                
write_ab_m.$(OBJ_EXT) : write_ab_m.f g
            PARAM.mod g
            PARAM1.mod g
            MATRIX.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            INDICES.mod g
            function.inc                                                
write_ab_m_var.$(OBJ_EXT) : write_ab_m_var.f g
            PARAM.mod g
            PARAM1.mod g
            MATRIX.mod g
            GEOMETRY.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            INDICES.mod g
            function.inc                                                
write_error.$(OBJ_EXT) : write_error.f g
            PARAM.mod g
            PARAM1.mod g
            FUNITS.mod 
write_header.$(OBJ_EXT) : write_header.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            OUTPUT.mod g
            FUNITS.mod g
            COMPAR.mod 
write_out0.$(OBJ_EXT) : write_out0.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            OUTPUT.mod g
            PHYSPROP.mod g
            GEOMETRY.mod g
            IC.mod g
            BC.mod g
            IS.mod g
            FLDVAR.mod g
            CONSTANT.mod g
            INDICES.mod g
            FUNITS.mod g
            TOLERANC.mod g
            SCALES.mod g
            SCALARS.mod g
            UR_FACS.mod g
            LEQSOL.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            SENDRECV.mod g
            function.inc                                                
write_out1.$(OBJ_EXT) : write_out1.f g
            PARAM.mod g
            PARAM1.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            RUN.mod g
            SCALARS.mod g
            FUNITS.mod g
            RXNS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod 
write_out3.$(OBJ_EXT) : write_out3.f g
            FUNITS.mod g
            COMPAR.mod 
write_res0.$(OBJ_EXT) : write_res0.f g
            PARAM.mod g
            PARAM1.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            RUN.mod g
            IC.mod g
            IS.mod g
            BC.mod g
            CONSTANT.mod g
            FUNITS.mod g
            OUTPUT.mod g
            SCALES.mod g
            SCALARS.mod g
            RXNS.mod g
            UR_FACS.mod g
            LEQSOL.mod g
            TOLERANC.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            SENDRECV.mod 
write_res1.$(OBJ_EXT) : write_res1.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            RUN.mod g
            SCALARS.mod g
            RXNS.mod g
            FUNITS.mod g
            OUTPUT.mod g
            ENERGY.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            SENDRECV.mod 
write_spx0.$(OBJ_EXT) : write_spx0.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            FUNITS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod 
write_spx1.$(OBJ_EXT) : write_spx1.f g
            PARAM.mod g
            PARAM1.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            RUN.mod g
            FUNITS.mod g
            SCALARS.mod g
            OUTPUT.mod g
            RXNS.mod g
            COMPAR.mod g
            MPI_UTILITY.mod g
            SENDRECV.mod 
write_table.$(OBJ_EXT) : write_table.f g
            PARAM.mod g
            PARAM1.mod g
            FUNITS.mod 
write_usr0.$(OBJ_EXT) : write_usr0.f 
write_usr1.$(OBJ_EXT) : write_usr1.f 
xerbla.$(OBJ_EXT) : xerbla.f g
            COMPAR.mod 
zero_array.$(OBJ_EXT) : zero_array.f g
            PARAM.mod g
            PARAM1.mod 
zero_norm_vel.$(OBJ_EXT) : zero_norm_vel.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            GEOMETRY.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            INDICES.mod g
            IS.mod g
            COMPAR.mod g
            function.inc                                                
add_part_to_link_list.$(OBJ_EXT) : ./cohesion/add_part_to_link_list.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/add_part_to_link_list.f 
calc_app_coh_force.$(OBJ_EXT) : ./cohesion/calc_app_coh_force.f g
            DISCRETELEMENT.mod g
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_app_coh_force.f 
calc_cap_coh_force.$(OBJ_EXT) : ./cohesion/calc_cap_coh_force.f g
            DISCRETELEMENT.mod g
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_cap_coh_force.f 
calc_cohesive_forces.$(OBJ_EXT) : ./cohesion/calc_cohesive_forces.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_cohesive_forces.f 
calc_esc_coh_force.$(OBJ_EXT) : ./cohesion/calc_esc_coh_force.f g
            DISCRETELEMENT.mod g
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_esc_coh_force.f 
calc_square_well.$(OBJ_EXT) : ./cohesion/calc_square_well.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_square_well.f 
calc_van_der_waals.$(OBJ_EXT) : ./cohesion/calc_van_der_waals.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_van_der_waals.f 
check_link.$(OBJ_EXT) : ./cohesion/check_link.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_link.f 
check_sw_wall_interaction.$(OBJ_EXT) : ./cohesion/check_sw_wall_interaction.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_sw_wall_interaction.f 
check_vdw_wall_interaction.$(OBJ_EXT) : ./cohesion/check_vdw_wall_interaction.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_vdw_wall_interaction.f 
initialize_cohesion_parameters.$(OBJ_EXT) : ./cohesion/initialize_cohesion_parameters.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/initialize_cohesion_parameters.f 
initialize_coh_int_search.$(OBJ_EXT) : ./cohesion/initialize_coh_int_search.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/initialize_coh_int_search.f 
linked_interaction_eval.$(OBJ_EXT) : ./cohesion/linked_interaction_eval.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/linked_interaction_eval.f 
remove_part_from_link_list.$(OBJ_EXT) : ./cohesion/remove_part_from_link_list.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/remove_part_from_link_list.f 
unlinked_interaction_eval.$(OBJ_EXT) : ./cohesion/unlinked_interaction_eval.f g
            DISCRETELEMENT.mod g
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/unlinked_interaction_eval.f 
update_search_grids.$(OBJ_EXT) : ./cohesion/update_search_grids.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/update_search_grids.f 
calc_force_des.$(OBJ_EXT) : ./des/calc_force_des.f g
            DISCRETELEMENT.mod g
            GEOMETRY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_force_des.f 
cfassign.$(OBJ_EXT) : ./des/cfassign.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfassign.f 
cffctow.$(OBJ_EXT) : ./des/cffctow.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffctow.f 
cffn.$(OBJ_EXT) : ./des/cffn.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffn.f 
cffnwall.$(OBJ_EXT) : ./des/cffnwall.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffnwall.f 
cfft.$(OBJ_EXT) : ./des/cfft.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfft.f 
cfftwall.$(OBJ_EXT) : ./des/cfftwall.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfftwall.f 
cfincrementaloverlaps.$(OBJ_EXT) : ./des/cfincrementaloverlaps.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfincrementaloverlaps.f 
cfnewvalues.$(OBJ_EXT) : ./des/cfnewvalues.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            COMPAR.mod g
            SENDRECV.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_G.mod g
            RXNS.mod g
            RUN.mod g
            GEOMETRY.mod g
            INDICES.mod g
            DRAG.mod g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfnewvalues.f 
cfnocontact.$(OBJ_EXT) : ./des/cfnocontact.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfnocontact.f 
cfnormal.$(OBJ_EXT) : ./des/cfnormal.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfnormal.f 
cfoutofbox.$(OBJ_EXT) : ./des/cfoutofbox.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfoutofbox.f 
cfperiodicwallneighbourx.$(OBJ_EXT) : ./des/cfperiodicwallneighbourx.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwallneighbourx.f 
cfperiodicwallneighboury.$(OBJ_EXT) : ./des/cfperiodicwallneighboury.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwallneighboury.f 
cfperiodicwallneighbourz.$(OBJ_EXT) : ./des/cfperiodicwallneighbourz.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwallneighbourz.f 
cfperiodicwallx.$(OBJ_EXT) : ./des/cfperiodicwallx.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwallx.f 
cfperiodicwally.$(OBJ_EXT) : ./des/cfperiodicwally.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwally.f 
cfperiodicwallz.$(OBJ_EXT) : ./des/cfperiodicwallz.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfperiodicwallz.f 
cfrelvel.$(OBJ_EXT) : ./des/cfrelvel.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfrelvel.f 
cfslide.$(OBJ_EXT) : ./des/cfslide.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslide.f 
cfslidewall.$(OBJ_EXT) : ./des/cfslidewall.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslidewall.f 
cftangent.$(OBJ_EXT) : ./des/cftangent.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cftangent.f 
cftotaloverlaps.$(OBJ_EXT) : ./des/cftotaloverlaps.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cftotaloverlaps.f 
cfupdateold.$(OBJ_EXT) : ./des/cfupdateold.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfupdateold.f 
cfvrn.$(OBJ_EXT) : ./des/cfvrn.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfvrn.f 
cfvrt.$(OBJ_EXT) : ./des/cfvrt.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfvrt.f 
cfwallcontact.$(OBJ_EXT) : ./des/cfwallcontact.f g
            DISCRETELEMENT.mod g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            RUN.mod g
            GEOMETRY.mod g
            MATRIX.mod g
            INDICES.mod g
            PHYSPROP.mod g
            DRAG.mod g
            CONSTANT.mod g
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallcontact.f 
cfwallposvel.$(OBJ_EXT) : ./des/cfwallposvel.f g
            DISCRETELEMENT.mod g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            RUN.mod g
            GEOMETRY.mod g
            MATRIX.mod g
            INDICES.mod g
            PHYSPROP.mod g
            DRAG.mod g
            CONSTANT.mod g
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallposvel.f 
des_calc_d.$(OBJ_EXT) : ./des/des_calc_d.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            INDICES.mod g
            PHYSPROP.mod g
            RUN.mod g
            SCALES.mod g
            COMPAR.mod g
            SENDRECV.mod g
            DISCRETELEMENT.mod g
            ep_s1.inc                                                    g
            fun_avg1.inc                                                 g
            function.inc                                                 g
            fun_avg2.inc                                                 g
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_calc_d.f 
des_granular_temperature.$(OBJ_EXT) : ./des/des_granular_temperature.f g
            DISCRETELEMENT.mod g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            RUN.mod g
            GEOMETRY.mod g
            MATRIX.mod g
            INDICES.mod g
            PHYSPROP.mod g
            DRAG.mod g
            CONSTANT.mod g
            COMPAR.mod g
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_granular_temperature.f 
des_init_namelist.$(OBJ_EXT) : ./des/des_init_namelist.f g
            DISCRETELEMENT.mod g
            des/desnamelist.inc                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_namelist.f 
des_inlet_outlet.$(OBJ_EXT) : ./des/des_inlet_outlet.f g
            DISCRETELEMENT.mod g
            GEOMETRY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_inlet_outlet.f 
des_time_march.$(OBJ_EXT) : ./des/des_time_march.f g
            PARAM.mod g
            PARAM1.mod g
            RUN.mod g
            OUTPUT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            GEOMETRY.mod g
            PGCOR.mod g
            PSCOR.mod g
            CONT.mod g
            COEFF.mod g
            TAU_G.mod g
            TAU_S.mod g
            VISC_G.mod g
            VISC_S.mod g
            FUNITS.mod g
            VSHEAR.mod g
            SCALARS.mod g
            DRAG.mod g
            RXNS.mod g
            COMPAR.mod g
            TIME_CPU.mod g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_time_march.f 
drag_fgs.$(OBJ_EXT) : ./des/drag_fgs.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_G.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            TAU_G.mod g
            BC.mod g
            COMPAR.mod g
            SENDRECV.mod g
            DISCRETELEMENT.mod g
            DRAG.mod g
            function.inc                                                 g
            ep_s1.inc                                                    g
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/drag_fgs.f 
gas_drag.$(OBJ_EXT) : ./des/gas_drag.f g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            MATRIX.mod g
            SCALES.mod g
            CONSTANT.mod g
            PHYSPROP.mod g
            FLDVAR.mod g
            VISC_G.mod g
            RXNS.mod g
            RUN.mod g
            TOLERANC.mod g
            GEOMETRY.mod g
            INDICES.mod g
            IS.mod g
            TAU_G.mod g
            BC.mod g
            COMPAR.mod g
            SENDRECV.mod g
            DISCRETELEMENT.mod g
            DRAG.mod g
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/gas_drag.f 
make_arrays_des.$(OBJ_EXT) : ./des/make_arrays_des.f g
            FUNITS.mod g
            COMPAR.mod g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/make_arrays_des.f 
neighbour.$(OBJ_EXT) : ./des/neighbour.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/neighbour.f 
nsquare.$(OBJ_EXT) : ./des/nsquare.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/nsquare.f 
octree.$(OBJ_EXT) : ./des/octree.f g
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/octree.f 
particles_in_cell.$(OBJ_EXT) : ./des/particles_in_cell.f g
            DISCRETELEMENT.mod g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            RUN.mod g
            GEOMETRY.mod g
            MATRIX.mod g
            INDICES.mod g
            PHYSPROP.mod g
            DRAG.mod g
            CONSTANT.mod g
            COMPAR.mod g
            SENDRECV.mod g
            function.inc                                                 g
            ep_s1.inc                                                    g
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/particles_in_cell.f 
periodic_wall_calc_force_des.$(OBJ_EXT) : ./des/periodic_wall_calc_force_des.f g
            DISCRETELEMENT.mod g
            GEOMETRY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/periodic_wall_calc_force_des.f 
pressure_drop.$(OBJ_EXT) : ./des/pressure_drop.f g
            DISCRETELEMENT.mod g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            RUN.mod g
            GEOMETRY.mod g
            MATRIX.mod g
            INDICES.mod g
            PHYSPROP.mod g
            DRAG.mod g
            CONSTANT.mod g
            COMPAR.mod g
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/pressure_drop.f 
print_vel.$(OBJ_EXT) : ./des/print_vel.f g
            DISCRETELEMENT.mod g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            RUN.mod g
            GEOMETRY.mod g
            MATRIX.mod g
            INDICES.mod g
            PHYSPROP.mod g
            DRAG.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/print_vel.f 
quadtree.$(OBJ_EXT) : ./des/quadtree.f g
            DISCRETELEMENT.mod g
            PARAM.mod g
            PARAM1.mod g
            PARALLEL.mod g
            FLDVAR.mod g
            RUN.mod g
            GEOMETRY.mod g
            MATRIX.mod g
            INDICES.mod g
            PHYSPROP.mod g
            DRAG.mod g
            CONSTANT.mod g
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/quadtree.f 
