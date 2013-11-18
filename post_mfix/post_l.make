.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
post_mfix : \
    ambm.mod \
    bc.mod \
    cdist.mod \
    constant.mod \
    cont.mod \
    correl.mod \
    des_bc.mod \
    drag.mod \
    energy.mod \
    fldvar.mod \
    funits.mod \
    geometry.mod \
    ic.mod \
    indices.mod \
    is.mod \
    kintheory.mod \
    leqsol.mod \
    machine.mod \
    stiff_chem.mod \
    mflux.mod \
    mfix_netcdf.mod \
    output.mod \
    parallel.mod \
    param1.mod \
    param.mod \
    parse.mod \
    pgcor.mod \
    physprop.mod \
    post3d.mod \
    pscor.mod \
    ps.mod \
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
    progress_bar.mod \
    quadric.mod \
    stl.mod \
    vtk.mod \
    cutcell.mod \
    polygon.mod \
    dashboard.mod \
    compar.mod \
    mpi.mod \
    dbg_util.mod \
    parallel_mpi.mod \
    debug.mod \
    gridmap.mod \
    mpi_utility.mod \
    sendrecv.mod \
    boundfunijk.mod \
    discretelement.mod \
    ghdtheory.mod \
    post_precision.mod \
    paralleldata.mod \
    usr_input.mod \
    mfix_pic.mod \
    qmom_kinetic_equation.mod \
    qmomk_parameters.mod \
    des_ic.mod \
    des_thermo.mod \
    des_rxns.mod \
    rxn_com.mod \
    machine.$(OBJ_EXT) \
    allocate_arrays.$(OBJ_EXT) \
    any_more_data.$(OBJ_EXT) \
    calc_cell2.$(OBJ_EXT) \
    calc_corr_01.$(OBJ_EXT) \
    calc_corr_type_1.$(OBJ_EXT) \
    calc_distance.$(OBJ_EXT) \
    calc_ep_g.$(OBJ_EXT) \
    calc_mu_s.$(OBJ_EXT) \
    calc_mw.$(OBJ_EXT) \
    calc_p_star.$(OBJ_EXT) \
    calc_quantities.$(OBJ_EXT) \
    calc_ro_g.$(OBJ_EXT) \
    calc_vol.$(OBJ_EXT) \
    check_data_03.$(OBJ_EXT) \
    check_data_04.$(OBJ_EXT) \
    check_data_05.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    compare.$(OBJ_EXT) \
    deallocate_arrays.$(OBJ_EXT) \
    cartesian_grid_init_namelist.$(OBJ_EXT) \
    eosg.$(OBJ_EXT) \
    error_routine.$(OBJ_EXT) \
    examine_data.$(OBJ_EXT) \
    exit.$(OBJ_EXT) \
    f_init_data.$(OBJ_EXT) \
    file_handle.$(OBJ_EXT) \
    finit.$(OBJ_EXT) \
    flow_gx.$(OBJ_EXT) \
    flow_gy.$(OBJ_EXT) \
    flow_gz.$(OBJ_EXT) \
    flow_sx.$(OBJ_EXT) \
    flow_sy.$(OBJ_EXT) \
    flow_sz.$(OBJ_EXT) \
    g_0.$(OBJ_EXT) \
    gas_flux.$(OBJ_EXT) \
    get_file_name.$(OBJ_EXT) \
    get_file_status.$(OBJ_EXT) \
    get_index.$(OBJ_EXT) \
    get_location.$(OBJ_EXT) \
    get_mu_s.$(OBJ_EXT) \
    get_same_time.$(OBJ_EXT) \
    get_selection.$(OBJ_EXT) \
    get_substr.$(OBJ_EXT) \
    granular_qty.$(OBJ_EXT) \
    header_main.$(OBJ_EXT) \
    ik_avg.$(OBJ_EXT) \
    ik_avg_out.$(OBJ_EXT) \
    init_namelist.$(OBJ_EXT) \
    interp_res.$(OBJ_EXT) \
    line_too_big.$(OBJ_EXT) \
    kintheory_energy_dissipation_ss.$(OBJ_EXT) \
    main_f.$(OBJ_EXT) \
    make_upper_case.$(OBJ_EXT) \
    open_file.$(OBJ_EXT) \
    out_from_res.$(OBJ_EXT) \
    out_from_spx.$(OBJ_EXT) \
    out_spec_time.$(OBJ_EXT) \
    out_time.$(OBJ_EXT) \
    parse_line.$(OBJ_EXT) \
    parse_rxn.$(OBJ_EXT) \
    print_out.$(OBJ_EXT) \
    read_namelist.$(OBJ_EXT) \
    read_res0.$(OBJ_EXT) \
    read_res1.$(OBJ_EXT) \
    read_spx0.$(OBJ_EXT) \
    read_spx1.$(OBJ_EXT) \
    remove_comment.$(OBJ_EXT) \
    res_from_spx.$(OBJ_EXT) \
    seek_comment.$(OBJ_EXT) \
    seek_end.$(OBJ_EXT) \
    seek_time.$(OBJ_EXT) \
    select_spx_rec.$(OBJ_EXT) \
    set_constants.$(OBJ_EXT) \
    set_dollar.$(OBJ_EXT) \
    set_geometry.$(OBJ_EXT) \
    set_increments.$(OBJ_EXT) \
    set_index1.$(OBJ_EXT) \
    set_index1a.$(OBJ_EXT) \
    set_max2.$(OBJ_EXT) \
    set_read_spx.$(OBJ_EXT) \
    shift_dxyz.$(OBJ_EXT) \
    sol_flux.$(OBJ_EXT) \
    strcmp.$(OBJ_EXT) \
    streqs.$(OBJ_EXT) \
    time_avg.$(OBJ_EXT) \
    usr_init_namelist.$(OBJ_EXT) \
    usr_post.$(OBJ_EXT) \
    usr_write_out1.$(OBJ_EXT) \
    write_out1.$(OBJ_EXT) \
    write_res0.$(OBJ_EXT) \
    write_res1.$(OBJ_EXT) \
    write_spx0.$(OBJ_EXT) \
    write_spx1.$(OBJ_EXT) \
    usr_input.$(OBJ_EXT) \
    write_error.$(OBJ_EXT) \
    transport_coeff_ghd.$(OBJ_EXT) \
    ghd.$(OBJ_EXT) \
    cooling_rate.$(OBJ_EXT) \
    cooling_rate_tc.$(OBJ_EXT) \
    pressure.$(OBJ_EXT) \
    bulk_viscosity.$(OBJ_EXT) \
    shear_viscosity.$(OBJ_EXT) \
    thermal_diffusivity.$(OBJ_EXT) \
    mass_mobility.$(OBJ_EXT) \
    thermal_conductivity.$(OBJ_EXT) \
    thermal_mobility.$(OBJ_EXT) \
    ordinary_diff.$(OBJ_EXT) \
    dufour_coeff.$(OBJ_EXT) \
    chi_ij_GHD.$(OBJ_EXT) \
    des_init_namelist.$(OBJ_EXT) \
    ornl_header.$(OBJ_EXT) \
    ornl_util.$(OBJ_EXT) \
    ornl_stats.$(OBJ_EXT) \
    ornl_stats_c.$(OBJ_EXT) \
    ornl_corr.$(OBJ_EXT) \
    ornl_corr_c.$(OBJ_EXT) \
    ornl_pca.$(OBJ_EXT) \
    ornl_ft.$(OBJ_EXT) \
    ornl_ft_c.$(OBJ_EXT) \
    ornl_filt.$(OBJ_EXT) \
    ornl_filt_c.$(OBJ_EXT) \
    ornl_zone.$(OBJ_EXT) \
    ornl_sym.$(OBJ_EXT) \
    qmomk_init_namelist.$(OBJ_EXT) \
    calc_trd_s.$(OBJ_EXT) \
    get_delh.$(OBJ_EXT) \
    define_quadrics.$(OBJ_EXT) \
    check_data_cartesian.$(OBJ_EXT) \
    get_poly_data.$(OBJ_EXT) \
    get_stl_data.$(OBJ_EXT) \
    get_connectivity.$(OBJ_EXT) \
    cut_cell_preprocessing.$(OBJ_EXT) \
    eval_usr_fct.$(OBJ_EXT) \
    allocate_cut_cell_arrays.$(OBJ_EXT) \
    allocate_dummy_cut_cell_arrays.$(OBJ_EXT) \
    calc_vort_out.$(OBJ_EXT) \
    deallocate_cut_cell_arrays.$(OBJ_EXT) \
    get_alpha.$(OBJ_EXT) \
    get_cut_cell_flags.$(OBJ_EXT) \
    get_cut_cell_volume_area.$(OBJ_EXT) \
    get_master.$(OBJ_EXT) \
    set_Odxyz.$(OBJ_EXT) \
    vtk_out.$(OBJ_EXT) \
    write_progress_bar.$(OBJ_EXT) \
    dmp_cartesian.$(OBJ_EXT) \
    set_geometry1.$(OBJ_EXT) \
    read_database.$(OBJ_EXT) \
    readTherm.$(OBJ_EXT) \
    get_values.$(OBJ_EXT) \
    get_bc_area.$(OBJ_EXT) \
    flow_to_vel.$(OBJ_EXT) \
    
	$(LINK_CMD) $(LINK_FLAGS) \
    ambm_mod.$(OBJ_EXT) \
    bc_mod.$(OBJ_EXT) \
    cdist_mod.$(OBJ_EXT) \
    constant_mod.$(OBJ_EXT) \
    cont_mod.$(OBJ_EXT) \
    correl_mod.$(OBJ_EXT) \
    des_bc_mod.$(OBJ_EXT) \
    drag_mod.$(OBJ_EXT) \
    energy_mod.$(OBJ_EXT) \
    fldvar_mod.$(OBJ_EXT) \
    funits_mod.$(OBJ_EXT) \
    geometry_mod.$(OBJ_EXT) \
    ic_mod.$(OBJ_EXT) \
    indices_mod.$(OBJ_EXT) \
    is_mod.$(OBJ_EXT) \
    kintheory_mod.$(OBJ_EXT) \
    leqsol_mod.$(OBJ_EXT) \
    machine_mod.$(OBJ_EXT) \
    machine.$(OBJ_EXT) \
    stiff_chem_mod.$(OBJ_EXT) \
    mflux_mod.$(OBJ_EXT) \
    mfix_netcdf_mod.$(OBJ_EXT) \
    output_mod.$(OBJ_EXT) \
    parallel_mod.$(OBJ_EXT) \
    param1_mod.$(OBJ_EXT) \
    param_mod.$(OBJ_EXT) \
    parse_mod.$(OBJ_EXT) \
    pgcor_mod.$(OBJ_EXT) \
    physprop_mod.$(OBJ_EXT) \
    post3d_mod.$(OBJ_EXT) \
    pscor_mod.$(OBJ_EXT) \
    ps_mod.$(OBJ_EXT) \
    residual_mod.$(OBJ_EXT) \
    run_mod.$(OBJ_EXT) \
    rxns_mod.$(OBJ_EXT) \
    scalars_mod.$(OBJ_EXT) \
    scales_mod.$(OBJ_EXT) \
    tau_g_mod.$(OBJ_EXT) \
    tau_s_mod.$(OBJ_EXT) \
    time_cpu_mod.$(OBJ_EXT) \
    tmp_array1_mod.$(OBJ_EXT) \
    tmp_array_mod.$(OBJ_EXT) \
    toleranc_mod.$(OBJ_EXT) \
    trace_mod.$(OBJ_EXT) \
    turb_mod.$(OBJ_EXT) \
    ur_facs_mod.$(OBJ_EXT) \
    usr_mod.$(OBJ_EXT) \
    visc_g_mod.$(OBJ_EXT) \
    visc_s_mod.$(OBJ_EXT) \
    vshear_mod.$(OBJ_EXT) \
    xsi_array_mod.$(OBJ_EXT) \
    allocate_arrays.$(OBJ_EXT) \
    any_more_data.$(OBJ_EXT) \
    calc_cell2.$(OBJ_EXT) \
    calc_corr_01.$(OBJ_EXT) \
    calc_corr_type_1.$(OBJ_EXT) \
    calc_distance.$(OBJ_EXT) \
    calc_ep_g.$(OBJ_EXT) \
    calc_mu_s.$(OBJ_EXT) \
    calc_mw.$(OBJ_EXT) \
    calc_p_star.$(OBJ_EXT) \
    calc_quantities.$(OBJ_EXT) \
    calc_ro_g.$(OBJ_EXT) \
    calc_vol.$(OBJ_EXT) \
    check_data_03.$(OBJ_EXT) \
    check_data_04.$(OBJ_EXT) \
    check_data_05.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    compare.$(OBJ_EXT) \
    deallocate_arrays.$(OBJ_EXT) \
    progress_bar_mod.$(OBJ_EXT) \
    quadric_mod.$(OBJ_EXT) \
    stl_mod.$(OBJ_EXT) \
    vtk_mod.$(OBJ_EXT) \
    cutcell_mod.$(OBJ_EXT) \
    polygon_mod.$(OBJ_EXT) \
    dashboard_mod.$(OBJ_EXT) \
    cartesian_grid_init_namelist.$(OBJ_EXT) \
    eosg.$(OBJ_EXT) \
    error_routine.$(OBJ_EXT) \
    examine_data.$(OBJ_EXT) \
    exit.$(OBJ_EXT) \
    f_init_data.$(OBJ_EXT) \
    file_handle.$(OBJ_EXT) \
    finit.$(OBJ_EXT) \
    flow_gx.$(OBJ_EXT) \
    flow_gy.$(OBJ_EXT) \
    flow_gz.$(OBJ_EXT) \
    flow_sx.$(OBJ_EXT) \
    flow_sy.$(OBJ_EXT) \
    flow_sz.$(OBJ_EXT) \
    g_0.$(OBJ_EXT) \
    gas_flux.$(OBJ_EXT) \
    get_file_name.$(OBJ_EXT) \
    get_file_status.$(OBJ_EXT) \
    get_index.$(OBJ_EXT) \
    get_location.$(OBJ_EXT) \
    get_mu_s.$(OBJ_EXT) \
    get_same_time.$(OBJ_EXT) \
    get_selection.$(OBJ_EXT) \
    get_substr.$(OBJ_EXT) \
    granular_qty.$(OBJ_EXT) \
    header_main.$(OBJ_EXT) \
    ik_avg.$(OBJ_EXT) \
    ik_avg_out.$(OBJ_EXT) \
    init_namelist.$(OBJ_EXT) \
    interp_res.$(OBJ_EXT) \
    kintheory_energy_dissipation_ss.$(OBJ_EXT) \
    line_too_big.$(OBJ_EXT) \
    main_f.$(OBJ_EXT) \
    make_upper_case.$(OBJ_EXT) \
    open_file.$(OBJ_EXT) \
    out_from_res.$(OBJ_EXT) \
    out_from_spx.$(OBJ_EXT) \
    out_spec_time.$(OBJ_EXT) \
    out_time.$(OBJ_EXT) \
    parse_line.$(OBJ_EXT) \
    parse_rxn.$(OBJ_EXT) \
    print_out.$(OBJ_EXT) \
    read_namelist.$(OBJ_EXT) \
    read_res0.$(OBJ_EXT) \
    read_res1.$(OBJ_EXT) \
    read_spx0.$(OBJ_EXT) \
    read_spx1.$(OBJ_EXT) \
    remove_comment.$(OBJ_EXT) \
    res_from_spx.$(OBJ_EXT) \
    seek_comment.$(OBJ_EXT) \
    seek_end.$(OBJ_EXT) \
    seek_time.$(OBJ_EXT) \
    select_spx_rec.$(OBJ_EXT) \
    set_constants.$(OBJ_EXT) \
    set_dollar.$(OBJ_EXT) \
    set_geometry.$(OBJ_EXT) \
    set_increments.$(OBJ_EXT) \
    set_index1.$(OBJ_EXT) \
    set_index1a.$(OBJ_EXT) \
    set_max2.$(OBJ_EXT) \
    set_read_spx.$(OBJ_EXT) \
    shift_dxyz.$(OBJ_EXT) \
    sol_flux.$(OBJ_EXT) \
    strcmp.$(OBJ_EXT) \
    streqs.$(OBJ_EXT) \
    time_avg.$(OBJ_EXT) \
    usr_init_namelist.$(OBJ_EXT) \
    usr_post.$(OBJ_EXT) \
    usr_write_out1.$(OBJ_EXT) \
    write_out1.$(OBJ_EXT) \
    write_res0.$(OBJ_EXT) \
    write_res1.$(OBJ_EXT) \
    write_spx0.$(OBJ_EXT) \
    write_spx1.$(OBJ_EXT) \
    usr_input.$(OBJ_EXT) \
    compar_mod.$(OBJ_EXT) \
    mpi_mod.$(OBJ_EXT) \
    dbg_util_mod.$(OBJ_EXT) \
    parallel_mpi_mod.$(OBJ_EXT) \
    debug_mod.$(OBJ_EXT) \
    gridmap_mod.$(OBJ_EXT) \
    mpi_utility_mod.$(OBJ_EXT) \
    sendrecv_mod.$(OBJ_EXT) \
    boundfunijk_mod.$(OBJ_EXT) \
    write_error.$(OBJ_EXT) \
    discretelement_mod.$(OBJ_EXT) \
    ghdtheory_mod.$(OBJ_EXT) \
    transport_coeff_ghd.$(OBJ_EXT) \
    ghd.$(OBJ_EXT) \
    cooling_rate.$(OBJ_EXT) \
    cooling_rate_tc.$(OBJ_EXT) \
    pressure.$(OBJ_EXT) \
    bulk_viscosity.$(OBJ_EXT) \
    shear_viscosity.$(OBJ_EXT) \
    thermal_diffusivity.$(OBJ_EXT) \
    mass_mobility.$(OBJ_EXT) \
    thermal_conductivity.$(OBJ_EXT) \
    thermal_mobility.$(OBJ_EXT) \
    ordinary_diff.$(OBJ_EXT) \
    dufour_coeff.$(OBJ_EXT) \
    chi_ij_GHD.$(OBJ_EXT) \
    des_init_namelist.$(OBJ_EXT) \
    ornl_header.$(OBJ_EXT) \
    ornl_util.$(OBJ_EXT) \
    ornl_stats.$(OBJ_EXT) \
    ornl_stats_c.$(OBJ_EXT) \
    ornl_corr.$(OBJ_EXT) \
    ornl_corr_c.$(OBJ_EXT) \
    ornl_pca.$(OBJ_EXT) \
    ornl_ft.$(OBJ_EXT) \
    ornl_ft_c.$(OBJ_EXT) \
    ornl_filt.$(OBJ_EXT) \
    ornl_filt_c.$(OBJ_EXT) \
    ornl_zone.$(OBJ_EXT) \
    ornl_sym.$(OBJ_EXT) \
    post_precision_mod.$(OBJ_EXT) \
    paralleldata_mod.$(OBJ_EXT) \
    usr_input_mod.$(OBJ_EXT) \
    mfix_pic_mod.$(OBJ_EXT) \
    qmom_kinetic_equation_mod.$(OBJ_EXT) \
    qmomk_parameters_mod.$(OBJ_EXT) \
    qmomk_init_namelist.$(OBJ_EXT) \
    des_ic_mod.$(OBJ_EXT) \
    des_thermo_mod.$(OBJ_EXT) \
    des_rxns_mod.$(OBJ_EXT) \
    calc_trd_s.$(OBJ_EXT) \
    get_delh.$(OBJ_EXT) \
    define_quadrics.$(OBJ_EXT) \
    check_data_cartesian.$(OBJ_EXT) \
    get_poly_data.$(OBJ_EXT) \
    get_stl_data.$(OBJ_EXT) \
    get_connectivity.$(OBJ_EXT) \
    cut_cell_preprocessing.$(OBJ_EXT) \
    eval_usr_fct.$(OBJ_EXT) \
    allocate_cut_cell_arrays.$(OBJ_EXT) \
    allocate_dummy_cut_cell_arrays.$(OBJ_EXT) \
    calc_vort_out.$(OBJ_EXT) \
    deallocate_cut_cell_arrays.$(OBJ_EXT) \
    get_alpha.$(OBJ_EXT) \
    get_cut_cell_flags.$(OBJ_EXT) \
    get_cut_cell_volume_area.$(OBJ_EXT) \
    get_master.$(OBJ_EXT) \
    set_Odxyz.$(OBJ_EXT) \
    vtk_out.$(OBJ_EXT) \
    write_progress_bar.$(OBJ_EXT) \
    dmp_cartesian.$(OBJ_EXT) \
    set_geometry1.$(OBJ_EXT) \
    read_database.$(OBJ_EXT) \
    readTherm.$(OBJ_EXT) \
    get_values.$(OBJ_EXT) \
    rxn_com_mod.$(OBJ_EXT) \
    get_bc_area.$(OBJ_EXT) \
    flow_to_vel.$(OBJ_EXT) \
  -o post_mfix $(LIB_FLAGS)
  
ambm.mod : ../model/ambm_mod.f \
            compar.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/ambm_mod.f 
bc.mod : ../model/bc_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/bc_mod.f 
cdist.mod : ../model/cdist_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cdist_mod.f 
constant.mod : ../model/constant_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/constant_mod.f 
cont.mod : ../model/cont_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cont_mod.f 
correl.mod : correl_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) correl_mod.f 
des_bc.mod : ../model/des/des_bc_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/des_bc_mod.f 
drag.mod : ../model/drag_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/drag_mod.f 
energy.mod : ../model/energy_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/energy_mod.f 
fldvar.mod : ../model/fldvar_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/fldvar_mod.f 
funits.mod : ../model/funits_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/funits_mod.f 
geometry.mod : ../model/geometry_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/geometry_mod.f 
ic.mod : ../model/ic_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/ic_mod.f 
indices.mod : ../model/indices_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/indices_mod.f 
is.mod : ../model/is_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/is_mod.f 
kintheory.mod : ../model/kintheory_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/kintheory_mod.f 
leqsol.mod : ../model/leqsol_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/leqsol_mod.f 
machine.mod : ../model/machine_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/machine_mod.f 
stiff_chem.mod : chem/stiff_chem_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) chem/stiff_chem_mod.f 
mflux.mod : ../model/mflux_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/mflux_mod.f 
mfix_netcdf.mod : ../model/mfix_netcdf_mod.f \
            MFIX_netcdf_constants.fi                                     \
            MFIX_netcdf_overloads.fi                                     \
            MFIX_netcdf_variables.fi                                     \
            MFIX_netcdf_misc.fi                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/mfix_netcdf_mod.f 
output.mod : ../model/output_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/output_mod.f 
parallel.mod : ../model/parallel_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/parallel_mod.f 
param1.mod : ../model/param1_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/param1_mod.f 
param.mod : ../model/param_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/param_mod.f 
parse.mod : ../model/parse_mod.f \
            param.mod \
            param1.mod \
            funits.mod \
            compar.mod \
            rxn_com.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/parse_mod.f 
pgcor.mod : ../model/pgcor_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/pgcor_mod.f 
physprop.mod : ../model/physprop_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/physprop_mod.f 
post3d.mod : post3d_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) post3d_mod.f 
pscor.mod : ../model/pscor_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/pscor_mod.f 
ps.mod : ../model/ps_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/ps_mod.f 
residual.mod : ../model/residual_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/residual_mod.f 
run.mod : ../model/run_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/run_mod.f 
rxns.mod : ../model/rxns_mod.f \
            param.mod \
            param1.mod \
            rxn_com.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/rxns_mod.f 
scalars.mod : ../model/scalars_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/scalars_mod.f 
scales.mod : ../model/scales_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/scales_mod.f 
tau_g.mod : ../model/tau_g_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/tau_g_mod.f 
tau_s.mod : ../model/tau_s_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/tau_s_mod.f 
time_cpu.mod : ../model/time_cpu_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/time_cpu_mod.f 
tmp_array1.mod : ../model/tmp_array1_mod.f \
            compar.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/tmp_array1_mod.f 
tmp_array.mod : ../model/tmp_array_mod.f \
            compar.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/tmp_array_mod.f 
toleranc.mod : ../model/toleranc_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/toleranc_mod.f 
trace.mod : ../model/trace_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/trace_mod.f 
turb.mod : ../model/turb_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/turb_mod.f 
ur_facs.mod : ../model/ur_facs_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/ur_facs_mod.f 
usr.mod : ../model/usr_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/usr_mod.f 
visc_g.mod : ../model/visc_g_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/visc_g_mod.f 
visc_s.mod : ../model/visc_s_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/visc_s_mod.f 
vshear.mod : ../model/vshear_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/vshear_mod.f 
xsi_array.mod : ../model/xsi_array_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/xsi_array_mod.f 
progress_bar.mod : ../model/cartesian_grid/progress_bar_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/progress_bar_mod.f 
quadric.mod : ../model/cartesian_grid/quadric_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/quadric_mod.f 
stl.mod : ../model/cartesian_grid/stl_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/stl_mod.f 
vtk.mod : ../model/cartesian_grid/vtk_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/vtk_mod.f 
cutcell.mod : ../model/cartesian_grid/cutcell_mod.f \
            param.mod \
            param1.mod \
            progress_bar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/cutcell_mod.f 
polygon.mod : ../model/cartesian_grid/polygon_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/polygon_mod.f 
dashboard.mod : ../model/cartesian_grid/dashboard_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/dashboard_mod.f 
compar.mod : ../model/dmp_modules/mpi_donothing/compar_mod.f \
            mpi.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/compar_mod.f 
mpi.mod : ../model/dmp_modules/mpi_donothing/mpi_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/mpi_mod.f 
dbg_util.mod : ../model/dmp_modules/mpi_donothing/dbg_util_mod.f \
            compar.mod \
            geometry.mod \
            parallel_mpi.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/dbg_util_mod.f 
parallel_mpi.mod : ../model/dmp_modules/mpi_donothing/parallel_mpi_mod.f \
            geometry.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/parallel_mpi_mod.f 
debug.mod : ../model/dmp_modules/mpi_donothing/debug_mod.f \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/debug_mod.f 
gridmap.mod : ../model/dmp_modules/mpi_donothing/gridmap_mod.f \
            mpi_utility.mod \
            parallel_mpi.mod \
            geometry.mod \
            sendrecv.mod \
            compar.mod \
            run.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/gridmap_mod.f 
mpi_utility.mod : ../model/dmp_modules/mpi_donothing/mpi_utility_mod.f \
            geometry.mod \
            compar.mod \
            parallel_mpi.mod \
            debug.mod \
            indices.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/mpi_utility_mod.f 
sendrecv.mod : ../model/dmp_modules/mpi_donothing/sendrecv_mod.f \
            parallel_mpi.mod \
            debug.mod \
            geometry.mod \
            compar.mod \
            indices.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/sendrecv_mod.f 
boundfunijk.mod : ../model/boundfunijk_mod.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/boundfunijk_mod.f 
discretelement.mod : ../model/des/discretelement_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/discretelement_mod.f 
ghdtheory.mod : ../model/GhdTheory/ghdtheory_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/ghdtheory_mod.f 
post_precision.mod : post_precision_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) post_precision_mod.f 
paralleldata.mod : paralleldata_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) paralleldata_mod.f 
usr_input.mod : usr_input_mod.f \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_input_mod.f 
mfix_pic.mod : ../model/des/mfix_pic_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/mfix_pic_mod.f 
qmom_kinetic_equation.mod : ../model/qmomk/qmom_kinetic_equation_mod.f \
            param.mod \
            param1.mod \
            qmomk_parameters.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/qmomk/qmom_kinetic_equation_mod.f 
qmomk_parameters.mod : ../model/qmomk/qmomk_parameters_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/qmomk/qmomk_parameters_mod.f 
des_ic.mod : ../model/des/des_ic_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/des_ic_mod.f 
des_thermo.mod : ../model/des/des_thermo_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/des_thermo_mod.f 
des_rxns.mod : ../model/des/des_rxns_mod.f \
            param.mod \
            rxn_com.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/des_rxns_mod.f 
rxn_com.mod : ../model/rxn_com_mod.f \
            param.mod \
            param1.mod \
            compar.mod \
            funits.mod \
            mfix_directory_path.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/rxn_com_mod.f 
machine.$(OBJ_EXT) : machine.f \
            machine.mod \
            param.mod \
            run.mod \
            funits.mod 
allocate_arrays.$(OBJ_EXT) : ../model/allocate_arrays.f \
            param.mod \
            param1.mod \
            ambm.mod \
            cont.mod \
            drag.mod \
            energy.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            physprop.mod \
            pscor.mod \
            residual.mod \
            rxns.mod \
            run.mod \
            scalars.mod \
            turb.mod \
            tau_g.mod \
            tau_s.mod \
            tmp_array.mod \
            tmp_array1.mod \
            trace.mod \
            visc_g.mod \
            visc_s.mod \
            xsi_array.mod \
            vshear.mod \
            mflux.mod \
            ghdtheory.mod \
            kintheory.mod \
            cdist.mod \
            des_rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/allocate_arrays.f 
any_more_data.$(OBJ_EXT) : any_more_data.f \
            param.mod \
            param1.mod 
calc_cell2.$(OBJ_EXT) : calc_cell2.f \
            param.mod \
            param1.mod 
calc_corr_01.$(OBJ_EXT) : calc_corr_01.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            correl.mod \
            compar.mod \
            xforms.inc                                                   \
            function.inc                                                
calc_corr_type_1.$(OBJ_EXT) : calc_corr_type_1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            post3d.mod \
            physprop.mod \
            correl.mod \
            compar.mod \
            xforms.inc                                                   \
            function.inc                                                
calc_distance.$(OBJ_EXT) : calc_distance.f 
calc_ep_g.$(OBJ_EXT) : calc_ep_g.f \
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
calc_mu_s.$(OBJ_EXT) : ../model/calc_mu_s.f \
            run.mod \
            vshear.mod \
            visc_s.mod \
            physprop.mod \
            constant.mod \
            fldvar.mod \
            compar.mod \
            indices.mod \
            geometry.mod \
            qmom_kinetic_equation.mod \
            param.mod \
            param1.mod \
            trace.mod \
            toleranc.mod \
            turb.mod \
            drag.mod \
            kintheory.mod \
            ur_facs.mod \
            cutcell.mod \
            parallel.mod \
            visc_g.mod \
            is.mod \
            sendrecv.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            s_pr1.inc                                                    \
            s_pr2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/calc_mu_s.f 
calc_mw.$(OBJ_EXT) : ../model/calc_mw.f \
            param.mod \
            param1.mod \
            toleranc.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/calc_mw.f 
calc_p_star.$(OBJ_EXT) : ../model/calc_p_star.f \
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
            run.mod \
            visc_s.mod \
            fldvar.mod \
            toleranc.mod \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/calc_p_star.f 
calc_quantities.$(OBJ_EXT) : calc_quantities.f \
            param.mod \
            param1.mod \
            post3d.mod 
calc_ro_g.$(OBJ_EXT) : calc_ro_g.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            function.inc                                                
calc_vol.$(OBJ_EXT) : calc_vol.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            fldvar.mod \
            physprop.mod \
            post3d.mod \
            compar.mod \
            function.inc                                                
check_data_03.$(OBJ_EXT) : ../model/check_data_03.f \
            param.mod \
            param1.mod \
            geometry.mod \
            bc.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/check_data_03.f 
check_data_04.$(OBJ_EXT) : ../model/check_data_04.f \
            param.mod \
            param1.mod \
            run.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            discretelement.mod \
            funits.mod \
            mfix_pic.mod \
            compar.mod \
            rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/check_data_04.f 
check_data_05.$(OBJ_EXT) : ../model/check_data_05.f \
            compar.mod \
            param.mod \
            param1.mod \
            physprop.mod \
            funits.mod \
            run.mod \
            indices.mod \
            rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/check_data_05.f 
check_one_axis.$(OBJ_EXT) : ../model/check_one_axis.f \
            param.mod \
            param1.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/check_one_axis.f 
compare.$(OBJ_EXT) : ../model/compare.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/compare.f 
deallocate_arrays.$(OBJ_EXT) : deallocate_arrays.f \
            param.mod \
            param1.mod \
            ambm.mod \
            cont.mod \
            drag.mod \
            energy.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            physprop.mod \
            pscor.mod \
            residual.mod \
            rxns.mod \
            run.mod \
            scalars.mod \
            tau_g.mod \
            tau_s.mod \
            tmp_array.mod \
            tmp_array1.mod \
            trace.mod \
            visc_g.mod \
            visc_s.mod \
            xsi_array.mod \
            vshear.mod \
            mflux.mod \
            kintheory.mod \
            ghdtheory.mod 
cartesian_grid_init_namelist.$(OBJ_EXT) : ../model/cartesian_grid/cartesian_grid_init_namelist.f \
            param1.mod \
            quadric.mod \
            cutcell.mod \
            polygon.mod \
            vtk.mod \
            progress_bar.mod \
            dashboard.mod \
            stl.mod \
            cartesian_grid/cartesian_grid_namelist.inc                  
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/cartesian_grid_init_namelist.f 
eosg.$(OBJ_EXT) : ../model/eosg.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            scales.mod \
            sc_p_g1.inc                                                  \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/eosg.f 
error_routine.$(OBJ_EXT) : ../model/error_routine.f \
            funits.mod \
            compar.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/error_routine.f 
examine_data.$(OBJ_EXT) : examine_data.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            indices.mod \
            run.mod \
            geometry.mod \
            post3d.mod \
            rxns.mod \
            scalars.mod \
            compar.mod \
            post_precision.mod \
            xforms.inc                                                   \
            function.inc                                                
exit.$(OBJ_EXT) : ../model/exit.f \
            funits.mod \
            compar.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/exit.f 
f_init_data.$(OBJ_EXT) : f_init_data.f \
            param.mod \
            param1.mod \
            geometry.mod \
            physprop.mod \
            constant.mod \
            xforms.inc                                                  
file_handle.$(OBJ_EXT) : file_handle.f \
            machine.mod \
            param.mod \
            param1.mod \
            post3d.mod \
            funits.mod \
            compar.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            function.inc                                                
finit.$(OBJ_EXT) : finit.f \
            param.mod \
            param1.mod \
            run.mod \
            post3d.mod \
            geometry.mod \
            indices.mod \
            fldvar.mod \
            physprop.mod \
            constant.mod \
            funits.mod \
            parallel_mpi.mod \
            gridmap.mod \
            cdist.mod \
            paralleldata.mod \
            scalars.mod \
            rxns.mod \
            drag.mod \
            energy.mod \
            xforms.inc                                                  
flow_gx.$(OBJ_EXT) : flow_gx.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            function.inc                                                
flow_gy.$(OBJ_EXT) : flow_gy.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            function.inc                                                
flow_gz.$(OBJ_EXT) : flow_gz.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            function.inc                                                
flow_sx.$(OBJ_EXT) : flow_sx.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
flow_sy.$(OBJ_EXT) : flow_sy.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
flow_sz.$(OBJ_EXT) : flow_sz.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
g_0.$(OBJ_EXT) : ../model/g_0.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            visc_s.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/g_0.f 
gas_flux.$(OBJ_EXT) : gas_flux.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            post3d.mod \
            physprop.mod \
            compar.mod \
            xforms.inc                                                   \
            function.inc                                                
get_file_name.$(OBJ_EXT) : get_file_name.f 
get_file_status.$(OBJ_EXT) : get_file_status.f 
get_index.$(OBJ_EXT) : get_index.f 
get_location.$(OBJ_EXT) : get_location.f \
            param.mod \
            param1.mod \
            geometry.mod \
            post3d.mod 
get_mu_s.$(OBJ_EXT) : get_mu_s.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            visc_s.mod \
            run.mod \
            constant.mod \
            funits.mod \
            compar.mod \
            xforms.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
get_same_time.$(OBJ_EXT) : get_same_time.f \
            param.mod \
            param1.mod 
get_selection.$(OBJ_EXT) : get_selection.f 
get_substr.$(OBJ_EXT) : get_substr.f 
granular_qty.$(OBJ_EXT) : granular_qty.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            visc_s.mod \
            run.mod \
            constant.mod \
            funits.mod \
            compar.mod \
            xforms.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
header_main.$(OBJ_EXT) : header_main.f \
            param.mod \
            param1.mod \
            post3d.mod 
ik_avg.$(OBJ_EXT) : ik_avg.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            geometry.mod \
            compar.mod \
            function.inc                                                
ik_avg_out.$(OBJ_EXT) : ik_avg_out.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            post3d.mod 
init_namelist.$(OBJ_EXT) : ../model/init_namelist.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            geometry.mod \
            ic.mod \
            bc.mod \
            ps.mod \
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
            scalars.mod \
            compar.mod \
            parallel.mod \
            cdist.mod \
            stiff_chem.mod \
            namelist.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/init_namelist.f 
interp_res.$(OBJ_EXT) : interp_res.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            sendrecv.mod \
            drag.mod \
            kintheory.mod \
            tmp_array.mod \
            xforms.inc                                                   \
            function.inc                                                
kintheory_energy_dissipation_ss.$(OBJ_EXT) : ../model/kintheory_energy_dissipation_ss.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            constant.mod \
            toleranc.mod \
            kintheory.mod \
            function.inc                                                   \
            fun_avg1.inc                                                   \
            fun_avg2.inc                                                   \
            ep_s1.inc                                                   \
            ep_s2.inc                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/kintheory_energy_dissipation_ss.f
line_too_big.$(OBJ_EXT) : ../model/line_too_big.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/line_too_big.f 
main_f.$(OBJ_EXT) : main_f.f \
            xforms.inc                                                  
make_upper_case.$(OBJ_EXT) : ../model/make_upper_case.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/make_upper_case.f 
open_file.$(OBJ_EXT) : ../model/open_file.f \
            cdist.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/open_file.f 
out_from_res.$(OBJ_EXT) : out_from_res.f \
            param.mod \
            param1.mod \
            run.mod \
            funits.mod \
            xforms.inc                                                  
out_from_spx.$(OBJ_EXT) : out_from_spx.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            funits.mod \
            post3d.mod \
            xforms.inc                                                  
out_spec_time.$(OBJ_EXT) : out_spec_time.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            physprop.mod \
            indices.mod \
            geometry.mod \
            compar.mod \
            function.inc                                                
out_time.$(OBJ_EXT) : out_time.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            physprop.mod \
            indices.mod \
            geometry.mod \
            compar.mod \
            function.inc                                                
parse_line.$(OBJ_EXT) : ../model/parse_line.f \
            compar.mod \
            des_rxns.mod \
            param.mod \
            param1.mod \
            parse.mod \
            rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/parse_line.f 
parse_rxn.$(OBJ_EXT) : ../model/parse_rxn.f \
            compar.mod \
            funits.mod \
            param.mod \
            param1.mod \
            parse.mod \
            rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/parse_rxn.f 
print_out.$(OBJ_EXT) : print_out.f \
            param.mod \
            param1.mod \
            post3d.mod \
            xforms.inc                                                  
read_namelist.$(OBJ_EXT) : ../model/read_namelist.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            geometry.mod \
            ic.mod \
            is.mod \
            bc.mod \
            ps.mod \
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
            scalars.mod \
            compar.mod \
            parallel.mod \
            discretelement.mod \
            mfix_pic.mod \
            usr.mod \
            des_bc.mod \
            des_ic.mod \
            des_thermo.mod \
            des_rxns.mod \
            stiff_chem.mod \
            cdist.mod \
            quadric.mod \
            cutcell.mod \
            vtk.mod \
            polygon.mod \
            dashboard.mod \
            stl.mod \
            qmom_kinetic_equation.mod \
            usrnlst.inc                                                  \
            namelist.inc                                                 \
            des/desnamelist.inc                                          \
            cartesian_grid/cartesian_grid_namelist.inc                   \
            qmomk/qmomknamelist.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/read_namelist.f 
read_res0.$(OBJ_EXT) : ../model/read_res0.f \
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
            scalars.mod \
            rxns.mod \
            compar.mod \
            mpi_utility.mod \
            fldvar.mod \
            stiff_chem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/read_res0.f 
read_res1.$(OBJ_EXT) : ../model/read_res1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            rxns.mod \
            scalars.mod \
            funits.mod \
            energy.mod \
            compar.mod \
            cdist.mod \
            mpi_utility.mod \
            sendrecv.mod \
            mfix_netcdf.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/read_res1.f 
read_spx0.$(OBJ_EXT) : read_spx0.f \
            param.mod \
            param1.mod \
            run.mod \
            funits.mod \
            post3d.mod 
read_spx1.$(OBJ_EXT) : read_spx1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            post3d.mod \
            scalars.mod \
            rxns.mod \
            machine.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
remove_comment.$(OBJ_EXT) : ../model/remove_comment.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/remove_comment.f 
res_from_spx.$(OBJ_EXT) : res_from_spx.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            funits.mod \
            post3d.mod \
            physprop.mod \
            fldvar.mod 
seek_comment.$(OBJ_EXT) : ../model/seek_comment.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/seek_comment.f 
seek_end.$(OBJ_EXT) : ../model/seek_end.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/seek_end.f 
seek_time.$(OBJ_EXT) : seek_time.f \
            param.mod \
            param1.mod \
            post3d.mod \
            funits.mod 
select_spx_rec.$(OBJ_EXT) : select_spx_rec.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            machine.mod \
            funits.mod \
            post3d.mod \
            physprop.mod \
            fldvar.mod \
            xforms.inc                                                  
set_constants.$(OBJ_EXT) : ../model/set_constants.f \
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
            drag.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_constants.f 
set_dollar.$(OBJ_EXT) : set_dollar.f 
set_geometry.$(OBJ_EXT) : ../model/set_geometry.f \
            param.mod \
            param1.mod \
            run.mod \
            geometry.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_geometry.f 
set_increments.$(OBJ_EXT) : ../model/set_increments.f \
            param.mod \
            param1.mod \
            indices.mod \
            geometry.mod \
            compar.mod \
            physprop.mod \
            fldvar.mod \
            funits.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_increments.f 
set_index1.$(OBJ_EXT) : ../model/set_index1.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            constant.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_index1.f 
set_index1a.$(OBJ_EXT) : ../model/set_index1a.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            boundfunijk.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_index1a.f 
set_max2.$(OBJ_EXT) : ../model/set_max2.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_max2.f 
set_read_spx.$(OBJ_EXT) : set_read_spx.f 
shift_dxyz.$(OBJ_EXT) : ../model/shift_dxyz.f \
            param.mod \
            param1.mod \
            geometry.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/shift_dxyz.f 
sol_flux.$(OBJ_EXT) : sol_flux.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            post3d.mod \
            physprop.mod \
            compar.mod \
            xforms.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
strcmp.$(OBJ_EXT) : strcmp.f 
streqs.$(OBJ_EXT) : streqs.f 
time_avg.$(OBJ_EXT) : time_avg.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            machine.mod \
            funits.mod \
            post3d.mod \
            physprop.mod \
            fldvar.mod \
            scalars.mod \
            rxns.mod \
            xforms.inc                                                  
usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f \
            usrnlst.inc                                                 
usr_post.$(OBJ_EXT) : usr_post.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            post3d.mod \
            physprop.mod \
            compar.mod \
            function.inc                                                
usr_write_out1.$(OBJ_EXT) : usr_write_out1.f 
write_out1.$(OBJ_EXT) : ../model/write_out1.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            rxns.mod \
            compar.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_out1.f 
write_res0.$(OBJ_EXT) : ../model/write_res0.f \
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
            scalars.mod \
            rxns.mod \
            ur_facs.mod \
            leqsol.mod \
            toleranc.mod \
            cdist.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            stiff_chem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_res0.f 
write_res1.$(OBJ_EXT) : ../model/write_res1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            scalars.mod \
            rxns.mod \
            funits.mod \
            output.mod \
            energy.mod \
            cdist.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            mfix_netcdf.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_res1.f 
write_spx0.$(OBJ_EXT) : ../model/write_spx0.f \
            param.mod \
            param1.mod \
            run.mod \
            funits.mod \
            cdist.mod \
            compar.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_spx0.f 
write_spx1.$(OBJ_EXT) : ../model/write_spx1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            scalars.mod \
            output.mod \
            rxns.mod \
            cdist.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            mfix_netcdf.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_spx1.f 
usr_input.$(OBJ_EXT) : usr_input.f \
            usr_input.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            compar.mod \
            constant.mod \
            post3d.mod \
            xforms.inc                                                   \
            function.inc                                                
write_error.$(OBJ_EXT) : ../model/write_error.f \
            param.mod \
            param1.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_error.f 
transport_coeff_ghd.$(OBJ_EXT) : ../model/GhdTheory/transport_coeff_ghd.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            visc_s.mod \
            ghdtheory.mod \
            physprop.mod \
            run.mod \
            constant.mod \
            toleranc.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/transport_coeff_ghd.f 
ghd.$(OBJ_EXT) : ../model/GhdTheory/ghd.f \
            drag.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/ghd.f 
cooling_rate.$(OBJ_EXT) : ../model/GhdTheory/cooling_rate.f \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/cooling_rate.f 
cooling_rate_tc.$(OBJ_EXT) : ../model/GhdTheory/cooling_rate_tc.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/cooling_rate_tc.f 
pressure.$(OBJ_EXT) : ../model/GhdTheory/pressure.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/pressure.f 
bulk_viscosity.$(OBJ_EXT) : ../model/GhdTheory/bulk_viscosity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/bulk_viscosity.f 
shear_viscosity.$(OBJ_EXT) : ../model/GhdTheory/shear_viscosity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/shear_viscosity.f 
thermal_diffusivity.$(OBJ_EXT) : ../model/GhdTheory/thermal_diffusivity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/thermal_diffusivity.f 
mass_mobility.$(OBJ_EXT) : ../model/GhdTheory/mass_mobility.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/mass_mobility.f 
thermal_conductivity.$(OBJ_EXT) : ../model/GhdTheory/thermal_conductivity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/thermal_conductivity.f 
thermal_mobility.$(OBJ_EXT) : ../model/GhdTheory/thermal_mobility.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/thermal_mobility.f 
ordinary_diff.$(OBJ_EXT) : ../model/GhdTheory/ordinary_diff.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/ordinary_diff.f 
dufour_coeff.$(OBJ_EXT) : ../model/GhdTheory/dufour_coeff.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/dufour_coeff.f 
chi_ij_GHD.$(OBJ_EXT) : ../model/GhdTheory/chi_ij_GHD.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/chi_ij_GHD.f 
des_init_namelist.$(OBJ_EXT) : ../model/des/des_init_namelist.f \
            param1.mod \
            discretelement.mod \
            mfix_pic.mod \
            des_bc.mod \
            des_ic.mod \
            des_thermo.mod \
            des_rxns.mod \
            des/desnamelist.inc                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/des_init_namelist.f 
ornl_header.$(OBJ_EXT) : ornl_header.f \
            geometry.mod \
            run.mod \
            usr_input.mod 
ornl_util.$(OBJ_EXT) : ornl_util.f 
ornl_stats.$(OBJ_EXT) : ornl_stats.f 
ornl_stats_c.$(OBJ_EXT) : ornl_stats_c.f \
            geometry.mod \
            fldvar.mod \
            run.mod \
            indices.mod \
            compar.mod \
            usr_input.mod \
            function.inc                                                
ornl_corr.$(OBJ_EXT) : ornl_corr.f 
ornl_corr_c.$(OBJ_EXT) : ornl_corr_c.f \
            geometry.mod \
            fldvar.mod \
            run.mod \
            indices.mod \
            compar.mod \
            usr_input.mod \
            function.inc                                                
ornl_pca.$(OBJ_EXT) : ornl_pca.f 
ornl_ft.$(OBJ_EXT) : ornl_ft.f 
ornl_ft_c.$(OBJ_EXT) : ornl_ft_c.f \
            usr_input.mod 
ornl_filt.$(OBJ_EXT) : ornl_filt.f 
ornl_filt_c.$(OBJ_EXT) : ornl_filt_c.f \
            usr_input.mod 
ornl_zone.$(OBJ_EXT) : ornl_zone.f \
            geometry.mod \
            fldvar.mod \
            run.mod \
            indices.mod \
            compar.mod \
            usr_input.mod \
            function.inc                                                
ornl_sym.$(OBJ_EXT) : ornl_sym.f 
qmomk_init_namelist.$(OBJ_EXT) : ../model/qmomk/qmomk_init_namelist.f \
            param1.mod \
            qmom_kinetic_equation.mod \
            qmomk/qmomknamelist.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/qmomk/qmomk_init_namelist.f 
calc_trd_s.$(OBJ_EXT) : ../model/calc_trd_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            sendrecv.mod \
            bc.mod \
            cutcell.mod \
            quadric.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/calc_trd_s.f 
get_delh.$(OBJ_EXT) : ../model/cartesian_grid/get_delh.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_delh.f 
define_quadrics.$(OBJ_EXT) : ../model/cartesian_grid/define_quadrics.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            vtk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/define_quadrics.f 
check_data_cartesian.$(OBJ_EXT) : ../model/cartesian_grid/check_data_cartesian.f \
            param.mod \
            param1.mod \
            constant.mod \
            run.mod \
            physprop.mod \
            indices.mod \
            scalars.mod \
            funits.mod \
            leqsol.mod \
            compar.mod \
            mpi_utility.mod \
            bc.mod \
            discretelement.mod \
            cutcell.mod \
            quadric.mod \
            vtk.mod \
            polygon.mod \
            dashboard.mod \
            stl.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/check_data_cartesian.f 
get_poly_data.$(OBJ_EXT) : ../model/cartesian_grid/get_poly_data.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            rxns.mod \
            compar.mod \
            mpi_utility.mod \
            progress_bar.mod \
            polygon.mod \
            parallel.mod \
            constant.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_poly_data.f 
get_stl_data.$(OBJ_EXT) : ../model/cartesian_grid/get_stl_data.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            rxns.mod \
            compar.mod \
            mpi_utility.mod \
            progress_bar.mod \
            stl.mod \
            vtk.mod \
            quadric.mod \
            constant.mod \
            bc.mod \
            cutcell.mod \
            parallel.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            sendrecv.mod \
            stl.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_stl_data.f 
get_connectivity.$(OBJ_EXT) : ../model/cartesian_grid/get_connectivity.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            polygon.mod \
            stl.mod \
            fldvar.mod \
            vtk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_connectivity.f 
cut_cell_preprocessing.$(OBJ_EXT) : ../model/cartesian_grid/cut_cell_preprocessing.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            vtk.mod \
            cdist.mod \
            fldvar.mod \
            polygon.mod \
            stl.mod \
            stl.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/cut_cell_preprocessing.f 
eval_usr_fct.$(OBJ_EXT) : ../model/cartesian_grid/eval_usr_fct.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            fldvar.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/eval_usr_fct.f 
allocate_cut_cell_arrays.$(OBJ_EXT) : ../model/cartesian_grid/allocate_cut_cell_arrays.f \
            param.mod \
            param1.mod \
            indices.mod \
            cutcell.mod \
            stl.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/allocate_cut_cell_arrays.f 
allocate_dummy_cut_cell_arrays.$(OBJ_EXT) : ../model/cartesian_grid/allocate_dummy_cut_cell_arrays.f \
            param.mod \
            param1.mod \
            indices.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/allocate_dummy_cut_cell_arrays.f 
calc_vort_out.$(OBJ_EXT) : ../model/cartesian_grid/calc_vort_out.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            fldvar.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/calc_vort_out.f 
deallocate_cut_cell_arrays.$(OBJ_EXT) : ../model/cartesian_grid/deallocate_cut_cell_arrays.f \
            param.mod \
            param1.mod \
            indices.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/deallocate_cut_cell_arrays.f 
get_alpha.$(OBJ_EXT) : ../model/cartesian_grid/get_alpha.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            bc.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_alpha.f 
get_cut_cell_flags.$(OBJ_EXT) : ../model/cartesian_grid/get_cut_cell_flags.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            vtk.mod \
            polygon.mod \
            stl.mod \
            physprop.mod \
            fldvar.mod \
            scalars.mod \
            funits.mod \
            rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_cut_cell_flags.f 
get_cut_cell_volume_area.$(OBJ_EXT) : ../model/cartesian_grid/get_cut_cell_volume_area.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            polygon.mod \
            stl.mod \
            bc.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_cut_cell_volume_area.f 
get_master.$(OBJ_EXT) : ../model/cartesian_grid/get_master.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            bc.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_master.f 
set_Odxyz.$(OBJ_EXT) : ../model/cartesian_grid/set_Odxyz.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            vtk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/set_Odxyz.f 
vtk_out.$(OBJ_EXT) : ../model/cartesian_grid/vtk_out.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            fldvar.mod \
            visc_s.mod \
            physprop.mod \
            pgcor.mod \
            vtk.mod \
            rxns.mod \
            output.mod \
            scalars.mod \
            mpi_utility.mod \
            parallel_mpi.mod \
            pscor.mod \
            discretelement.mod \
            mfix_pic.mod \
            cdist.mod \
            polygon.mod \
            stl.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/vtk_out.f 
write_progress_bar.$(OBJ_EXT) : ../model/cartesian_grid/write_progress_bar.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            rxns.mod \
            compar.mod \
            mpi_utility.mod \
            progress_bar.mod \
            parallel.mod \
            sendrecv.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/write_progress_bar.f 
dmp_cartesian.$(OBJ_EXT) : ../model/cartesian_grid/dmp_cartesian.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/dmp_cartesian.f 
set_geometry1.$(OBJ_EXT) : ../model/set_geometry1.f \
            param.mod \
            param1.mod \
            parallel.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_geometry1.f 
read_database.$(OBJ_EXT) : ../model/read_database.f \
            param.mod \
            param1.mod \
            physprop.mod \
            constant.mod \
            compar.mod \
            rxns.mod \
            funits.mod \
            discretelement.mod \
            des_rxns.mod \
            mfix_directory_path.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/read_database.f 
readTherm.$(OBJ_EXT) : ../model/thermochemical/readTherm.f \
            physprop.mod \
            des_rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/thermochemical/readTherm.f 
get_values.$(OBJ_EXT) : ../model/thermochemical/get_values.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/thermochemical/get_values.f 
get_bc_area.$(OBJ_EXT) : ../model/get_bc_area.f \
            param.mod \
            param1.mod \
            geometry.mod \
            bc.mod \
            compar.mod \
            parallel.mod \
            indices.mod \
            sendrecv.mod \
            mpi_utility.mod \
            cutcell.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/get_bc_area.f 
flow_to_vel.$(OBJ_EXT) : ../model/flow_to_vel.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            bc.mod \
            scales.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            discretelement.mod \
            mfix_pic.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/flow_to_vel.f 
