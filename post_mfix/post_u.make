.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
post_mfix : \
    AMBM.mod \
    BC.mod \
    CDIST.mod \
    CONSTANT.mod \
    CONT.mod \
    CORREL.mod \
    DES_BC.mod \
    DRAG.mod \
    ENERGY.mod \
    FLDVAR.mod \
    FUNITS.mod \
    GEOMETRY.mod \
    IC.mod \
    INDICES.mod \
    IS.mod \
    KINTHEORY.mod \
    LEQSOL.mod \
    MACHINE.mod \
    STIFF_CHEM.mod \
    MFLUX.mod \
    MFIX_NETCDF.mod \
    OUTPUT.mod \
    PARALLEL.mod \
    PARAM1.mod \
    PARAM.mod \
    PARSE.mod \
    PGCOR.mod \
    PHYSPROP.mod \
    POST3D.mod \
    PSCOR.mod \
    PS.mod \
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
    PROGRESS_BAR.mod \
    QUADRIC.mod \
    STL.mod \
    VTK.mod \
    CUTCELL.mod \
    POLYGON.mod \
    DASHBOARD.mod \
    COMPAR.mod \
    MPI.mod \
    DBG_UTIL.mod \
    PARALLEL_MPI.mod \
    DEBUG.mod \
    GRIDMAP.mod \
    MPI_UTILITY.mod \
    SENDRECV.mod \
    BOUNDFUNIJK.mod \
    DISCRETELEMENT.mod \
    GHDTHEORY.mod \
    POST_PRECISION.mod \
    PARALLELDATA.mod \
    USR_INPUT.mod \
    MFIX_PIC.mod \
    QMOM_KINETIC_EQUATION.mod \
    QMOMK_PARAMETERS.mod \
    DES_IC.mod \
    DES_THERMO.mod \
    DES_RXNS.mod \
    RXN_COM.mod \
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
  
AMBM.mod : ../model/ambm_mod.f \
            COMPAR.mod \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/ambm_mod.f 
BC.mod : ../model/bc_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/bc_mod.f 
CDIST.mod : ../model/cdist_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cdist_mod.f 
CONSTANT.mod : ../model/constant_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/constant_mod.f 
CONT.mod : ../model/cont_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cont_mod.f 
CORREL.mod : correl_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) correl_mod.f 
DES_BC.mod : ../model/des/des_bc_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/des_bc_mod.f 
DRAG.mod : ../model/drag_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/drag_mod.f 
ENERGY.mod : ../model/energy_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/energy_mod.f 
FLDVAR.mod : ../model/fldvar_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/fldvar_mod.f 
FUNITS.mod : ../model/funits_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/funits_mod.f 
GEOMETRY.mod : ../model/geometry_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/geometry_mod.f 
IC.mod : ../model/ic_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/ic_mod.f 
INDICES.mod : ../model/indices_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/indices_mod.f 
IS.mod : ../model/is_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/is_mod.f 
KINTHEORY.mod : ../model/kintheory_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/kintheory_mod.f 
LEQSOL.mod : ../model/leqsol_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/leqsol_mod.f 
MACHINE.mod : ../model/machine_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/machine_mod.f 
STIFF_CHEM.mod : chem/stiff_chem_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) chem/stiff_chem_mod.f 
MFLUX.mod : ../model/mflux_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/mflux_mod.f 
MFIX_NETCDF.mod : ../model/mfix_netcdf_mod.f \
            MFIX_netcdf_constants.fi                                     \
            MFIX_netcdf_overloads.fi                                     \
            MFIX_netcdf_variables.fi                                     \
            MFIX_netcdf_misc.fi                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/mfix_netcdf_mod.f 
OUTPUT.mod : ../model/output_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/output_mod.f 
PARALLEL.mod : ../model/parallel_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/parallel_mod.f 
PARAM1.mod : ../model/param1_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/param1_mod.f 
PARAM.mod : ../model/param_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/param_mod.f 
PARSE.mod : ../model/parse_mod.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod \
            COMPAR.mod \
            RXN_COM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/parse_mod.f 
PGCOR.mod : ../model/pgcor_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/pgcor_mod.f 
PHYSPROP.mod : ../model/physprop_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/physprop_mod.f 
POST3D.mod : post3d_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) post3d_mod.f 
PSCOR.mod : ../model/pscor_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/pscor_mod.f 
PS.mod : ../model/ps_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/ps_mod.f 
RESIDUAL.mod : ../model/residual_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/residual_mod.f 
RUN.mod : ../model/run_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/run_mod.f 
RXNS.mod : ../model/rxns_mod.f \
            PARAM.mod \
            PARAM1.mod \
            RXN_COM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/rxns_mod.f 
SCALARS.mod : ../model/scalars_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/scalars_mod.f 
SCALES.mod : ../model/scales_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/scales_mod.f 
TAU_G.mod : ../model/tau_g_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/tau_g_mod.f 
TAU_S.mod : ../model/tau_s_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/tau_s_mod.f 
TIME_CPU.mod : ../model/time_cpu_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/time_cpu_mod.f 
TMP_ARRAY1.mod : ../model/tmp_array1_mod.f \
            COMPAR.mod \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/tmp_array1_mod.f 
TMP_ARRAY.mod : ../model/tmp_array_mod.f \
            COMPAR.mod \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/tmp_array_mod.f 
TOLERANC.mod : ../model/toleranc_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/toleranc_mod.f 
TRACE.mod : ../model/trace_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/trace_mod.f 
TURB.mod : ../model/turb_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/turb_mod.f 
UR_FACS.mod : ../model/ur_facs_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/ur_facs_mod.f 
USR.mod : ../model/usr_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/usr_mod.f 
VISC_G.mod : ../model/visc_g_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/visc_g_mod.f 
VISC_S.mod : ../model/visc_s_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/visc_s_mod.f 
VSHEAR.mod : ../model/vshear_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/vshear_mod.f 
XSI_ARRAY.mod : ../model/xsi_array_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/xsi_array_mod.f 
PROGRESS_BAR.mod : ../model/cartesian_grid/progress_bar_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/progress_bar_mod.f 
QUADRIC.mod : ../model/cartesian_grid/quadric_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/quadric_mod.f 
STL.mod : ../model/cartesian_grid/stl_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/stl_mod.f 
VTK.mod : ../model/cartesian_grid/vtk_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/vtk_mod.f 
CUTCELL.mod : ../model/cartesian_grid/cutcell_mod.f \
            PARAM.mod \
            PARAM1.mod \
            PROGRESS_BAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/cutcell_mod.f 
POLYGON.mod : ../model/cartesian_grid/polygon_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/polygon_mod.f 
DASHBOARD.mod : ../model/cartesian_grid/dashboard_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/dashboard_mod.f 
COMPAR.mod : ../model/dmp_modules/mpi_donothing/compar_mod.f \
            MPI.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/compar_mod.f 
MPI.mod : ../model/dmp_modules/mpi_donothing/mpi_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/mpi_mod.f 
DBG_UTIL.mod : ../model/dmp_modules/mpi_donothing/dbg_util_mod.f \
            COMPAR.mod \
            GEOMETRY.mod \
            PARALLEL_MPI.mod \
            INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/dbg_util_mod.f 
PARALLEL_MPI.mod : ../model/dmp_modules/mpi_donothing/parallel_mpi_mod.f \
            GEOMETRY.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/parallel_mpi_mod.f 
DEBUG.mod : ../model/dmp_modules/mpi_donothing/debug_mod.f \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/debug_mod.f 
GRIDMAP.mod : ../model/dmp_modules/mpi_donothing/gridmap_mod.f \
            MPI_UTILITY.mod \
            PARALLEL_MPI.mod \
            GEOMETRY.mod \
            SENDRECV.mod \
            COMPAR.mod \
            RUN.mod \
            INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/gridmap_mod.f 
MPI_UTILITY.mod : ../model/dmp_modules/mpi_donothing/mpi_utility_mod.f \
            GEOMETRY.mod \
            COMPAR.mod \
            PARALLEL_MPI.mod \
            DEBUG.mod \
            INDICES.mod \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/mpi_utility_mod.f 
SENDRECV.mod : ../model/dmp_modules/mpi_donothing/sendrecv_mod.f \
            PARALLEL_MPI.mod \
            DEBUG.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/sendrecv_mod.f 
BOUNDFUNIJK.mod : ../model/boundfunijk_mod.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/boundfunijk_mod.f 
DISCRETELEMENT.mod : ../model/des/discretelement_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/discretelement_mod.f 
GHDTHEORY.mod : ../model/GhdTheory/ghdtheory_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/ghdtheory_mod.f 
POST_PRECISION.mod : post_precision_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) post_precision_mod.f 
PARALLELDATA.mod : paralleldata_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) paralleldata_mod.f 
USR_INPUT.mod : usr_input_mod.f \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_input_mod.f 
MFIX_PIC.mod : ../model/des/mfix_pic_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/mfix_pic_mod.f 
QMOM_KINETIC_EQUATION.mod : ../model/qmomk/qmom_kinetic_equation_mod.f \
            PARAM.mod \
            PARAM1.mod \
            QMOMK_PARAMETERS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/qmomk/qmom_kinetic_equation_mod.f 
QMOMK_PARAMETERS.mod : ../model/qmomk/qmomk_parameters_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/qmomk/qmomk_parameters_mod.f 
DES_IC.mod : ../model/des/des_ic_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/des_ic_mod.f 
DES_THERMO.mod : ../model/des/des_thermo_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/des_thermo_mod.f 
DES_RXNS.mod : ../model/des/des_rxns_mod.f \
            PARAM.mod \
            RXN_COM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/des_rxns_mod.f 
RXN_COM.mod : ../model/rxn_com_mod.f \
            PARAM.mod \
            PARAM1.mod \
            COMPAR.mod \
            FUNITS.mod \
            mfix_directory_path.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/rxn_com_mod.f 
machine.$(OBJ_EXT) : machine.f \
            MACHINE.mod \
            PARAM.mod \
            RUN.mod \
            FUNITS.mod 
allocate_arrays.$(OBJ_EXT) : ../model/allocate_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            AMBM.mod \
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
            RUN.mod \
            SCALARS.mod \
            TURB.mod \
            TAU_G.mod \
            TAU_S.mod \
            TMP_ARRAY.mod \
            TMP_ARRAY1.mod \
            TRACE.mod \
            VISC_G.mod \
            VISC_S.mod \
            XSI_ARRAY.mod \
            VSHEAR.mod \
            MFLUX.mod \
            GHDTHEORY.mod \
            KINTHEORY.mod \
            CDIST.mod \
            DES_RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/allocate_arrays.f 
any_more_data.$(OBJ_EXT) : any_more_data.f \
            PARAM.mod \
            PARAM1.mod 
calc_cell2.$(OBJ_EXT) : calc_cell2.f \
            PARAM.mod \
            PARAM1.mod 
calc_corr_01.$(OBJ_EXT) : calc_corr_01.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            CORREL.mod \
            COMPAR.mod \
            xforms.inc                                                   \
            function.inc                                                
calc_corr_type_1.$(OBJ_EXT) : calc_corr_type_1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            POST3D.mod \
            PHYSPROP.mod \
            CORREL.mod \
            COMPAR.mod \
            xforms.inc                                                   \
            function.inc                                                
calc_distance.$(OBJ_EXT) : calc_distance.f 
calc_ep_g.$(OBJ_EXT) : calc_ep_g.f \
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
calc_mu_s.$(OBJ_EXT) : ../model/calc_mu_s.f \
            RUN.mod \
            VSHEAR.mod \
            VISC_S.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            COMPAR.mod \
            INDICES.mod \
            GEOMETRY.mod \
            QMOM_KINETIC_EQUATION.mod \
            PARAM.mod \
            PARAM1.mod \
            TRACE.mod \
            TOLERANC.mod \
            TURB.mod \
            DRAG.mod \
            KINTHEORY.mod \
            UR_FACS.mod \
            CUTCELL.mod \
            PARALLEL.mod \
            VISC_G.mod \
            IS.mod \
            SENDRECV.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            s_pr1.inc                                                    \
            s_pr2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/calc_mu_s.f 
calc_mw.$(OBJ_EXT) : ../model/calc_mw.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/calc_mw.f 
calc_p_star.$(OBJ_EXT) : ../model/calc_p_star.f \
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
            RUN.mod \
            VISC_S.mod \
            FLDVAR.mod \
            TOLERANC.mod \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/calc_p_star.f 
calc_quantities.$(OBJ_EXT) : calc_quantities.f \
            PARAM.mod \
            PARAM1.mod \
            POST3D.mod 
calc_ro_g.$(OBJ_EXT) : calc_ro_g.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            function.inc                                                
calc_vol.$(OBJ_EXT) : calc_vol.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            POST3D.mod \
            COMPAR.mod \
            function.inc                                                
check_data_03.$(OBJ_EXT) : ../model/check_data_03.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            BC.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/check_data_03.f 
check_data_04.$(OBJ_EXT) : ../model/check_data_04.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            DISCRETELEMENT.mod \
            FUNITS.mod \
            MFIX_PIC.mod \
            COMPAR.mod \
            RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/check_data_04.f 
check_data_05.$(OBJ_EXT) : ../model/check_data_05.f \
            COMPAR.mod \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FUNITS.mod \
            RUN.mod \
            INDICES.mod \
            RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/check_data_05.f 
check_one_axis.$(OBJ_EXT) : ../model/check_one_axis.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/check_one_axis.f 
compare.$(OBJ_EXT) : ../model/compare.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/compare.f 
deallocate_arrays.$(OBJ_EXT) : deallocate_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            AMBM.mod \
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
            RUN.mod \
            SCALARS.mod \
            TAU_G.mod \
            TAU_S.mod \
            TMP_ARRAY.mod \
            TMP_ARRAY1.mod \
            TRACE.mod \
            VISC_G.mod \
            VISC_S.mod \
            XSI_ARRAY.mod \
            VSHEAR.mod \
            MFLUX.mod \
            KINTHEORY.mod \
            GHDTHEORY.mod 
cartesian_grid_init_namelist.$(OBJ_EXT) : ../model/cartesian_grid/cartesian_grid_init_namelist.f \
            PARAM1.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            POLYGON.mod \
            VTK.mod \
            PROGRESS_BAR.mod \
            DASHBOARD.mod \
            STL.mod \
            cartesian_grid/cartesian_grid_namelist.inc                  
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/cartesian_grid_init_namelist.f 
eosg.$(OBJ_EXT) : ../model/eosg.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            SCALES.mod \
            sc_p_g1.inc                                                  \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/eosg.f 
error_routine.$(OBJ_EXT) : ../model/error_routine.f \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/error_routine.f 
examine_data.$(OBJ_EXT) : examine_data.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            INDICES.mod \
            RUN.mod \
            GEOMETRY.mod \
            POST3D.mod \
            RXNS.mod \
            SCALARS.mod \
            COMPAR.mod \
            POST_PRECISION.mod \
            xforms.inc                                                   \
            function.inc                                                
exit.$(OBJ_EXT) : ../model/exit.f \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/exit.f 
f_init_data.$(OBJ_EXT) : f_init_data.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            xforms.inc                                                  
file_handle.$(OBJ_EXT) : file_handle.f \
            MACHINE.mod \
            PARAM.mod \
            PARAM1.mod \
            POST3D.mod \
            FUNITS.mod \
            COMPAR.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            function.inc                                                
finit.$(OBJ_EXT) : finit.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            POST3D.mod \
            GEOMETRY.mod \
            INDICES.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            FUNITS.mod \
            PARALLEL_MPI.mod \
            GRIDMAP.mod \
            CDIST.mod \
            PARALLELDATA.mod \
            SCALARS.mod \
            RXNS.mod \
            DRAG.mod \
            ENERGY.mod \
            xforms.inc                                                  
flow_gx.$(OBJ_EXT) : flow_gx.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            function.inc                                                
flow_gy.$(OBJ_EXT) : flow_gy.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            function.inc                                                
flow_gz.$(OBJ_EXT) : flow_gz.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            function.inc                                                
flow_sx.$(OBJ_EXT) : flow_sx.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
flow_sy.$(OBJ_EXT) : flow_sy.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
flow_sz.$(OBJ_EXT) : flow_sz.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
g_0.$(OBJ_EXT) : ../model/g_0.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            VISC_S.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/g_0.f 
gas_flux.$(OBJ_EXT) : gas_flux.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            POST3D.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            xforms.inc                                                   \
            function.inc                                                
get_file_name.$(OBJ_EXT) : get_file_name.f 
get_file_status.$(OBJ_EXT) : get_file_status.f 
get_index.$(OBJ_EXT) : get_index.f 
get_location.$(OBJ_EXT) : get_location.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            POST3D.mod 
get_mu_s.$(OBJ_EXT) : get_mu_s.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            VISC_S.mod \
            RUN.mod \
            CONSTANT.mod \
            FUNITS.mod \
            COMPAR.mod \
            xforms.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
get_same_time.$(OBJ_EXT) : get_same_time.f \
            PARAM.mod \
            PARAM1.mod 
get_selection.$(OBJ_EXT) : get_selection.f 
get_substr.$(OBJ_EXT) : get_substr.f 
granular_qty.$(OBJ_EXT) : granular_qty.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            VISC_S.mod \
            RUN.mod \
            CONSTANT.mod \
            FUNITS.mod \
            COMPAR.mod \
            xforms.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
header_main.$(OBJ_EXT) : header_main.f \
            PARAM.mod \
            PARAM1.mod \
            POST3D.mod 
ik_avg.$(OBJ_EXT) : ik_avg.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            function.inc                                                
ik_avg_out.$(OBJ_EXT) : ik_avg_out.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            INDICES.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            POST3D.mod 
init_namelist.$(OBJ_EXT) : ../model/init_namelist.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            IC.mod \
            BC.mod \
            PS.mod \
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
            PARALLEL.mod \
            CDIST.mod \
            STIFF_CHEM.mod \
            namelist.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/init_namelist.f 
interp_res.$(OBJ_EXT) : interp_res.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            ENERGY.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            POST3D.mod \
            RUN.mod \
            RXNS.mod \
            SCALARS.mod \
            FUNITS.mod \
            COMPAR.mod \
            GRIDMAP.mod \
            TMP_ARRAY.mod \
            xforms.inc                                                   \
            function.inc                                                
kintheory_energy_dissipation_ss.$(OBJ_EXT) : ../model/kintheory_energy_dissipation_ss.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            KINTHEORY.mod \
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
            CDIST.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/open_file.f 
out_from_res.$(OBJ_EXT) : out_from_res.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            FUNITS.mod \
            xforms.inc                                                  
out_from_spx.$(OBJ_EXT) : out_from_spx.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            FUNITS.mod \
            POST3D.mod \
            xforms.inc                                                  
out_spec_time.$(OBJ_EXT) : out_spec_time.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            PHYSPROP.mod \
            INDICES.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            function.inc                                                
out_time.$(OBJ_EXT) : out_time.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            PHYSPROP.mod \
            INDICES.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            function.inc                                                
parse_line.$(OBJ_EXT) : ../model/parse_line.f \
            COMPAR.mod \
            DES_RXNS.mod \
            PARAM.mod \
            PARAM1.mod \
            PARSE.mod \
            RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/parse_line.f 
parse_rxn.$(OBJ_EXT) : ../model/parse_rxn.f \
            COMPAR.mod \
            FUNITS.mod \
            PARAM.mod \
            PARAM1.mod \
            PARSE.mod \
            RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/parse_rxn.f 
print_out.$(OBJ_EXT) : print_out.f \
            PARAM.mod \
            PARAM1.mod \
            POST3D.mod \
            xforms.inc                                                  
read_namelist.$(OBJ_EXT) : ../model/read_namelist.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            IC.mod \
            IS.mod \
            BC.mod \
            PS.mod \
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
            PARALLEL.mod \
            DISCRETELEMENT.mod \
            MFIX_PIC.mod \
            USR.mod \
            DES_BC.mod \
            DES_IC.mod \
            DES_THERMO.mod \
            DES_RXNS.mod \
            STIFF_CHEM.mod \
            CDIST.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            VTK.mod \
            POLYGON.mod \
            DASHBOARD.mod \
            STL.mod \
            QMOM_KINETIC_EQUATION.mod \
            usrnlst.inc                                                  \
            namelist.inc                                                 \
            des/desnamelist.inc                                          \
            cartesian_grid/cartesian_grid_namelist.inc                   \
            qmomk/qmomknamelist.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/read_namelist.f 
read_res0.$(OBJ_EXT) : ../model/read_res0.f \
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
            RXNS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            FLDVAR.mod \
            STIFF_CHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/read_res0.f 
read_res1.$(OBJ_EXT) : ../model/read_res1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            RXNS.mod \
            SCALARS.mod \
            FUNITS.mod \
            ENERGY.mod \
            COMPAR.mod \
            CDIST.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            MFIX_NETCDF.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/read_res1.f 
read_spx0.$(OBJ_EXT) : read_spx0.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            FUNITS.mod \
            POST3D.mod 
read_spx1.$(OBJ_EXT) : read_spx1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            FUNITS.mod \
            POST3D.mod \
            SCALARS.mod \
            RXNS.mod \
            MACHINE.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
remove_comment.$(OBJ_EXT) : ../model/remove_comment.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/remove_comment.f 
res_from_spx.$(OBJ_EXT) : res_from_spx.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            FUNITS.mod \
            POST3D.mod \
            PHYSPROP.mod \
            FLDVAR.mod 
seek_comment.$(OBJ_EXT) : ../model/seek_comment.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/seek_comment.f 
seek_end.$(OBJ_EXT) : ../model/seek_end.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/seek_end.f 
seek_time.$(OBJ_EXT) : seek_time.f \
            PARAM.mod \
            PARAM1.mod \
            POST3D.mod \
            FUNITS.mod 
select_spx_rec.$(OBJ_EXT) : select_spx_rec.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            MACHINE.mod \
            FUNITS.mod \
            POST3D.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            xforms.inc                                                  
set_constants.$(OBJ_EXT) : ../model/set_constants.f \
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
            DRAG.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_constants.f 
set_dollar.$(OBJ_EXT) : set_dollar.f 
set_geometry.$(OBJ_EXT) : ../model/set_geometry.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            GEOMETRY.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_geometry.f 
set_increments.$(OBJ_EXT) : ../model/set_increments.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            FUNITS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_increments.f 
set_index1.$(OBJ_EXT) : ../model/set_index1.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            CONSTANT.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_index1.f 
set_index1a.$(OBJ_EXT) : ../model/set_index1a.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            BOUNDFUNIJK.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_index1a.f 
set_max2.$(OBJ_EXT) : ../model/set_max2.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_max2.f 
set_read_spx.$(OBJ_EXT) : set_read_spx.f 
shift_dxyz.$(OBJ_EXT) : ../model/shift_dxyz.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/shift_dxyz.f 
sol_flux.$(OBJ_EXT) : sol_flux.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            POST3D.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            xforms.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
strcmp.$(OBJ_EXT) : strcmp.f 
streqs.$(OBJ_EXT) : streqs.f 
time_avg.$(OBJ_EXT) : time_avg.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            MACHINE.mod \
            FUNITS.mod \
            POST3D.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            SCALARS.mod \
            RXNS.mod \
            xforms.inc                                                  
usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f \
            usrnlst.inc                                                 
usr_post.$(OBJ_EXT) : usr_post.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            POST3D.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            function.inc                                                
usr_write_out1.$(OBJ_EXT) : usr_write_out1.f 
write_out1.$(OBJ_EXT) : ../model/write_out1.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            RXNS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_out1.f 
write_res0.$(OBJ_EXT) : ../model/write_res0.f \
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
            RXNS.mod \
            UR_FACS.mod \
            LEQSOL.mod \
            TOLERANC.mod \
            CDIST.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            STIFF_CHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_res0.f 
write_res1.$(OBJ_EXT) : ../model/write_res1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            SCALARS.mod \
            RXNS.mod \
            FUNITS.mod \
            OUTPUT.mod \
            ENERGY.mod \
            CDIST.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            MFIX_NETCDF.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_res1.f 
write_spx0.$(OBJ_EXT) : ../model/write_spx0.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            FUNITS.mod \
            CDIST.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_spx0.f 
write_spx1.$(OBJ_EXT) : ../model/write_spx1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            FUNITS.mod \
            SCALARS.mod \
            OUTPUT.mod \
            RXNS.mod \
            CDIST.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            MFIX_NETCDF.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_spx1.f 
usr_input.$(OBJ_EXT) : usr_input.f \
            USR_INPUT.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            COMPAR.mod \
            CONSTANT.mod \
            POST3D.mod \
            xforms.inc                                                   \
            function.inc                                                
write_error.$(OBJ_EXT) : ../model/write_error.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_error.f 
transport_coeff_ghd.$(OBJ_EXT) : ../model/GhdTheory/transport_coeff_ghd.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            VISC_S.mod \
            GHDTHEORY.mod \
            PHYSPROP.mod \
            RUN.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/transport_coeff_ghd.f 
ghd.$(OBJ_EXT) : ../model/GhdTheory/ghd.f \
            DRAG.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/GhdTheory/ghd.f 
cooling_rate.$(OBJ_EXT) : ../model/GhdTheory/cooling_rate.f \
            COMPAR.mod 
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
            PARAM1.mod \
            DISCRETELEMENT.mod \
            MFIX_PIC.mod \
            DES_BC.mod \
            DES_IC.mod \
            DES_THERMO.mod \
            DES_RXNS.mod \
            des/desnamelist.inc                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/des/des_init_namelist.f 
ornl_header.$(OBJ_EXT) : ornl_header.f \
            GEOMETRY.mod \
            RUN.mod \
            USR_INPUT.mod 
ornl_util.$(OBJ_EXT) : ornl_util.f 
ornl_stats.$(OBJ_EXT) : ornl_stats.f 
ornl_stats_c.$(OBJ_EXT) : ornl_stats_c.f \
            GEOMETRY.mod \
            FLDVAR.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            USR_INPUT.mod \
            function.inc                                                
ornl_corr.$(OBJ_EXT) : ornl_corr.f 
ornl_corr_c.$(OBJ_EXT) : ornl_corr_c.f \
            GEOMETRY.mod \
            FLDVAR.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            USR_INPUT.mod \
            function.inc                                                
ornl_pca.$(OBJ_EXT) : ornl_pca.f 
ornl_ft.$(OBJ_EXT) : ornl_ft.f 
ornl_ft_c.$(OBJ_EXT) : ornl_ft_c.f \
            USR_INPUT.mod 
ornl_filt.$(OBJ_EXT) : ornl_filt.f 
ornl_filt_c.$(OBJ_EXT) : ornl_filt_c.f \
            USR_INPUT.mod 
ornl_zone.$(OBJ_EXT) : ornl_zone.f \
            GEOMETRY.mod \
            FLDVAR.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            USR_INPUT.mod \
            function.inc                                                
ornl_sym.$(OBJ_EXT) : ornl_sym.f 
qmomk_init_namelist.$(OBJ_EXT) : ../model/qmomk/qmomk_init_namelist.f \
            PARAM1.mod \
            QMOM_KINETIC_EQUATION.mod \
            qmomk/qmomknamelist.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/qmomk/qmomk_init_namelist.f 
calc_trd_s.$(OBJ_EXT) : ../model/calc_trd_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            SENDRECV.mod \
            BC.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/calc_trd_s.f 
get_delh.$(OBJ_EXT) : ../model/cartesian_grid/get_delh.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_delh.f 
define_quadrics.$(OBJ_EXT) : ../model/cartesian_grid/define_quadrics.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            VTK.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/define_quadrics.f 
check_data_cartesian.$(OBJ_EXT) : ../model/cartesian_grid/check_data_cartesian.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            RUN.mod \
            PHYSPROP.mod \
            INDICES.mod \
            SCALARS.mod \
            FUNITS.mod \
            LEQSOL.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            BC.mod \
            DISCRETELEMENT.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            VTK.mod \
            POLYGON.mod \
            DASHBOARD.mod \
            STL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/check_data_cartesian.f 
get_poly_data.$(OBJ_EXT) : ../model/cartesian_grid/get_poly_data.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            RXNS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PROGRESS_BAR.mod \
            POLYGON.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_poly_data.f 
get_stl_data.$(OBJ_EXT) : ../model/cartesian_grid/get_stl_data.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            RXNS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PROGRESS_BAR.mod \
            STL.mod \
            VTK.mod \
            QUADRIC.mod \
            CONSTANT.mod \
            BC.mod \
            CUTCELL.mod \
            PARALLEL.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            SENDRECV.mod \
            STL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_stl_data.f 
get_connectivity.$(OBJ_EXT) : ../model/cartesian_grid/get_connectivity.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            POLYGON.mod \
            STL.mod \
            FLDVAR.mod \
            VTK.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_connectivity.f 
cut_cell_preprocessing.$(OBJ_EXT) : ../model/cartesian_grid/cut_cell_preprocessing.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            VTK.mod \
            CDIST.mod \
            FLDVAR.mod \
            POLYGON.mod \
            STL.mod \
            STL.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/cut_cell_preprocessing.f 
eval_usr_fct.$(OBJ_EXT) : ../model/cartesian_grid/eval_usr_fct.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            FLDVAR.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/eval_usr_fct.f 
allocate_cut_cell_arrays.$(OBJ_EXT) : ../model/cartesian_grid/allocate_cut_cell_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            CUTCELL.mod \
            STL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/allocate_cut_cell_arrays.f 
allocate_dummy_cut_cell_arrays.$(OBJ_EXT) : ../model/cartesian_grid/allocate_dummy_cut_cell_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/allocate_dummy_cut_cell_arrays.f 
calc_vort_out.$(OBJ_EXT) : ../model/cartesian_grid/calc_vort_out.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            FLDVAR.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/calc_vort_out.f 
deallocate_cut_cell_arrays.$(OBJ_EXT) : ../model/cartesian_grid/deallocate_cut_cell_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/deallocate_cut_cell_arrays.f 
get_alpha.$(OBJ_EXT) : ../model/cartesian_grid/get_alpha.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            BC.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_alpha.f 
get_cut_cell_flags.$(OBJ_EXT) : ../model/cartesian_grid/get_cut_cell_flags.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            VTK.mod \
            POLYGON.mod \
            STL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            SCALARS.mod \
            FUNITS.mod \
            RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_cut_cell_flags.f 
get_cut_cell_volume_area.$(OBJ_EXT) : ../model/cartesian_grid/get_cut_cell_volume_area.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            POLYGON.mod \
            STL.mod \
            BC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_cut_cell_volume_area.f 
get_master.$(OBJ_EXT) : ../model/cartesian_grid/get_master.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            BC.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/get_master.f 
set_Odxyz.$(OBJ_EXT) : ../model/cartesian_grid/set_Odxyz.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            VTK.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/set_Odxyz.f 
vtk_out.$(OBJ_EXT) : ../model/cartesian_grid/vtk_out.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            FLDVAR.mod \
            VISC_S.mod \
            PHYSPROP.mod \
            PGCOR.mod \
            VTK.mod \
            RXNS.mod \
            OUTPUT.mod \
            SCALARS.mod \
            MPI_UTILITY.mod \
            PARALLEL_MPI.mod \
            PSCOR.mod \
            DISCRETELEMENT.mod \
            MFIX_PIC.mod \
            CDIST.mod \
            POLYGON.mod \
            STL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/vtk_out.f 
write_progress_bar.$(OBJ_EXT) : ../model/cartesian_grid/write_progress_bar.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            RXNS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PROGRESS_BAR.mod \
            PARALLEL.mod \
            SENDRECV.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/write_progress_bar.f 
dmp_cartesian.$(OBJ_EXT) : ../model/cartesian_grid/dmp_cartesian.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cartesian_grid/dmp_cartesian.f 
set_geometry1.$(OBJ_EXT) : ../model/set_geometry1.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_geometry1.f 
read_database.$(OBJ_EXT) : ../model/read_database.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            COMPAR.mod \
            RXNS.mod \
            FUNITS.mod \
            DISCRETELEMENT.mod \
            DES_RXNS.mod \
            mfix_directory_path.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/read_database.f 
readTherm.$(OBJ_EXT) : ../model/thermochemical/readTherm.f \
            PHYSPROP.mod \
            DES_RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/thermochemical/readTherm.f 
get_values.$(OBJ_EXT) : ../model/thermochemical/get_values.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/thermochemical/get_values.f 
get_bc_area.$(OBJ_EXT) : ../model/get_bc_area.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            BC.mod \
            COMPAR.mod \
            PARALLEL.mod \
            INDICES.mod \
            SENDRECV.mod \
            MPI_UTILITY.mod \
            CUTCELL.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/get_bc_area.f 
flow_to_vel.$(OBJ_EXT) : ../model/flow_to_vel.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            BC.mod \
            SCALES.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod \
            MFIX_PIC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/flow_to_vel.f 
