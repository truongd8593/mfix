.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
post_mfix : \
    ambm.mod \
    bc.mod \
    coeff.mod \
    constant.mod \
    cont.mod \
    correl.mod \
    drag.mod \
    energy.mod \
    fldvar.mod \
    funits.mod \
    geometry.mod \
    ic.mod \
    indices.mod \
    is.mod \
    leqsol.mod \
    machine.mod \
    output.mod \
    parallel.mod \
    param1.mod \
    param.mod \
    parse.mod \
    pgcor.mod \
    physprop.mod \
    post3d.mod \
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
    ur_facs.mod \
    visc_g.mod \
    visc_s.mod \
    vshear.mod \
    xsi_array.mod \
    write_error.$(OBJ_EXT)                    \
    compar.mod                                \
    dbg_util.mod                              \
    debug.mod                                 \
    gridmap.mod                               \
    mpi.mod                                   \
    mpi_utility.mod                           \
    parallel_mpi.mod                          \
    sendrecv.mod                              \
    allocate_arrays.$(OBJ_EXT) \
    any_more_data.$(OBJ_EXT) \
    boundfunijk.mod                           \
    calc_cell2.$(OBJ_EXT) \
    calc_corr_01.$(OBJ_EXT) \
    calc_corr_type_1.$(OBJ_EXT) \
    calc_distance.$(OBJ_EXT) \
    calc_ep_g.$(OBJ_EXT) \
    calc_mu_s.$(OBJ_EXT) \
    calc_mw.$(OBJ_EXT) \
    calc_quantities.$(OBJ_EXT) \
    calc_ro_g.$(OBJ_EXT) \
    calc_vol.$(OBJ_EXT) \
    check_data_03.$(OBJ_EXT) \
    check_data_04.$(OBJ_EXT) \
    check_data_05.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    compare.$(OBJ_EXT) \
    deallocate_arrays.$(OBJ_EXT) \
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
    machine.$(OBJ_EXT) \
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
    usr_init_namelist.$(OBJ_EXT) \
    usr_input.$(OBJ_EXT) \
    usr_post.$(OBJ_EXT) \
    usr_write_out1.$(OBJ_EXT) \
    write_out1.$(OBJ_EXT) \
    write_res0.$(OBJ_EXT) \
    write_res1.$(OBJ_EXT) \
    write_spx0.$(OBJ_EXT) \
    write_spx1.$(OBJ_EXT) \
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
    ornl_sym.$(OBJ_EXT) 
	$(LINK_CMD) $(LINK_FLAGS) \
    ambm_mod.$(OBJ_EXT) \
    bc_mod.$(OBJ_EXT) \
    boundfunijk_mod.$(OBJ_EXT)                           \
    coeff_mod.$(OBJ_EXT) \
    constant_mod.$(OBJ_EXT) \
    cont_mod.$(OBJ_EXT) \
    correl_mod.$(OBJ_EXT) \
    drag_mod.$(OBJ_EXT) \
    energy_mod.$(OBJ_EXT) \
    fldvar_mod.$(OBJ_EXT) \
    funits_mod.$(OBJ_EXT) \
    geometry_mod.$(OBJ_EXT) \
    ic_mod.$(OBJ_EXT) \
    indices_mod.$(OBJ_EXT) \
    is_mod.$(OBJ_EXT) \
    leqsol_mod.$(OBJ_EXT) \
    machine_mod.$(OBJ_EXT) \
    output_mod.$(OBJ_EXT) \
    parallel_mod.$(OBJ_EXT) \
    param1_mod.$(OBJ_EXT) \
    param_mod.$(OBJ_EXT) \
    parse_mod.$(OBJ_EXT) \
    pgcor_mod.$(OBJ_EXT) \
    physprop_mod.$(OBJ_EXT) \
    post3d_mod.$(OBJ_EXT) \
    pscor_mod.$(OBJ_EXT) \
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
    ur_facs_mod.$(OBJ_EXT) \
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
    calc_quantities.$(OBJ_EXT) \
    calc_ro_g.$(OBJ_EXT) \
    calc_vol.$(OBJ_EXT) \
    check_data_03.$(OBJ_EXT) \
    check_data_04.$(OBJ_EXT) \
    check_data_05.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    compare.$(OBJ_EXT) \
    deallocate_arrays.$(OBJ_EXT) \
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
    machine.$(OBJ_EXT) \
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
    usr_input.$(OBJ_EXT) \
    usr_post.$(OBJ_EXT) \
    usr_write_out1.$(OBJ_EXT) \
    write_error.$(OBJ_EXT)                               \
    write_out1.$(OBJ_EXT) \
    write_res0.$(OBJ_EXT) \
    write_res1.$(OBJ_EXT) \
    write_spx0.$(OBJ_EXT) \
    write_spx1.$(OBJ_EXT) \
    ornl_header.$(OBJ_EXT) \
    ornl_util.$(OBJ_EXT) \
    ornl_stats.$(OBJ_EXT) \
    ornl_stats_c.$(OBJ_EXT) \
    ornl_corr.$(OBJ_EXT) \
    ornl_corr_c.$(OBJ_EXT) \
    ornl_pca.$(OBJ_EXT) \
    ornl_ft.$(OBJ_EXT) \
    ornl_ft_c.$(OBJ_EXT) \
    ornl_zone.$(OBJ_EXT) \
    ornl_filt.$(OBJ_EXT) \
    ornl_filt_c.$(OBJ_EXT) \
    ornl_sym.$(OBJ_EXT) \
    debug_mod.$(OBJ_EXT)                           \
    compar_mod.$(OBJ_EXT)                          \
    dbg_util_mod.$(OBJ_EXT)                        \
    gridmap_mod.$(OBJ_EXT)                         \
    mpi_mod.$(OBJ_EXT)                             \
    mpi_utility_mod.$(OBJ_EXT)                     \
    parallel_mpi_mod.$(OBJ_EXT)                    \
    sendrecv_mod.$(OBJ_EXT)                        \
  -o post_mfix $(LIB_FLAGS)
  
ambm.mod : ../model/ambm_mod.f \
            param.mod \
            param1.mod \
            mpi_utility.mod
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/ambm_mod.f 
bc.mod : ../model/bc_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/bc_mod.f 
coeff.mod : ../model/coeff_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/coeff_mod.f 
constant.mod : ../model/constant_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/constant_mod.f 
cont.mod : ../model/cont_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/cont_mod.f 
correl.mod : correl_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) correl_mod.f 
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
leqsol.mod : ../model/leqsol_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/leqsol_mod.f 
machine.mod : ../model/machine_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/machine_mod.f 
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
            param1.mod 
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
            param1.mod 
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
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/tmp_array1_mod.f 
tmp_array.mod : ../model/tmp_array_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/tmp_array_mod.f 
toleranc.mod : ../model/toleranc_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/toleranc_mod.f 
trace.mod : ../model/trace_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/trace_mod.f 
ur_facs.mod : ../model/ur_facs_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/ur_facs_mod.f 
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
allocate_arrays.$(OBJ_EXT) : ../model/allocate_arrays.f \
            param.mod \
            param1.mod \
            ambm.mod \
            coeff.mod \
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
            scalars.mod \
            tau_g.mod \
            tau_s.mod \
            tmp_array.mod \
            tmp_array1.mod \
            trace.mod \
            visc_g.mod \
            visc_s.mod \
            xsi_array.mod \
            vshear.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/allocate_arrays.f 
any_more_data.$(OBJ_EXT) : any_more_data.f 
calc_cell2.$(OBJ_EXT) : calc_cell2.f 
calc_corr_01.$(OBJ_EXT) : calc_corr_01.f \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            correl.mod \
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
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
calc_mu_s.$(OBJ_EXT) : ../model/calc_mu_s.f \
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
            vshear.mod \
            s_pr1.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            s_pr2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/calc_mu_s.f 
calc_mw.$(OBJ_EXT) : ../model/calc_mw.f \
            param.mod \
            param1.mod \
            toleranc.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/calc_mw.f 
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
            function.inc                                                
calc_vol.$(OBJ_EXT) : calc_vol.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            fldvar.mod \
            physprop.mod \
            post3d.mod \
            function.inc                                                
check_data_03.$(OBJ_EXT) : ../model/check_data_03.f \
            param.mod \
            param1.mod \
            geometry.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/check_data_03.f 
check_data_04.$(OBJ_EXT) : ../model/check_data_04.f \
            param.mod \
            param1.mod \
            run.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/check_data_04.f 
check_data_05.$(OBJ_EXT) : ../model/check_data_05.f \
            param.mod \
            param1.mod \
            physprop.mod \
            funits.mod \
            run.mod \
            indices.mod 
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
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/compare.f 
deallocate_arrays.$(OBJ_EXT) : deallocate_arrays.f \
            param.mod \
            param1.mod \
            ambm.mod \
            coeff.mod \
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
            scalars.mod \
            tau_g.mod \
            tau_s.mod \
            tmp_array.mod \
            tmp_array1.mod \
            trace.mod \
            visc_g.mod \
            visc_s.mod \
            xsi_array.mod \
            vshear.mod 
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
            funits.mod 
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
            scalars.mod \
            xforms.inc                                                   \
            function.inc                                                
exit.$(OBJ_EXT) : ../model/exit.f 
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
            parallel_mpi.mod          \
            xforms.inc                                                  
flow_gx.$(OBJ_EXT) : flow_gx.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            function.inc                                                
flow_gy.$(OBJ_EXT) : flow_gy.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            function.inc                                                
flow_gz.$(OBJ_EXT) : flow_gz.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            function.inc                                                
flow_sx.$(OBJ_EXT) : flow_sx.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
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
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
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
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/init_namelist.f 
interp_res.$(OBJ_EXT) : interp_res.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            energy.mod \
            physprop.mod \
            fldvar.mod \
            post3d.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            xforms.inc                                                   \
            function.inc                                                
line_too_big.$(OBJ_EXT) : ../model/line_too_big.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/line_too_big.f 
machine.$(OBJ_EXT) : machine.f \
            machine.mod \
            param.mod \
            run.mod \
            funits.mod 
main_f.$(OBJ_EXT) : main_f.f \
            xforms.inc                                                  
make_upper_case.$(OBJ_EXT) : ../model/make_upper_case.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/make_upper_case.f 
open_file.$(OBJ_EXT) : ../model/open_file.f 
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
            function.inc                                                
out_time.$(OBJ_EXT) : out_time.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            physprop.mod \
            indices.mod \
            geometry.mod \
            function.inc                                                
parse_line.$(OBJ_EXT) : ../model/parse_line.f \
            param.mod \
            param1.mod \
            parse.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/parse_line.f 
parse_rxn.$(OBJ_EXT) : ../model/parse_rxn.f \
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
            usrnlst.inc                                                  \
            namelist.inc                                                
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
            scalars.mod \
            ur_facs.mod \
            toleranc.mod \
            leqsol.mod \
            tmp_array.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/read_res0.f 
read_res1.$(OBJ_EXT) : ../model/read_res1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            energy.mod \
            tmp_array.mod 
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
            scalars.mod 
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
            drag.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_constants.f 
set_dollar.$(OBJ_EXT) : set_dollar.f 
set_geometry.$(OBJ_EXT) : ../model/set_geometry.f \
            param.mod \
            param1.mod \
            run.mod \
            geometry.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_geometry.f 
set_increments.$(OBJ_EXT) : ../model/set_increments.f \
            param.mod \
            param1.mod \
            indices.mod \
            geometry.mod \
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
            boundfunijk.mod                              \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/set_index1a.f 
set_max2.$(OBJ_EXT) : ../model/set_max2.f \
            param.mod \
            param1.mod \
            geometry.mod 
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
            xforms.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
strcmp.$(OBJ_EXT) : strcmp.f 
streqs.$(OBJ_EXT) : streqs.f 
time_avg.$(OBJ_EXT) : time_avg.f \
            param.mod
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
            function.inc                                                
usr_write_out1.$(OBJ_EXT) : usr_write_out1.f 
write_out1.$(OBJ_EXT) : ../model/write_out1.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            scalars.mod \
            funits.mod 
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
            ur_facs.mod \
            leqsol.mod \
            toleranc.mod \
            tmp_array.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_res0.f 
write_res1.$(OBJ_EXT) : ../model/write_res1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            output.mod \
            energy.mod \
            tmp_array.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_res1.f 
write_spx0.$(OBJ_EXT) : ../model/write_spx0.f \
            param.mod \
            param1.mod \
            run.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_spx0.f 
write_spx1.$(OBJ_EXT) : ../model/write_spx1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            output.mod \
            scalars.mod \
            tmp_array.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_spx1.f 
usr_input.$(OBJ_EXT) : usr_input.f \
            physprop.mod  \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            compar.mod \
            constant.mod \
            post3d.mod
compar.mod : ../model/dmp_modules/mpi_donothing/compar_mod.f \
            mpi.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/compar_mod.f 
mpi.mod : ../model/dmp_modules/mpi_donothing/mpi_mod.f \
            mpif.h                                                      
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
            dbg_util.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/debug_mod.f 
gridmap.mod : ../model/dmp_modules/mpi_donothing/gridmap_mod.f \
            mpi_utility.mod \
            parallel_mpi.mod \
            geometry.mod \
            sendrecv.mod \
            compar.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/gridmap_mod.f 
mpi_utility.mod : ../model/dmp_modules/mpi_donothing/mpi_utility_mod.f \
            geometry.mod \
            compar.mod \
            parallel_mpi.mod \
            debug.mod \
            indices.mod \
            funits.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/dmp_modules/mpi_donothing/mpi_utility_mod.f 
sendrecv.mod : ../model/dmp_modules/mpi_donothing/sendrecv_mod.f \
            parallel_mpi.mod \
            debug.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            mpi.mod \
            function.inc                                                
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
write_error.$(OBJ_EXT) : ../model/write_error.f \
            param.mod \
            param1.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ../model/write_error.f 
