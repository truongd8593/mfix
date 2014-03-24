.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
$(EXEC_FILE) : \
    $(DPO)ambm.mod \
    $(DPO)bc.mod \
    $(DPO)boundfunijk3.mod \
    $(DPO)boundfunijk.mod \
    $(DPO)cdist.mod \
    $(DPO)check.mod \
    $(DPO)chischeme.mod \
    $(DPO)coeff.mod \
    $(DPO)constant.mod \
    $(DPO)cont.mod \
    $(DPO)corner.mod \
    $(DPO)dbg.mod \
    $(DPO)drag.mod \
    $(DPO)energy.mod \
    $(DPO)error_manager.mod \
    $(DPO)fldvar.mod \
    $(DPO)function.mod \
    $(DPO)funits.mod \
    $(DPO)geometry.mod \
    $(DPO)ic.mod \
    $(DPO)indices.mod \
    $(DPO)is.mod \
    $(DPO)kintheory.mod \
    $(DPO)leqsol.mod \
    $(DPO)machine.mod \
    $(DPO)matrix.mod \
    $(DPO)mfix_netcdf.mod \
    $(DPO)mflux.mod \
    $(DPO)mms.mod \
    $(DPO)output.mod \
    $(DPO)parallel.mod \
    $(DPO)param1.mod \
    $(DPO)param.mod \
    $(DPO)parse.mod \
    $(DPO)pgcor.mod \
    $(DPO)physprop.mod \
    $(DPO)pscor.mod \
    $(DPO)ps.mod \
    $(DPO)residual.mod \
    $(DPO)run.mod \
    $(DPO)rxn_com.mod \
    $(DPO)rxns.mod \
    $(DPO)scalars.mod \
    $(DPO)scales.mod \
    $(DPO)tau_g.mod \
    $(DPO)tau_s.mod \
    $(DPO)time_cpu.mod \
    $(DPO)tmp_array1.mod \
    $(DPO)tmp_array.mod \
    $(DPO)toleranc.mod \
    $(DPO)trace.mod \
    $(DPO)turb.mod \
    $(DPO)ur_facs.mod \
    $(DPO)usr.mod \
    $(DPO)visc_g.mod \
    $(DPO)visc_s.mod \
    $(DPO)vshear.mod \
    $(DPO)xsi_array.mod \
    $(DPO)cutcell.mod \
    $(DPO)dashboard.mod \
    $(DPO)polygon.mod \
    $(DPO)progress_bar.mod \
    $(DPO)quadric.mod \
    $(DPO)stl.mod \
    $(DPO)vtk.mod \
    $(DPO)stiff_chem_dbg.mod \
    $(DPO)stiff_chem_maps.mod \
    $(DPO)stiff_chem.mod \
    $(DPO)stiff_chem_stats.mod \
    $(DPO)des_bc.mod \
    $(DPO)des_cluster.mod \
    $(DPO)desgrid.mod \
    $(DPO)des_ic.mod \
    $(DPO)des_linked_list_data.mod \
    $(DPO)des_linked_list_funcs.mod \
    $(DPO)desmpi.mod \
    $(DPO)desmpi_wrapper.mod \
    $(DPO)des_rxns.mod \
    $(DPO)des_stl_functions.mod \
    $(DPO)des_thermo.mod \
    $(DPO)discretelement.mod \
    $(DPO)interpolation.mod \
    $(DPO)mfix_pic.mod \
    $(DPO)mppic_wallbc.mod \
    $(DPO)randomno.mod \
    $(DPO)sendrecvnode.mod \
    $(DPO)softspring_funcs_cutcell.mod \
    $(DPO)compar.mod \
    $(DPO)dbg_util.mod \
    $(DPO)debug.mod \
    $(DPO)gridmap.mod \
    $(DPO)mpi.mod \
    $(DPO)mpi_utility.mod \
    $(DPO)parallel_mpi.mod \
    $(DPO)sendrecv3.mod \
    $(DPO)sendrecv.mod \
    $(DPO)ghdtheory.mod \
    $(DPO)qmomk_bc.mod \
    $(DPO)qmomk_collision.mod \
    $(DPO)qmomk_fluxes.mod \
    $(DPO)qmom_kinetic_equation.mod \
    $(DPO)qmomk_parameters.mod \
    $(DPO)qmomk_quadrature.mod \
    $(DPO)qmomk_tools.mod \
    $(DPO)accum_resid.$(OBJ_EXT) \
    $(DPO)adjust_a_u_g.$(OBJ_EXT) \
    $(DPO)adjust_a_u_s.$(OBJ_EXT) \
    $(DPO)adjust_a_v_g.$(OBJ_EXT) \
    $(DPO)adjust_a_v_s.$(OBJ_EXT) \
    $(DPO)adjust_a_w_g.$(OBJ_EXT) \
    $(DPO)adjust_a_w_s.$(OBJ_EXT) \
    $(DPO)adjust_dt.$(OBJ_EXT) \
    $(DPO)adjust_eps.$(OBJ_EXT) \
    $(DPO)adjust_leq.$(OBJ_EXT) \
    $(DPO)adjust_rop.$(OBJ_EXT) \
    $(DPO)adjust_theta.$(OBJ_EXT) \
    $(DPO)allocate_arrays.$(OBJ_EXT) \
    $(DPO)bc_phi.$(OBJ_EXT) \
    $(DPO)bc_theta.$(OBJ_EXT) \
    $(DPO)b_m_p_star.$(OBJ_EXT) \
    $(DPO)bound_x.$(OBJ_EXT) \
    $(DPO)calc_cell.$(OBJ_EXT) \
    $(DPO)calc_coeff.$(OBJ_EXT) \
    $(DPO)calc_d.$(OBJ_EXT) \
    $(DPO)calc_dif_g.$(OBJ_EXT) \
    $(DPO)calc_dif_s.$(OBJ_EXT) \
    $(DPO)calc_drag.$(OBJ_EXT) \
    $(DPO)calc_e.$(OBJ_EXT) \
    $(DPO)calc_gama.$(OBJ_EXT) \
    $(DPO)calc_grbdry.$(OBJ_EXT) \
    $(DPO)calc_h.$(OBJ_EXT) \
    $(DPO)calc_k_cp.$(OBJ_EXT) \
    $(DPO)calc_k_g.$(OBJ_EXT) \
    $(DPO)calc_k_s.$(OBJ_EXT) \
    $(DPO)calc_mflux.$(OBJ_EXT) \
    $(DPO)calc_mu_g.$(OBJ_EXT) \
    $(DPO)calc_mu_s.$(OBJ_EXT) \
    $(DPO)calc_mw.$(OBJ_EXT) \
    $(DPO)calc_outflow.$(OBJ_EXT) \
    $(DPO)calc_p_star.$(OBJ_EXT) \
    $(DPO)calc_resid.$(OBJ_EXT) \
    $(DPO)calc_s_ddot_s.$(OBJ_EXT) \
    $(DPO)calc_trd_g.$(OBJ_EXT) \
    $(DPO)calc_trd_s.$(OBJ_EXT) \
    $(DPO)calc_u_friction.$(OBJ_EXT) \
    $(DPO)calc_vol_fr.$(OBJ_EXT) \
    $(DPO)calc_xsi.$(OBJ_EXT) \
    $(DPO)cal_d.$(OBJ_EXT) \
    $(DPO)check_ab_m.$(OBJ_EXT) \
    $(DPO)check_convergence.$(OBJ_EXT) \
    $(DPO)check_data_03.$(OBJ_EXT) \
    $(DPO)check_data_06.$(OBJ_EXT) \
    $(DPO)check_data_07.$(OBJ_EXT) \
    $(DPO)check_data_08.$(OBJ_EXT) \
    $(DPO)check_data_09.$(OBJ_EXT) \
    $(DPO)check_data_10.$(OBJ_EXT) \
    $(DPO)check_data_20.$(OBJ_EXT) \
    $(DPO)check_data_30.$(OBJ_EXT) \
    $(DPO)check_mass_balance.$(OBJ_EXT) \
    $(DPO)check_one_axis.$(OBJ_EXT) \
    $(DPO)check_plane.$(OBJ_EXT) \
    $(DPO)cn_extrapol.$(OBJ_EXT) \
    $(DPO)compare.$(OBJ_EXT) \
    $(DPO)conv_dif_phi.$(OBJ_EXT) \
    $(DPO)conv_dif_u_g.$(OBJ_EXT) \
    $(DPO)conv_dif_u_s.$(OBJ_EXT) \
    $(DPO)conv_dif_v_g.$(OBJ_EXT) \
    $(DPO)conv_dif_v_s.$(OBJ_EXT) \
    $(DPO)conv_dif_w_g.$(OBJ_EXT) \
    $(DPO)conv_dif_w_s.$(OBJ_EXT) \
    $(DPO)conv_pp_g.$(OBJ_EXT) \
    $(DPO)conv_rop.$(OBJ_EXT) \
    $(DPO)conv_rop_g.$(OBJ_EXT) \
    $(DPO)conv_rop_s.$(OBJ_EXT) \
    $(DPO)conv_source_epp.$(OBJ_EXT) \
    $(DPO)copy_a.$(OBJ_EXT) \
    $(DPO)corner.$(OBJ_EXT) \
    $(DPO)correct_0.$(OBJ_EXT) \
    $(DPO)correct_1.$(OBJ_EXT) \
    $(DPO)dgtsl.$(OBJ_EXT) \
    $(DPO)dif_u_is.$(OBJ_EXT) \
    $(DPO)dif_v_is.$(OBJ_EXT) \
    $(DPO)dif_w_is.$(OBJ_EXT) \
    $(DPO)discretize.$(OBJ_EXT) \
    $(DPO)display_resid.$(OBJ_EXT) \
    $(DPO)drag_gs.$(OBJ_EXT) \
    $(DPO)drag_ss.$(OBJ_EXT) \
    $(DPO)eosg.$(OBJ_EXT) \
    $(DPO)eoss.$(OBJ_EXT) \
    $(DPO)equal.$(OBJ_EXT) \
    $(DPO)error_routine.$(OBJ_EXT) \
    $(DPO)exchange.$(OBJ_EXT) \
    $(DPO)exit.$(OBJ_EXT) \
    $(DPO)flow_to_vel.$(OBJ_EXT) \
    $(DPO)g_0.$(OBJ_EXT) \
    $(DPO)get_bc_area.$(OBJ_EXT) \
    $(DPO)get_data.$(OBJ_EXT) \
    $(DPO)get_eq.$(OBJ_EXT) \
    $(DPO)get_flow_bc.$(OBJ_EXT) \
    $(DPO)get_hloss.$(OBJ_EXT) \
    $(DPO)get_is.$(OBJ_EXT) \
    $(DPO)get_philoss.$(OBJ_EXT) \
    $(DPO)get_ps.$(OBJ_EXT) \
    $(DPO)get_smass.$(OBJ_EXT) \
    $(DPO)get_stats.$(OBJ_EXT) \
    $(DPO)get_walls_bc.$(OBJ_EXT) \
    $(DPO)in_bin_512.$(OBJ_EXT) \
    $(DPO)in_bin_512i.$(OBJ_EXT) \
    $(DPO)init_ab_m.$(OBJ_EXT) \
    $(DPO)init_fvars.$(OBJ_EXT) \
    $(DPO)init_namelist.$(OBJ_EXT) \
    $(DPO)init_resid.$(OBJ_EXT) \
    $(DPO)iterate.$(OBJ_EXT) \
    $(DPO)k_epsilon_prop.$(OBJ_EXT) \
    $(DPO)kintheory_drag_ss.$(OBJ_EXT) \
    $(DPO)kintheory_energy_dissipation_ss.$(OBJ_EXT) \
    $(DPO)kintheory_u_s.$(OBJ_EXT) \
    $(DPO)kintheory_v_s.$(OBJ_EXT) \
    $(DPO)kintheory_w_s.$(OBJ_EXT) \
    $(DPO)leq_bicgs.$(OBJ_EXT) \
    $(DPO)leq_bicgst.$(OBJ_EXT) \
    $(DPO)leq_cg.$(OBJ_EXT) \
    $(DPO)leq_gmres.$(OBJ_EXT) \
    $(DPO)leq_sor.$(OBJ_EXT) \
    $(DPO)line_too_big.$(OBJ_EXT) \
    $(DPO)location_check.$(OBJ_EXT) \
    $(DPO)location.$(OBJ_EXT) \
    $(DPO)machine.$(OBJ_EXT) \
    $(DPO)make_upper_case.$(OBJ_EXT) \
    $(DPO)mark_phase_4_cor.$(OBJ_EXT) \
    $(DPO)mfix.$(OBJ_EXT) \
    $(DPO)mod_bc_i.$(OBJ_EXT) \
    $(DPO)mod_bc_j.$(OBJ_EXT) \
    $(DPO)mod_bc_k.$(OBJ_EXT) \
    $(DPO)open_file.$(OBJ_EXT) \
    $(DPO)open_files.$(OBJ_EXT) \
    $(DPO)out_array_c.$(OBJ_EXT) \
    $(DPO)out_array.$(OBJ_EXT) \
    $(DPO)out_array_kc.$(OBJ_EXT) \
    $(DPO)out_array_k.$(OBJ_EXT) \
    $(DPO)out_bin_512.$(OBJ_EXT) \
    $(DPO)out_bin_512i.$(OBJ_EXT) \
    $(DPO)out_bin_512r.$(OBJ_EXT) \
    $(DPO)out_bin_r.$(OBJ_EXT) \
    $(DPO)parse_line.$(OBJ_EXT) \
    $(DPO)parse_resid_string.$(OBJ_EXT) \
    $(DPO)parse_rxn.$(OBJ_EXT) \
    $(DPO)partial_elim.$(OBJ_EXT) \
    $(DPO)physical_prop.$(OBJ_EXT) \
    $(DPO)read_database.$(OBJ_EXT) \
    $(DPO)read_namelist.$(OBJ_EXT) \
    $(DPO)read_res0.$(OBJ_EXT) \
    $(DPO)read_res1.$(OBJ_EXT) \
    $(DPO)remove_comment.$(OBJ_EXT) \
    $(DPO)reset_new.$(OBJ_EXT) \
    $(DPO)rrates0.$(OBJ_EXT) \
    $(DPO)rrates.$(OBJ_EXT) \
    $(DPO)rrates_init.$(OBJ_EXT) \
    $(DPO)scalar_prop.$(OBJ_EXT) \
    $(DPO)seek_comment.$(OBJ_EXT) \
    $(DPO)seek_end.$(OBJ_EXT) \
    $(DPO)set_bc0.$(OBJ_EXT) \
    $(DPO)set_bc1.$(OBJ_EXT) \
    $(DPO)set_bc_flow.$(OBJ_EXT) \
    $(DPO)set_constants.$(OBJ_EXT) \
    $(DPO)set_constprop.$(OBJ_EXT) \
    $(DPO)set_flags.$(OBJ_EXT) \
    $(DPO)set_fluidbed_p.$(OBJ_EXT) \
    $(DPO)set_geometry1.$(OBJ_EXT) \
    $(DPO)set_geometry.$(OBJ_EXT) \
    $(DPO)set_icbc_flags.$(OBJ_EXT) \
    $(DPO)set_ic.$(OBJ_EXT) \
    $(DPO)set_increments3.$(OBJ_EXT) \
    $(DPO)set_increments.$(OBJ_EXT) \
    $(DPO)set_index1a3.$(OBJ_EXT) \
    $(DPO)set_index1a.$(OBJ_EXT) \
    $(DPO)set_index1.$(OBJ_EXT) \
    $(DPO)set_l_scale.$(OBJ_EXT) \
    $(DPO)set_max2.$(OBJ_EXT) \
    $(DPO)set_mw_mix_g.$(OBJ_EXT) \
    $(DPO)set_outflow.$(OBJ_EXT) \
    $(DPO)set_ps.$(OBJ_EXT) \
    $(DPO)set_ro_g.$(OBJ_EXT) \
    $(DPO)set_ro_s.$(OBJ_EXT) \
    $(DPO)set_wall_bc.$(OBJ_EXT) \
    $(DPO)shift_dxyz.$(OBJ_EXT) \
    $(DPO)solve_continuity.$(OBJ_EXT) \
    $(DPO)solve_energy_eq.$(OBJ_EXT) \
    $(DPO)solve_epp.$(OBJ_EXT) \
    $(DPO)solve_granular_energy.$(OBJ_EXT) \
    $(DPO)solve_k_epsilon_eq.$(OBJ_EXT) \
    $(DPO)solve_lin_eq.$(OBJ_EXT) \
    $(DPO)solve_pp_g.$(OBJ_EXT) \
    $(DPO)solve_scalar_eq.$(OBJ_EXT) \
    $(DPO)solve_species_eq.$(OBJ_EXT) \
    $(DPO)solve_vel_star.$(OBJ_EXT) \
    $(DPO)source_granular_energy.$(OBJ_EXT) \
    $(DPO)source_phi.$(OBJ_EXT) \
    $(DPO)source_pp_g.$(OBJ_EXT) \
    $(DPO)source_rop_g.$(OBJ_EXT) \
    $(DPO)source_rop_s.$(OBJ_EXT) \
    $(DPO)source_u_g.$(OBJ_EXT) \
    $(DPO)source_u_s.$(OBJ_EXT) \
    $(DPO)source_v_g.$(OBJ_EXT) \
    $(DPO)source_v_s.$(OBJ_EXT) \
    $(DPO)source_w_g.$(OBJ_EXT) \
    $(DPO)source_w_s.$(OBJ_EXT) \
    $(DPO)tau_u_g.$(OBJ_EXT) \
    $(DPO)tau_u_s.$(OBJ_EXT) \
    $(DPO)tau_v_g.$(OBJ_EXT) \
    $(DPO)tau_v_s.$(OBJ_EXT) \
    $(DPO)tau_w_g.$(OBJ_EXT) \
    $(DPO)tau_w_s.$(OBJ_EXT) \
    $(DPO)test_lin_eq.$(OBJ_EXT) \
    $(DPO)time_march.$(OBJ_EXT) \
    $(DPO)transfer.$(OBJ_EXT) \
    $(DPO)transport_prop.$(OBJ_EXT) \
    $(DPO)undef_2_0.$(OBJ_EXT) \
    $(DPO)under_relax.$(OBJ_EXT) \
    $(DPO)update_old.$(OBJ_EXT) \
    $(DPO)usr0.$(OBJ_EXT) \
    $(DPO)usr1.$(OBJ_EXT) \
    $(DPO)usr2.$(OBJ_EXT) \
    $(DPO)usr3.$(OBJ_EXT) \
    $(DPO)usr_init_namelist.$(OBJ_EXT) \
    $(DPO)usr_rates.$(OBJ_EXT) \
    $(DPO)usr_write_out0.$(OBJ_EXT) \
    $(DPO)usr_write_out1.$(OBJ_EXT) \
    $(DPO)utilities.$(OBJ_EXT) \
    $(DPO)vavg_u_g.$(OBJ_EXT) \
    $(DPO)vavg_u_s.$(OBJ_EXT) \
    $(DPO)vavg_v_g.$(OBJ_EXT) \
    $(DPO)vavg_v_s.$(OBJ_EXT) \
    $(DPO)vavg_w_g.$(OBJ_EXT) \
    $(DPO)vavg_w_s.$(OBJ_EXT) \
    $(DPO)vf_gs_x.$(OBJ_EXT) \
    $(DPO)vf_gs_y.$(OBJ_EXT) \
    $(DPO)vf_gs_z.$(OBJ_EXT) \
    $(DPO)vtc_scalar.$(OBJ_EXT) \
    $(DPO)write_ab_m.$(OBJ_EXT) \
    $(DPO)write_ab_m_var.$(OBJ_EXT) \
    $(DPO)write_error.$(OBJ_EXT) \
    $(DPO)write_header.$(OBJ_EXT) \
    $(DPO)write_out0.$(OBJ_EXT) \
    $(DPO)write_out1.$(OBJ_EXT) \
    $(DPO)write_out3.$(OBJ_EXT) \
    $(DPO)write_res0.$(OBJ_EXT) \
    $(DPO)write_res1.$(OBJ_EXT) \
    $(DPO)write_spx0.$(OBJ_EXT) \
    $(DPO)write_spx1.$(OBJ_EXT) \
    $(DPO)write_table.$(OBJ_EXT) \
    $(DPO)write_usr0.$(OBJ_EXT) \
    $(DPO)write_usr1.$(OBJ_EXT) \
    $(DPO)xerbla.$(OBJ_EXT) \
    $(DPO)zero_array.$(OBJ_EXT) \
    $(DPO)zero_norm_vel.$(OBJ_EXT) \
    $(DPO)allocate_cut_cell_arrays.$(OBJ_EXT) \
    $(DPO)allocate_dummy_cut_cell_arrays.$(OBJ_EXT) \
    $(DPO)calc_vort_out.$(OBJ_EXT) \
    $(DPO)cartesian_grid_init_namelist.$(OBJ_EXT) \
    $(DPO)CG_set_bc0.$(OBJ_EXT) \
    $(DPO)CG_set_outflow.$(OBJ_EXT) \
    $(DPO)CG_source_u_g.$(OBJ_EXT) \
    $(DPO)CG_source_u_s.$(OBJ_EXT) \
    $(DPO)CG_source_v_g.$(OBJ_EXT) \
    $(DPO)CG_source_v_s.$(OBJ_EXT) \
    $(DPO)CG_source_w_g.$(OBJ_EXT) \
    $(DPO)CG_source_w_s.$(OBJ_EXT) \
    $(DPO)check_data_cartesian.$(OBJ_EXT) \
    $(DPO)cut_cell_preprocessing.$(OBJ_EXT) \
    $(DPO)deallocate_cut_cell_arrays.$(OBJ_EXT) \
    $(DPO)define_quadrics.$(OBJ_EXT) \
    $(DPO)dmp_cartesian.$(OBJ_EXT) \
    $(DPO)eval_usr_fct.$(OBJ_EXT) \
    $(DPO)get_alpha.$(OBJ_EXT) \
    $(DPO)get_connectivity.$(OBJ_EXT) \
    $(DPO)get_cut_cell_flags.$(OBJ_EXT) \
    $(DPO)get_cut_cell_volume_area.$(OBJ_EXT) \
    $(DPO)get_delh.$(OBJ_EXT) \
    $(DPO)get_master.$(OBJ_EXT) \
    $(DPO)get_poly_data.$(OBJ_EXT) \
    $(DPO)get_stl_data.$(OBJ_EXT) \
    $(DPO)set_Odxyz.$(OBJ_EXT) \
    $(DPO)update_dashboard.$(OBJ_EXT) \
    $(DPO)vtk_out.$(OBJ_EXT) \
    $(DPO)write_progress_bar.$(OBJ_EXT) \
    $(DPO)check_axis.$(OBJ_EXT) \
    $(DPO)check_bc_geometry.$(OBJ_EXT) \
    $(DPO)check_bc_inflow.$(OBJ_EXT) \
    $(DPO)check_bc_outflow.$(OBJ_EXT) \
    $(DPO)check_bc_walls.$(OBJ_EXT) \
    $(DPO)check_boundary_conditions.$(OBJ_EXT) \
    $(DPO)check_dmp_prereqs.$(OBJ_EXT) \
    $(DPO)check_gas_phase.$(OBJ_EXT) \
    $(DPO)check_geometry.$(OBJ_EXT) \
    $(DPO)check_geometry_prereqs.$(OBJ_EXT) \
    $(DPO)check_ic_common_discrete.$(OBJ_EXT) \
    $(DPO)check_initial_conditions_dem.$(OBJ_EXT) \
    $(DPO)check_initial_conditions.$(OBJ_EXT) \
    $(DPO)check_initial_conditions_mppic.$(OBJ_EXT) \
    $(DPO)check_internal_surfaces.$(OBJ_EXT) \
    $(DPO)check_numerics.$(OBJ_EXT) \
    $(DPO)check_output_control.$(OBJ_EXT) \
    $(DPO)check_point_sources.$(OBJ_EXT) \
    $(DPO)check_run_control.$(OBJ_EXT) \
    $(DPO)check_solids_common_all.$(OBJ_EXT) \
    $(DPO)check_solids_common_discrete.$(OBJ_EXT) \
    $(DPO)check_solids_continuum.$(OBJ_EXT) \
    $(DPO)check_solids_des.$(OBJ_EXT) \
    $(DPO)check_solids_model_prereqs.$(OBJ_EXT) \
    $(DPO)check_solids_mppic.$(OBJ_EXT) \
    $(DPO)check_solids_phases.$(OBJ_EXT) \
    $(DPO)check_data_odepack.$(OBJ_EXT) \
    $(DPO)stiff_chem_rrates.$(OBJ_EXT) \
    $(DPO)calc_force_des_cutcell.$(OBJ_EXT) \
    $(DPO)calc_force_des.$(OBJ_EXT) \
    $(DPO)calc_rrate_des.$(OBJ_EXT) \
    $(DPO)calc_thermo_des.$(OBJ_EXT) \
    $(DPO)cfassign.$(OBJ_EXT) \
    $(DPO)cffctowall.$(OBJ_EXT) \
    $(DPO)cffctow.$(OBJ_EXT) \
    $(DPO)cfnewvalues.$(OBJ_EXT) \
    $(DPO)cfrelvel.$(OBJ_EXT) \
    $(DPO)cfslide.$(OBJ_EXT) \
    $(DPO)cfslidewall.$(OBJ_EXT) \
    $(DPO)cfupdateold.$(OBJ_EXT) \
    $(DPO)cfwallcontact.$(OBJ_EXT) \
    $(DPO)cfwallposvel.$(OBJ_EXT) \
    $(DPO)check_des_bc.$(OBJ_EXT) \
    $(DPO)check_des_cohesion.$(OBJ_EXT) \
    $(DPO)check_des_collision.$(OBJ_EXT) \
    $(DPO)check_des_coupling.$(OBJ_EXT) \
    $(DPO)check_des_data.$(OBJ_EXT) \
    $(DPO)check_des_energy.$(OBJ_EXT) \
    $(DPO)check_des_geometry.$(OBJ_EXT) \
    $(DPO)check_des_hybrid.$(OBJ_EXT) \
    $(DPO)check_des_mppic.$(OBJ_EXT) \
    $(DPO)check_des_rxns.$(OBJ_EXT) \
    $(DPO)check_des_thermo.$(OBJ_EXT) \
    $(DPO)des_allocate_arrays.$(OBJ_EXT) \
    $(DPO)des_check_particle.$(OBJ_EXT) \
    $(DPO)des_cluster_identification.$(OBJ_EXT) \
    $(DPO)des_functions.$(OBJ_EXT) \
    $(DPO)des_granular_temperature.$(OBJ_EXT) \
    $(DPO)des_init_arrays.$(OBJ_EXT) \
    $(DPO)des_init_bc.$(OBJ_EXT) \
    $(DPO)des_init_namelist.$(OBJ_EXT) \
    $(DPO)des_mass_inlet.$(OBJ_EXT) \
    $(DPO)des_physical_prop.$(OBJ_EXT) \
    $(DPO)des_reaction_model.$(OBJ_EXT) \
    $(DPO)des_rrates0.$(OBJ_EXT) \
    $(DPO)des_set_ic.$(OBJ_EXT) \
    $(DPO)des_thermo_cond.$(OBJ_EXT) \
    $(DPO)des_thermo_conv.$(OBJ_EXT) \
    $(DPO)des_thermo_newvalues.$(OBJ_EXT) \
    $(DPO)des_thermo_rad.$(OBJ_EXT) \
    $(DPO)des_time_march.$(OBJ_EXT) \
    $(DPO)des_wallbc_preprocessing.$(OBJ_EXT) \
    $(DPO)drag_fgs.$(OBJ_EXT) \
    $(DPO)gas_drag.$(OBJ_EXT) \
    $(DPO)generate_particle_config.$(OBJ_EXT) \
    $(DPO)make_arrays_des.$(OBJ_EXT) \
    $(DPO)mppic_routines.$(OBJ_EXT) \
    $(DPO)neighbour.$(OBJ_EXT) \
    $(DPO)nsquare.$(OBJ_EXT) \
    $(DPO)particles_in_cell.$(OBJ_EXT) \
    $(DPO)read_des_restart.$(OBJ_EXT) \
    $(DPO)solid_drag.$(OBJ_EXT) \
    $(DPO)usr0_des.$(OBJ_EXT) \
    $(DPO)usr1_des.$(OBJ_EXT) \
    $(DPO)usr2_des.$(OBJ_EXT) \
    $(DPO)usr3_des.$(OBJ_EXT) \
    $(DPO)usr4_des.$(OBJ_EXT) \
    $(DPO)usr_rates_des.$(OBJ_EXT) \
    $(DPO)walledgecontact.$(OBJ_EXT) \
    $(DPO)wallfacecontact.$(OBJ_EXT) \
    $(DPO)wallnodecontact.$(OBJ_EXT) \
    $(DPO)write_des_data.$(OBJ_EXT) \
    $(DPO)write_des_restart.$(OBJ_EXT) \
    $(DPO)gaussj.$(OBJ_EXT) \
    $(DPO)odeint.$(OBJ_EXT) \
    $(DPO)rkck.$(OBJ_EXT) \
    $(DPO)rkqs.$(OBJ_EXT) \
    $(DPO)source_population_eq.$(OBJ_EXT) \
    $(DPO)usr_dqmom.$(OBJ_EXT) \
    $(DPO)adjust_eps_ghd.$(OBJ_EXT) \
    $(DPO)bulk_viscosity.$(OBJ_EXT) \
    $(DPO)calc_d_ghd.$(OBJ_EXT) \
    $(DPO)calc_external_forces.$(OBJ_EXT) \
    $(DPO)calc_nflux.$(OBJ_EXT) \
    $(DPO)chi_ij_GHD.$(OBJ_EXT) \
    $(DPO)cooling_rate.$(OBJ_EXT) \
    $(DPO)cooling_rate_tc.$(OBJ_EXT) \
    $(DPO)dufour_coeff.$(OBJ_EXT) \
    $(DPO)ghd.$(OBJ_EXT) \
    $(DPO)ghdmassflux.$(OBJ_EXT) \
    $(DPO)mass_mobility.$(OBJ_EXT) \
    $(DPO)ordinary_diff.$(OBJ_EXT) \
    $(DPO)partial_elim_ghd.$(OBJ_EXT) \
    $(DPO)pressure.$(OBJ_EXT) \
    $(DPO)shear_viscosity.$(OBJ_EXT) \
    $(DPO)source_ghd_granular_energy.$(OBJ_EXT) \
    $(DPO)thermal_conductivity.$(OBJ_EXT) \
    $(DPO)thermal_diffusivity.$(OBJ_EXT) \
    $(DPO)thermal_mobility.$(OBJ_EXT) \
    $(DPO)transport_coeff_ghd.$(OBJ_EXT) \
    $(DPO)qmomk_allocate_arrays.$(OBJ_EXT) \
    $(DPO)qmomk_gas_drag.$(OBJ_EXT) \
    $(DPO)qmomk_init_bc.$(OBJ_EXT) \
    $(DPO)qmomk_initial_conditions.$(OBJ_EXT) \
    $(DPO)qmomk_init_namelist.$(OBJ_EXT) \
    $(DPO)qmomk_make_arrays.$(OBJ_EXT) \
    $(DPO)qmomk_read_restart.$(OBJ_EXT) \
    $(DPO)qmomk_set_bc.$(OBJ_EXT) \
    $(DPO)qmomk_time_march.$(OBJ_EXT) \
    $(DPO)qmomk_write_restart.$(OBJ_EXT) \
    $(DPO)get_values.$(OBJ_EXT) \
    $(DPO)readTherm.$(OBJ_EXT) \
    $(DPO)blas90.a $(DPO)odepack.a $(DPO)dgtsv90.a
	$(LINK_CMD) $(LINK_FLAGS) \
    $(DPO)accum_resid.$(OBJ_EXT) \
    $(DPO)adjust_a_u_g.$(OBJ_EXT) \
    $(DPO)adjust_a_u_s.$(OBJ_EXT) \
    $(DPO)adjust_a_v_g.$(OBJ_EXT) \
    $(DPO)adjust_a_v_s.$(OBJ_EXT) \
    $(DPO)adjust_a_w_g.$(OBJ_EXT) \
    $(DPO)adjust_a_w_s.$(OBJ_EXT) \
    $(DPO)adjust_dt.$(OBJ_EXT) \
    $(DPO)adjust_eps.$(OBJ_EXT) \
    $(DPO)adjust_leq.$(OBJ_EXT) \
    $(DPO)adjust_rop.$(OBJ_EXT) \
    $(DPO)adjust_theta.$(OBJ_EXT) \
    $(DPO)allocate_arrays.$(OBJ_EXT) \
    $(DPO)ambm_mod.$(OBJ_EXT) \
    $(DPO)bc_mod.$(OBJ_EXT) \
    $(DPO)bc_phi.$(OBJ_EXT) \
    $(DPO)bc_theta.$(OBJ_EXT) \
    $(DPO)b_m_p_star.$(OBJ_EXT) \
    $(DPO)boundfunijk3_mod.$(OBJ_EXT) \
    $(DPO)boundfunijk_mod.$(OBJ_EXT) \
    $(DPO)bound_x.$(OBJ_EXT) \
    $(DPO)calc_cell.$(OBJ_EXT) \
    $(DPO)calc_coeff.$(OBJ_EXT) \
    $(DPO)calc_d.$(OBJ_EXT) \
    $(DPO)calc_dif_g.$(OBJ_EXT) \
    $(DPO)calc_dif_s.$(OBJ_EXT) \
    $(DPO)calc_drag.$(OBJ_EXT) \
    $(DPO)calc_e.$(OBJ_EXT) \
    $(DPO)calc_gama.$(OBJ_EXT) \
    $(DPO)calc_grbdry.$(OBJ_EXT) \
    $(DPO)calc_h.$(OBJ_EXT) \
    $(DPO)calc_k_cp.$(OBJ_EXT) \
    $(DPO)calc_k_g.$(OBJ_EXT) \
    $(DPO)calc_k_s.$(OBJ_EXT) \
    $(DPO)calc_mflux.$(OBJ_EXT) \
    $(DPO)calc_mu_g.$(OBJ_EXT) \
    $(DPO)calc_mu_s.$(OBJ_EXT) \
    $(DPO)calc_mw.$(OBJ_EXT) \
    $(DPO)calc_outflow.$(OBJ_EXT) \
    $(DPO)calc_p_star.$(OBJ_EXT) \
    $(DPO)calc_resid.$(OBJ_EXT) \
    $(DPO)calc_s_ddot_s.$(OBJ_EXT) \
    $(DPO)calc_trd_g.$(OBJ_EXT) \
    $(DPO)calc_trd_s.$(OBJ_EXT) \
    $(DPO)calc_u_friction.$(OBJ_EXT) \
    $(DPO)calc_vol_fr.$(OBJ_EXT) \
    $(DPO)calc_xsi.$(OBJ_EXT) \
    $(DPO)cal_d.$(OBJ_EXT) \
    $(DPO)cdist_mod.$(OBJ_EXT) \
    $(DPO)check_ab_m.$(OBJ_EXT) \
    $(DPO)check_convergence.$(OBJ_EXT) \
    $(DPO)check_data_03.$(OBJ_EXT) \
    $(DPO)check_data_06.$(OBJ_EXT) \
    $(DPO)check_data_07.$(OBJ_EXT) \
    $(DPO)check_data_08.$(OBJ_EXT) \
    $(DPO)check_data_09.$(OBJ_EXT) \
    $(DPO)check_data_10.$(OBJ_EXT) \
    $(DPO)check_data_20.$(OBJ_EXT) \
    $(DPO)check_data_30.$(OBJ_EXT) \
    $(DPO)check_mass_balance.$(OBJ_EXT) \
    $(DPO)check_mod.$(OBJ_EXT) \
    $(DPO)check_one_axis.$(OBJ_EXT) \
    $(DPO)check_plane.$(OBJ_EXT) \
    $(DPO)chischeme_mod.$(OBJ_EXT) \
    $(DPO)cn_extrapol.$(OBJ_EXT) \
    $(DPO)coeff_mod.$(OBJ_EXT) \
    $(DPO)compare.$(OBJ_EXT) \
    $(DPO)constant_mod.$(OBJ_EXT) \
    $(DPO)cont_mod.$(OBJ_EXT) \
    $(DPO)conv_dif_phi.$(OBJ_EXT) \
    $(DPO)conv_dif_u_g.$(OBJ_EXT) \
    $(DPO)conv_dif_u_s.$(OBJ_EXT) \
    $(DPO)conv_dif_v_g.$(OBJ_EXT) \
    $(DPO)conv_dif_v_s.$(OBJ_EXT) \
    $(DPO)conv_dif_w_g.$(OBJ_EXT) \
    $(DPO)conv_dif_w_s.$(OBJ_EXT) \
    $(DPO)conv_pp_g.$(OBJ_EXT) \
    $(DPO)conv_rop.$(OBJ_EXT) \
    $(DPO)conv_rop_g.$(OBJ_EXT) \
    $(DPO)conv_rop_s.$(OBJ_EXT) \
    $(DPO)conv_source_epp.$(OBJ_EXT) \
    $(DPO)copy_a.$(OBJ_EXT) \
    $(DPO)corner.$(OBJ_EXT) \
    $(DPO)corner_mod.$(OBJ_EXT) \
    $(DPO)correct_0.$(OBJ_EXT) \
    $(DPO)correct_1.$(OBJ_EXT) \
    $(DPO)dbg_mod.$(OBJ_EXT) \
    $(DPO)dgtsl.$(OBJ_EXT) \
    $(DPO)dif_u_is.$(OBJ_EXT) \
    $(DPO)dif_v_is.$(OBJ_EXT) \
    $(DPO)dif_w_is.$(OBJ_EXT) \
    $(DPO)discretize.$(OBJ_EXT) \
    $(DPO)display_resid.$(OBJ_EXT) \
    $(DPO)drag_gs.$(OBJ_EXT) \
    $(DPO)drag_mod.$(OBJ_EXT) \
    $(DPO)drag_ss.$(OBJ_EXT) \
    $(DPO)energy_mod.$(OBJ_EXT) \
    $(DPO)eosg.$(OBJ_EXT) \
    $(DPO)eoss.$(OBJ_EXT) \
    $(DPO)equal.$(OBJ_EXT) \
    $(DPO)error_manager_mod.$(OBJ_EXT) \
    $(DPO)error_routine.$(OBJ_EXT) \
    $(DPO)exchange.$(OBJ_EXT) \
    $(DPO)exit.$(OBJ_EXT) \
    $(DPO)fldvar_mod.$(OBJ_EXT) \
    $(DPO)flow_to_vel.$(OBJ_EXT) \
    $(DPO)function_mod.$(OBJ_EXT) \
    $(DPO)funits_mod.$(OBJ_EXT) \
    $(DPO)g_0.$(OBJ_EXT) \
    $(DPO)geometry_mod.$(OBJ_EXT) \
    $(DPO)get_bc_area.$(OBJ_EXT) \
    $(DPO)get_data.$(OBJ_EXT) \
    $(DPO)get_eq.$(OBJ_EXT) \
    $(DPO)get_flow_bc.$(OBJ_EXT) \
    $(DPO)get_hloss.$(OBJ_EXT) \
    $(DPO)get_is.$(OBJ_EXT) \
    $(DPO)get_philoss.$(OBJ_EXT) \
    $(DPO)get_ps.$(OBJ_EXT) \
    $(DPO)get_smass.$(OBJ_EXT) \
    $(DPO)get_stats.$(OBJ_EXT) \
    $(DPO)get_walls_bc.$(OBJ_EXT) \
    $(DPO)ic_mod.$(OBJ_EXT) \
    $(DPO)in_bin_512.$(OBJ_EXT) \
    $(DPO)in_bin_512i.$(OBJ_EXT) \
    $(DPO)indices_mod.$(OBJ_EXT) \
    $(DPO)init_ab_m.$(OBJ_EXT) \
    $(DPO)init_fvars.$(OBJ_EXT) \
    $(DPO)init_namelist.$(OBJ_EXT) \
    $(DPO)init_resid.$(OBJ_EXT) \
    $(DPO)is_mod.$(OBJ_EXT) \
    $(DPO)iterate.$(OBJ_EXT) \
    $(DPO)k_epsilon_prop.$(OBJ_EXT) \
    $(DPO)kintheory_drag_ss.$(OBJ_EXT) \
    $(DPO)kintheory_energy_dissipation_ss.$(OBJ_EXT) \
    $(DPO)kintheory_mod.$(OBJ_EXT) \
    $(DPO)kintheory_u_s.$(OBJ_EXT) \
    $(DPO)kintheory_v_s.$(OBJ_EXT) \
    $(DPO)kintheory_w_s.$(OBJ_EXT) \
    $(DPO)leq_bicgs.$(OBJ_EXT) \
    $(DPO)leq_bicgst.$(OBJ_EXT) \
    $(DPO)leq_cg.$(OBJ_EXT) \
    $(DPO)leq_gmres.$(OBJ_EXT) \
    $(DPO)leqsol_mod.$(OBJ_EXT) \
    $(DPO)leq_sor.$(OBJ_EXT) \
    $(DPO)line_too_big.$(OBJ_EXT) \
    $(DPO)location_check.$(OBJ_EXT) \
    $(DPO)location.$(OBJ_EXT) \
    $(DPO)machine.$(OBJ_EXT) \
    $(DPO)machine_mod.$(OBJ_EXT) \
    $(DPO)make_upper_case.$(OBJ_EXT) \
    $(DPO)mark_phase_4_cor.$(OBJ_EXT) \
    $(DPO)matrix_mod.$(OBJ_EXT) \
    $(DPO)mfix.$(OBJ_EXT) \
    $(DPO)mfix_netcdf_mod.$(OBJ_EXT) \
    $(DPO)mflux_mod.$(OBJ_EXT) \
    $(DPO)mms_mod.$(OBJ_EXT) \
    $(DPO)mod_bc_i.$(OBJ_EXT) \
    $(DPO)mod_bc_j.$(OBJ_EXT) \
    $(DPO)mod_bc_k.$(OBJ_EXT) \
    $(DPO)open_file.$(OBJ_EXT) \
    $(DPO)open_files.$(OBJ_EXT) \
    $(DPO)out_array_c.$(OBJ_EXT) \
    $(DPO)out_array.$(OBJ_EXT) \
    $(DPO)out_array_kc.$(OBJ_EXT) \
    $(DPO)out_array_k.$(OBJ_EXT) \
    $(DPO)out_bin_512.$(OBJ_EXT) \
    $(DPO)out_bin_512i.$(OBJ_EXT) \
    $(DPO)out_bin_512r.$(OBJ_EXT) \
    $(DPO)out_bin_r.$(OBJ_EXT) \
    $(DPO)output_mod.$(OBJ_EXT) \
    $(DPO)parallel_mod.$(OBJ_EXT) \
    $(DPO)param1_mod.$(OBJ_EXT) \
    $(DPO)param_mod.$(OBJ_EXT) \
    $(DPO)parse_line.$(OBJ_EXT) \
    $(DPO)parse_mod.$(OBJ_EXT) \
    $(DPO)parse_resid_string.$(OBJ_EXT) \
    $(DPO)parse_rxn.$(OBJ_EXT) \
    $(DPO)partial_elim.$(OBJ_EXT) \
    $(DPO)pgcor_mod.$(OBJ_EXT) \
    $(DPO)physical_prop.$(OBJ_EXT) \
    $(DPO)physprop_mod.$(OBJ_EXT) \
    $(DPO)pscor_mod.$(OBJ_EXT) \
    $(DPO)ps_mod.$(OBJ_EXT) \
    $(DPO)read_database.$(OBJ_EXT) \
    $(DPO)read_namelist.$(OBJ_EXT) \
    $(DPO)read_res0.$(OBJ_EXT) \
    $(DPO)read_res1.$(OBJ_EXT) \
    $(DPO)remove_comment.$(OBJ_EXT) \
    $(DPO)reset_new.$(OBJ_EXT) \
    $(DPO)residual_mod.$(OBJ_EXT) \
    $(DPO)rrates0.$(OBJ_EXT) \
    $(DPO)rrates.$(OBJ_EXT) \
    $(DPO)rrates_init.$(OBJ_EXT) \
    $(DPO)run_mod.$(OBJ_EXT) \
    $(DPO)rxn_com_mod.$(OBJ_EXT) \
    $(DPO)rxns_mod.$(OBJ_EXT) \
    $(DPO)scalar_prop.$(OBJ_EXT) \
    $(DPO)scalars_mod.$(OBJ_EXT) \
    $(DPO)scales_mod.$(OBJ_EXT) \
    $(DPO)seek_comment.$(OBJ_EXT) \
    $(DPO)seek_end.$(OBJ_EXT) \
    $(DPO)set_bc0.$(OBJ_EXT) \
    $(DPO)set_bc1.$(OBJ_EXT) \
    $(DPO)set_bc_flow.$(OBJ_EXT) \
    $(DPO)set_constants.$(OBJ_EXT) \
    $(DPO)set_constprop.$(OBJ_EXT) \
    $(DPO)set_flags.$(OBJ_EXT) \
    $(DPO)set_fluidbed_p.$(OBJ_EXT) \
    $(DPO)set_geometry1.$(OBJ_EXT) \
    $(DPO)set_geometry.$(OBJ_EXT) \
    $(DPO)set_icbc_flags.$(OBJ_EXT) \
    $(DPO)set_ic.$(OBJ_EXT) \
    $(DPO)set_increments3.$(OBJ_EXT) \
    $(DPO)set_increments.$(OBJ_EXT) \
    $(DPO)set_index1a3.$(OBJ_EXT) \
    $(DPO)set_index1a.$(OBJ_EXT) \
    $(DPO)set_index1.$(OBJ_EXT) \
    $(DPO)set_l_scale.$(OBJ_EXT) \
    $(DPO)set_max2.$(OBJ_EXT) \
    $(DPO)set_mw_mix_g.$(OBJ_EXT) \
    $(DPO)set_outflow.$(OBJ_EXT) \
    $(DPO)set_ps.$(OBJ_EXT) \
    $(DPO)set_ro_g.$(OBJ_EXT) \
    $(DPO)set_ro_s.$(OBJ_EXT) \
    $(DPO)set_wall_bc.$(OBJ_EXT) \
    $(DPO)shift_dxyz.$(OBJ_EXT) \
    $(DPO)solve_continuity.$(OBJ_EXT) \
    $(DPO)solve_energy_eq.$(OBJ_EXT) \
    $(DPO)solve_epp.$(OBJ_EXT) \
    $(DPO)solve_granular_energy.$(OBJ_EXT) \
    $(DPO)solve_k_epsilon_eq.$(OBJ_EXT) \
    $(DPO)solve_lin_eq.$(OBJ_EXT) \
    $(DPO)solve_pp_g.$(OBJ_EXT) \
    $(DPO)solve_scalar_eq.$(OBJ_EXT) \
    $(DPO)solve_species_eq.$(OBJ_EXT) \
    $(DPO)solve_vel_star.$(OBJ_EXT) \
    $(DPO)source_granular_energy.$(OBJ_EXT) \
    $(DPO)source_phi.$(OBJ_EXT) \
    $(DPO)source_pp_g.$(OBJ_EXT) \
    $(DPO)source_rop_g.$(OBJ_EXT) \
    $(DPO)source_rop_s.$(OBJ_EXT) \
    $(DPO)source_u_g.$(OBJ_EXT) \
    $(DPO)source_u_s.$(OBJ_EXT) \
    $(DPO)source_v_g.$(OBJ_EXT) \
    $(DPO)source_v_s.$(OBJ_EXT) \
    $(DPO)source_w_g.$(OBJ_EXT) \
    $(DPO)source_w_s.$(OBJ_EXT) \
    $(DPO)tau_g_mod.$(OBJ_EXT) \
    $(DPO)tau_s_mod.$(OBJ_EXT) \
    $(DPO)tau_u_g.$(OBJ_EXT) \
    $(DPO)tau_u_s.$(OBJ_EXT) \
    $(DPO)tau_v_g.$(OBJ_EXT) \
    $(DPO)tau_v_s.$(OBJ_EXT) \
    $(DPO)tau_w_g.$(OBJ_EXT) \
    $(DPO)tau_w_s.$(OBJ_EXT) \
    $(DPO)test_lin_eq.$(OBJ_EXT) \
    $(DPO)time_cpu_mod.$(OBJ_EXT) \
    $(DPO)time_march.$(OBJ_EXT) \
    $(DPO)tmp_array1_mod.$(OBJ_EXT) \
    $(DPO)tmp_array_mod.$(OBJ_EXT) \
    $(DPO)toleranc_mod.$(OBJ_EXT) \
    $(DPO)trace_mod.$(OBJ_EXT) \
    $(DPO)transfer.$(OBJ_EXT) \
    $(DPO)transport_prop.$(OBJ_EXT) \
    $(DPO)turb_mod.$(OBJ_EXT) \
    $(DPO)undef_2_0.$(OBJ_EXT) \
    $(DPO)under_relax.$(OBJ_EXT) \
    $(DPO)update_old.$(OBJ_EXT) \
    $(DPO)ur_facs_mod.$(OBJ_EXT) \
    $(DPO)usr0.$(OBJ_EXT) \
    $(DPO)usr1.$(OBJ_EXT) \
    $(DPO)usr2.$(OBJ_EXT) \
    $(DPO)usr3.$(OBJ_EXT) \
    $(DPO)usr_init_namelist.$(OBJ_EXT) \
    $(DPO)usr_mod.$(OBJ_EXT) \
    $(DPO)usr_rates.$(OBJ_EXT) \
    $(DPO)usr_write_out0.$(OBJ_EXT) \
    $(DPO)usr_write_out1.$(OBJ_EXT) \
    $(DPO)utilities.$(OBJ_EXT) \
    $(DPO)vavg_u_g.$(OBJ_EXT) \
    $(DPO)vavg_u_s.$(OBJ_EXT) \
    $(DPO)vavg_v_g.$(OBJ_EXT) \
    $(DPO)vavg_v_s.$(OBJ_EXT) \
    $(DPO)vavg_w_g.$(OBJ_EXT) \
    $(DPO)vavg_w_s.$(OBJ_EXT) \
    $(DPO)vf_gs_x.$(OBJ_EXT) \
    $(DPO)vf_gs_y.$(OBJ_EXT) \
    $(DPO)vf_gs_z.$(OBJ_EXT) \
    $(DPO)visc_g_mod.$(OBJ_EXT) \
    $(DPO)visc_s_mod.$(OBJ_EXT) \
    $(DPO)vshear_mod.$(OBJ_EXT) \
    $(DPO)vtc_scalar.$(OBJ_EXT) \
    $(DPO)write_ab_m.$(OBJ_EXT) \
    $(DPO)write_ab_m_var.$(OBJ_EXT) \
    $(DPO)write_error.$(OBJ_EXT) \
    $(DPO)write_header.$(OBJ_EXT) \
    $(DPO)write_out0.$(OBJ_EXT) \
    $(DPO)write_out1.$(OBJ_EXT) \
    $(DPO)write_out3.$(OBJ_EXT) \
    $(DPO)write_res0.$(OBJ_EXT) \
    $(DPO)write_res1.$(OBJ_EXT) \
    $(DPO)write_spx0.$(OBJ_EXT) \
    $(DPO)write_spx1.$(OBJ_EXT) \
    $(DPO)write_table.$(OBJ_EXT) \
    $(DPO)write_usr0.$(OBJ_EXT) \
    $(DPO)write_usr1.$(OBJ_EXT) \
    $(DPO)xerbla.$(OBJ_EXT) \
    $(DPO)xsi_array_mod.$(OBJ_EXT) \
    $(DPO)zero_array.$(OBJ_EXT) \
    $(DPO)zero_norm_vel.$(OBJ_EXT) \
    $(DPO)allocate_cut_cell_arrays.$(OBJ_EXT) \
    $(DPO)allocate_dummy_cut_cell_arrays.$(OBJ_EXT) \
    $(DPO)calc_vort_out.$(OBJ_EXT) \
    $(DPO)cartesian_grid_init_namelist.$(OBJ_EXT) \
    $(DPO)CG_set_bc0.$(OBJ_EXT) \
    $(DPO)CG_set_outflow.$(OBJ_EXT) \
    $(DPO)CG_source_u_g.$(OBJ_EXT) \
    $(DPO)CG_source_u_s.$(OBJ_EXT) \
    $(DPO)CG_source_v_g.$(OBJ_EXT) \
    $(DPO)CG_source_v_s.$(OBJ_EXT) \
    $(DPO)CG_source_w_g.$(OBJ_EXT) \
    $(DPO)CG_source_w_s.$(OBJ_EXT) \
    $(DPO)check_data_cartesian.$(OBJ_EXT) \
    $(DPO)cutcell_mod.$(OBJ_EXT) \
    $(DPO)cut_cell_preprocessing.$(OBJ_EXT) \
    $(DPO)dashboard_mod.$(OBJ_EXT) \
    $(DPO)deallocate_cut_cell_arrays.$(OBJ_EXT) \
    $(DPO)define_quadrics.$(OBJ_EXT) \
    $(DPO)dmp_cartesian.$(OBJ_EXT) \
    $(DPO)eval_usr_fct.$(OBJ_EXT) \
    $(DPO)get_alpha.$(OBJ_EXT) \
    $(DPO)get_connectivity.$(OBJ_EXT) \
    $(DPO)get_cut_cell_flags.$(OBJ_EXT) \
    $(DPO)get_cut_cell_volume_area.$(OBJ_EXT) \
    $(DPO)get_delh.$(OBJ_EXT) \
    $(DPO)get_master.$(OBJ_EXT) \
    $(DPO)get_poly_data.$(OBJ_EXT) \
    $(DPO)get_stl_data.$(OBJ_EXT) \
    $(DPO)polygon_mod.$(OBJ_EXT) \
    $(DPO)progress_bar_mod.$(OBJ_EXT) \
    $(DPO)quadric_mod.$(OBJ_EXT) \
    $(DPO)set_Odxyz.$(OBJ_EXT) \
    $(DPO)stl_mod.$(OBJ_EXT) \
    $(DPO)update_dashboard.$(OBJ_EXT) \
    $(DPO)vtk_mod.$(OBJ_EXT) \
    $(DPO)vtk_out.$(OBJ_EXT) \
    $(DPO)write_progress_bar.$(OBJ_EXT) \
    $(DPO)check_axis.$(OBJ_EXT) \
    $(DPO)check_bc_geometry.$(OBJ_EXT) \
    $(DPO)check_bc_inflow.$(OBJ_EXT) \
    $(DPO)check_bc_outflow.$(OBJ_EXT) \
    $(DPO)check_bc_walls.$(OBJ_EXT) \
    $(DPO)check_boundary_conditions.$(OBJ_EXT) \
    $(DPO)check_dmp_prereqs.$(OBJ_EXT) \
    $(DPO)check_gas_phase.$(OBJ_EXT) \
    $(DPO)check_geometry.$(OBJ_EXT) \
    $(DPO)check_geometry_prereqs.$(OBJ_EXT) \
    $(DPO)check_ic_common_discrete.$(OBJ_EXT) \
    $(DPO)check_initial_conditions_dem.$(OBJ_EXT) \
    $(DPO)check_initial_conditions.$(OBJ_EXT) \
    $(DPO)check_initial_conditions_mppic.$(OBJ_EXT) \
    $(DPO)check_internal_surfaces.$(OBJ_EXT) \
    $(DPO)check_numerics.$(OBJ_EXT) \
    $(DPO)check_output_control.$(OBJ_EXT) \
    $(DPO)check_point_sources.$(OBJ_EXT) \
    $(DPO)check_run_control.$(OBJ_EXT) \
    $(DPO)check_solids_common_all.$(OBJ_EXT) \
    $(DPO)check_solids_common_discrete.$(OBJ_EXT) \
    $(DPO)check_solids_continuum.$(OBJ_EXT) \
    $(DPO)check_solids_des.$(OBJ_EXT) \
    $(DPO)check_solids_model_prereqs.$(OBJ_EXT) \
    $(DPO)check_solids_mppic.$(OBJ_EXT) \
    $(DPO)check_solids_phases.$(OBJ_EXT) \
    $(DPO)check_data_odepack.$(OBJ_EXT) \
    $(DPO)stiff_chem_dbg_mod.$(OBJ_EXT) \
    $(DPO)stiff_chem_maps_mod.$(OBJ_EXT) \
    $(DPO)stiff_chem_mod.$(OBJ_EXT) \
    $(DPO)stiff_chem_rrates.$(OBJ_EXT) \
    $(DPO)stiff_chem_stats_mod.$(OBJ_EXT) \
    $(DPO)calc_force_des_cutcell.$(OBJ_EXT) \
    $(DPO)calc_force_des.$(OBJ_EXT) \
    $(DPO)calc_rrate_des.$(OBJ_EXT) \
    $(DPO)calc_thermo_des.$(OBJ_EXT) \
    $(DPO)cfassign.$(OBJ_EXT) \
    $(DPO)cffctowall.$(OBJ_EXT) \
    $(DPO)cffctow.$(OBJ_EXT) \
    $(DPO)cfnewvalues.$(OBJ_EXT) \
    $(DPO)cfrelvel.$(OBJ_EXT) \
    $(DPO)cfslide.$(OBJ_EXT) \
    $(DPO)cfslidewall.$(OBJ_EXT) \
    $(DPO)cfupdateold.$(OBJ_EXT) \
    $(DPO)cfwallcontact.$(OBJ_EXT) \
    $(DPO)cfwallposvel.$(OBJ_EXT) \
    $(DPO)check_des_bc.$(OBJ_EXT) \
    $(DPO)check_des_cohesion.$(OBJ_EXT) \
    $(DPO)check_des_collision.$(OBJ_EXT) \
    $(DPO)check_des_coupling.$(OBJ_EXT) \
    $(DPO)check_des_data.$(OBJ_EXT) \
    $(DPO)check_des_energy.$(OBJ_EXT) \
    $(DPO)check_des_geometry.$(OBJ_EXT) \
    $(DPO)check_des_hybrid.$(OBJ_EXT) \
    $(DPO)check_des_mppic.$(OBJ_EXT) \
    $(DPO)check_des_rxns.$(OBJ_EXT) \
    $(DPO)check_des_thermo.$(OBJ_EXT) \
    $(DPO)des_allocate_arrays.$(OBJ_EXT) \
    $(DPO)des_bc_mod.$(OBJ_EXT) \
    $(DPO)des_check_particle.$(OBJ_EXT) \
    $(DPO)des_cluster_identification.$(OBJ_EXT) \
    $(DPO)des_cluster_mod.$(OBJ_EXT) \
    $(DPO)des_functions.$(OBJ_EXT) \
    $(DPO)des_granular_temperature.$(OBJ_EXT) \
    $(DPO)desgrid_mod.$(OBJ_EXT) \
    $(DPO)des_ic_mod.$(OBJ_EXT) \
    $(DPO)des_init_arrays.$(OBJ_EXT) \
    $(DPO)des_init_bc.$(OBJ_EXT) \
    $(DPO)des_init_namelist.$(OBJ_EXT) \
    $(DPO)des_linked_list_data_mod.$(OBJ_EXT) \
    $(DPO)des_linked_list_funcs_mod.$(OBJ_EXT) \
    $(DPO)des_mass_inlet.$(OBJ_EXT) \
    $(DPO)desmpi_mod.$(OBJ_EXT) \
    $(DPO)desmpi_wrapper_mod.$(OBJ_EXT) \
    $(DPO)des_physical_prop.$(OBJ_EXT) \
    $(DPO)des_reaction_model.$(OBJ_EXT) \
    $(DPO)des_rrates0.$(OBJ_EXT) \
    $(DPO)des_rxns_mod.$(OBJ_EXT) \
    $(DPO)des_set_ic.$(OBJ_EXT) \
    $(DPO)des_stl_functions_mod.$(OBJ_EXT) \
    $(DPO)des_thermo_cond.$(OBJ_EXT) \
    $(DPO)des_thermo_conv.$(OBJ_EXT) \
    $(DPO)des_thermo_mod.$(OBJ_EXT) \
    $(DPO)des_thermo_newvalues.$(OBJ_EXT) \
    $(DPO)des_thermo_rad.$(OBJ_EXT) \
    $(DPO)des_time_march.$(OBJ_EXT) \
    $(DPO)des_wallbc_preprocessing.$(OBJ_EXT) \
    $(DPO)discretelement_mod.$(OBJ_EXT) \
    $(DPO)drag_fgs.$(OBJ_EXT) \
    $(DPO)gas_drag.$(OBJ_EXT) \
    $(DPO)generate_particle_config.$(OBJ_EXT) \
    $(DPO)interpolation_mod.$(OBJ_EXT) \
    $(DPO)make_arrays_des.$(OBJ_EXT) \
    $(DPO)mfix_pic_mod.$(OBJ_EXT) \
    $(DPO)mppic_routines.$(OBJ_EXT) \
    $(DPO)mppic_wallbc_mod.$(OBJ_EXT) \
    $(DPO)neighbour.$(OBJ_EXT) \
    $(DPO)nsquare.$(OBJ_EXT) \
    $(DPO)particles_in_cell.$(OBJ_EXT) \
    $(DPO)randomno_mod.$(OBJ_EXT) \
    $(DPO)read_des_restart.$(OBJ_EXT) \
    $(DPO)sendrecvnode_mod.$(OBJ_EXT) \
    $(DPO)softspring_funcs_cutcell_mod.$(OBJ_EXT) \
    $(DPO)solid_drag.$(OBJ_EXT) \
    $(DPO)usr0_des.$(OBJ_EXT) \
    $(DPO)usr1_des.$(OBJ_EXT) \
    $(DPO)usr2_des.$(OBJ_EXT) \
    $(DPO)usr3_des.$(OBJ_EXT) \
    $(DPO)usr4_des.$(OBJ_EXT) \
    $(DPO)usr_rates_des.$(OBJ_EXT) \
    $(DPO)walledgecontact.$(OBJ_EXT) \
    $(DPO)wallfacecontact.$(OBJ_EXT) \
    $(DPO)wallnodecontact.$(OBJ_EXT) \
    $(DPO)write_des_data.$(OBJ_EXT) \
    $(DPO)write_des_restart.$(OBJ_EXT) \
    $(DPO)compar_mod.$(OBJ_EXT) \
    $(DPO)dbg_util_mod.$(OBJ_EXT) \
    $(DPO)debug_mod.$(OBJ_EXT) \
    $(DPO)gridmap_mod.$(OBJ_EXT) \
    $(DPO)mpi_mod.$(OBJ_EXT) \
    $(DPO)mpi_utility_mod.$(OBJ_EXT) \
    $(DPO)parallel_mpi_mod.$(OBJ_EXT) \
    $(DPO)sendrecv3_mod.$(OBJ_EXT) \
    $(DPO)sendrecv_mod.$(OBJ_EXT) \
    $(DPO)gaussj.$(OBJ_EXT) \
    $(DPO)odeint.$(OBJ_EXT) \
    $(DPO)rkck.$(OBJ_EXT) \
    $(DPO)rkqs.$(OBJ_EXT) \
    $(DPO)source_population_eq.$(OBJ_EXT) \
    $(DPO)usr_dqmom.$(OBJ_EXT) \
    $(DPO)adjust_eps_ghd.$(OBJ_EXT) \
    $(DPO)bulk_viscosity.$(OBJ_EXT) \
    $(DPO)calc_d_ghd.$(OBJ_EXT) \
    $(DPO)calc_external_forces.$(OBJ_EXT) \
    $(DPO)calc_nflux.$(OBJ_EXT) \
    $(DPO)chi_ij_GHD.$(OBJ_EXT) \
    $(DPO)cooling_rate.$(OBJ_EXT) \
    $(DPO)cooling_rate_tc.$(OBJ_EXT) \
    $(DPO)dufour_coeff.$(OBJ_EXT) \
    $(DPO)ghd.$(OBJ_EXT) \
    $(DPO)ghdmassflux.$(OBJ_EXT) \
    $(DPO)ghdtheory_mod.$(OBJ_EXT) \
    $(DPO)mass_mobility.$(OBJ_EXT) \
    $(DPO)ordinary_diff.$(OBJ_EXT) \
    $(DPO)partial_elim_ghd.$(OBJ_EXT) \
    $(DPO)pressure.$(OBJ_EXT) \
    $(DPO)shear_viscosity.$(OBJ_EXT) \
    $(DPO)source_ghd_granular_energy.$(OBJ_EXT) \
    $(DPO)thermal_conductivity.$(OBJ_EXT) \
    $(DPO)thermal_diffusivity.$(OBJ_EXT) \
    $(DPO)thermal_mobility.$(OBJ_EXT) \
    $(DPO)transport_coeff_ghd.$(OBJ_EXT) \
    $(DPO)qmomk_allocate_arrays.$(OBJ_EXT) \
    $(DPO)qmomk_bc_mod.$(OBJ_EXT) \
    $(DPO)qmomk_collision_mod.$(OBJ_EXT) \
    $(DPO)qmomk_fluxes_mod.$(OBJ_EXT) \
    $(DPO)qmomk_gas_drag.$(OBJ_EXT) \
    $(DPO)qmom_kinetic_equation_mod.$(OBJ_EXT) \
    $(DPO)qmomk_init_bc.$(OBJ_EXT) \
    $(DPO)qmomk_initial_conditions.$(OBJ_EXT) \
    $(DPO)qmomk_init_namelist.$(OBJ_EXT) \
    $(DPO)qmomk_make_arrays.$(OBJ_EXT) \
    $(DPO)qmomk_parameters_mod.$(OBJ_EXT) \
    $(DPO)qmomk_quadrature_mod.$(OBJ_EXT) \
    $(DPO)qmomk_read_restart.$(OBJ_EXT) \
    $(DPO)qmomk_set_bc.$(OBJ_EXT) \
    $(DPO)qmomk_time_march.$(OBJ_EXT) \
    $(DPO)qmomk_tools_mod.$(OBJ_EXT) \
    $(DPO)qmomk_write_restart.$(OBJ_EXT) \
    $(DPO)get_values.$(OBJ_EXT) \
    $(DPO)readTherm.$(OBJ_EXT) \
  -o $(EXEC_FILE) $(LIB_FLAGS)
  
$(DPO)blas90.a : $(DPO)BLAS.o
	ar cr $(DPO)blas90.a $(DPO)BLAS.o
$(DPO)BLAS.o : BLAS.F
	$(FORTRAN_CMD) $(FORT_FLAGS) BLAS.F -o $(DPO)BLAS.o
$(DPO)dgtsv90.a : $(DPO)DGTSV.o
	ar cr $(DPO)dgtsv90.a $(DPO)DGTSV.o
$(DPO)DGTSV.o : DGTSV.F
	$(FORTRAN_CMD) $(FORT_FLAGS) DGTSV.F -o $(DPO)DGTSV.o
$(DPO)odepack.a : $(DPO)ODEPACK.o
	ar cr $(DPO)odepack.a $(DPO)ODEPACK.o
$(DPO)ODEPACK.o : ODEPACK.F
	$(FORTRAN_CMD) $(FORT_FLAGS3) ODEPACK.F -o $(DPO)ODEPACK.o
$(DPO)ambm.mod : ambm_mod.f \
            $(DPO)compar.mod \
            $(DPO)funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ambm_mod.f  -o $(DPO)ambm_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)bc.mod : bc_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) bc_mod.f  -o $(DPO)bc_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)boundfunijk3.mod : boundfunijk3_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk3_mod.f  -o $(DPO)boundfunijk3_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)boundfunijk.mod : boundfunijk_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk_mod.f  -o $(DPO)boundfunijk_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)cdist.mod : cdist_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) cdist_mod.f  -o $(DPO)cdist_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)check.mod : check_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_mod.f  -o $(DPO)check_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)chischeme.mod : chischeme_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) chischeme_mod.f  -o $(DPO)chischeme_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)coeff.mod : coeff_mod.f \
            $(DPO)param.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)visc_g.mod \
            $(DPO)discretelement.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)param1.mod \
            $(DPO)mms.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) coeff_mod.f  -o $(DPO)coeff_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)constant.mod : constant_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) constant_mod.f  -o $(DPO)constant_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)cont.mod : cont_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) cont_mod.f  -o $(DPO)cont_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)corner.mod : corner_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) corner_mod.f  -o $(DPO)corner_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)dbg.mod : dbg_mod.f \
            $(DPO)param1.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)mflux.mod \
            $(DPO)parallel.mod \
            $(DPO)param.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)sendrecv.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) dbg_mod.f  -o $(DPO)dbg_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)drag.mod : drag_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) drag_mod.f  -o $(DPO)drag_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)energy.mod : energy_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) energy_mod.f  -o $(DPO)energy_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)error_manager.mod : error_manager_mod.f \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)param1.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) error_manager_mod.f  -o $(DPO)error_manager_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)fldvar.mod : fldvar_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) fldvar_mod.f  -o $(DPO)fldvar_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)function.mod : function_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) function_mod.f  -o $(DPO)function_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)funits.mod : funits_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) funits_mod.f  -o $(DPO)funits_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)geometry.mod : geometry_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) geometry_mod.f  -o $(DPO)geometry_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)ic.mod : ic_mod.f \
            $(DPO)param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ic_mod.f  -o $(DPO)ic_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)indices.mod : indices_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) indices_mod.f  -o $(DPO)indices_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)is.mod : is_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) is_mod.f  -o $(DPO)is_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory.mod : kintheory_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_mod.f  -o $(DPO)kintheory_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)leqsol.mod : leqsol_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) leqsol_mod.f  -o $(DPO)leqsol_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)machine.mod : machine_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) machine_mod.f  -o $(DPO)machine_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)matrix.mod : matrix_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) matrix_mod.f  -o $(DPO)matrix_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)mfix_netcdf.mod : mfix_netcdf_mod.f \
            MFIX_netcdf_constants.fi                                     \
            MFIX_netcdf_overloads.fi                                     \
            MFIX_netcdf_variables.fi                                     \
            MFIX_netcdf_misc.fi                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) mfix_netcdf_mod.f  -o $(DPO)mfix_netcdf_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)mflux.mod : mflux_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) mflux_mod.f  -o $(DPO)mflux_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)mms.mod : mms_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)fldvar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) mms_mod.f  -o $(DPO)mms_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)output.mod : output_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) output_mod.f  -o $(DPO)output_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)parallel.mod : parallel_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parallel_mod.f  -o $(DPO)parallel_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)param1.mod : param1_mod.f \
            $(DPO)param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) param1_mod.f  -o $(DPO)param1_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)param.mod : param_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) param_mod.f  -o $(DPO)param_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)parse.mod : parse_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)rxn_com.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_mod.f  -o $(DPO)parse_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)pgcor.mod : pgcor_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pgcor_mod.f  -o $(DPO)pgcor_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)physprop.mod : physprop_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) physprop_mod.f  -o $(DPO)physprop_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)pscor.mod : pscor_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pscor_mod.f  -o $(DPO)pscor_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)ps.mod : ps_mod.f \
            $(DPO)param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ps_mod.f  -o $(DPO)ps_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)residual.mod : residual_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) residual_mod.f  -o $(DPO)residual_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)run.mod : run_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) run_mod.f  -o $(DPO)run_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)rxn_com.mod : rxn_com_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)compar.mod \
            $(DPO)funits.mod \
            mfix_directory_path.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) rxn_com_mod.f  -o $(DPO)rxn_com_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)rxns.mod : rxns_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)rxn_com.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) rxns_mod.f  -o $(DPO)rxns_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)scalars.mod : scalars_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scalars_mod.f  -o $(DPO)scalars_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)scales.mod : scales_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scales_mod.f  -o $(DPO)scales_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_g.mod : tau_g_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_g_mod.f  -o $(DPO)tau_g_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_s.mod : tau_s_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_s_mod.f  -o $(DPO)tau_s_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)time_cpu.mod : time_cpu_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) time_cpu_mod.f  -o $(DPO)time_cpu_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)tmp_array1.mod : tmp_array1_mod.f \
            $(DPO)compar.mod \
            $(DPO)funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array1_mod.f  -o $(DPO)tmp_array1_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)tmp_array.mod : tmp_array_mod.f \
            $(DPO)compar.mod \
            $(DPO)funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array_mod.f  -o $(DPO)tmp_array_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)toleranc.mod : toleranc_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) toleranc_mod.f  -o $(DPO)toleranc_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)trace.mod : trace_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) trace_mod.f  -o $(DPO)trace_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)turb.mod : turb_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) turb_mod.f  -o $(DPO)turb_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)ur_facs.mod : ur_facs_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ur_facs_mod.f  -o $(DPO)ur_facs_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)usr.mod : usr_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_mod.f  -o $(DPO)usr_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)visc_g.mod : visc_g_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_g_mod.f  -o $(DPO)visc_g_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)visc_s.mod : visc_s_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_s_mod.f  -o $(DPO)visc_s_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)vshear.mod : vshear_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) vshear_mod.f  -o $(DPO)vshear_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)xsi_array.mod : xsi_array_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) xsi_array_mod.f  -o $(DPO)xsi_array_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)cutcell.mod : ./cartesian_grid/cutcell_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)progress_bar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cutcell_mod.f  -o $(DPO)cutcell_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)dashboard.mod : ./cartesian_grid/dashboard_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/dashboard_mod.f  -o $(DPO)dashboard_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)polygon.mod : ./cartesian_grid/polygon_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/polygon_mod.f  -o $(DPO)polygon_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)progress_bar.mod : ./cartesian_grid/progress_bar_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/progress_bar_mod.f  -o $(DPO)progress_bar_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)quadric.mod : ./cartesian_grid/quadric_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/quadric_mod.f  -o $(DPO)quadric_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)stl.mod : ./cartesian_grid/stl_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/stl_mod.f  -o $(DPO)stl_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)vtk.mod : ./cartesian_grid/vtk_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/vtk_mod.f  -o $(DPO)vtk_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)stiff_chem_dbg.mod : ./chem/stiff_chem_dbg_mod.f \
            $(DPO)compar.mod \
            $(DPO)stiff_chem_stats.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)rxns.mod \
            $(DPO)indices.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/stiff_chem_dbg_mod.f  -o $(DPO)stiff_chem_dbg_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)stiff_chem_maps.mod : ./chem/stiff_chem_maps_mod.f \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod \
            $(DPO)run.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/stiff_chem_maps_mod.f  -o $(DPO)stiff_chem_maps_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)stiff_chem.mod : ./chem/stiff_chem_mod.f \
            $(DPO)stiff_chem_maps.mod \
            $(DPO)funits.mod \
            $(DPO)output.mod \
            $(DPO)run.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)stiff_chem_dbg.mod \
            $(DPO)stiff_chem_stats.mod \
            $(DPO)rxns.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/stiff_chem_mod.f  -o $(DPO)stiff_chem_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)stiff_chem_stats.mod : ./chem/stiff_chem_stats_mod.f \
            $(DPO)compar.mod \
            $(DPO)output.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/stiff_chem_stats_mod.f  -o $(DPO)stiff_chem_stats_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)des_bc.mod : ./des/des_bc_mod.f \
            $(DPO)param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_bc_mod.f  -o $(DPO)des_bc_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)des_cluster.mod : ./des/des_cluster_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)compar.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)discretelement.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)desmpi_wrapper.mod \
            $(DPO)desmpi.mod \
            $(DPO)parallel.mod \
            $(DPO)sendrecv.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_cluster_mod.f  -o $(DPO)des_cluster_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)desgrid.mod : ./des/desgrid_mod.f \
            $(DPO)param1.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)constant.mod \
            $(DPO)desmpi_wrapper.mod \
            $(DPO)des_thermo.mod \
            des/desgrid_functions.inc                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desgrid_mod.f  -o $(DPO)desgrid_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)des_ic.mod : ./des/des_ic_mod.f \
            $(DPO)param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_ic_mod.f  -o $(DPO)des_ic_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)des_linked_list_data.mod : ./des/des_linked_list_data_mod.f \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_linked_list_data_mod.f  -o $(DPO)des_linked_list_data_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)des_linked_list_funcs.mod : ./des/des_linked_list_funcs_mod.f \
            $(DPO)des_linked_list_data.mod \
            $(DPO)discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_linked_list_funcs_mod.f  -o $(DPO)des_linked_list_funcs_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)desmpi.mod : ./des/desmpi_mod.f \
            $(DPO)parallel_mpi.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)discretelement.mod \
            $(DPO)desgrid.mod \
            $(DPO)compar.mod \
            $(DPO)physprop.mod \
            $(DPO)sendrecv.mod \
            $(DPO)des_bc.mod \
            $(DPO)desmpi_wrapper.mod \
            $(DPO)sendrecvnode.mod \
            $(DPO)mfix_pic.mod \
            des/desgrid_functions.inc                                    \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desmpi_mod.f  -o $(DPO)desmpi_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)desmpi_wrapper.mod : ./des/desmpi_wrapper_mod.f \
            $(DPO)parallel_mpi.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)compar.mod \
            $(DPO)funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desmpi_wrapper_mod.f  -o $(DPO)desmpi_wrapper_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)des_rxns.mod : ./des/des_rxns_mod.f \
            $(DPO)param.mod \
            $(DPO)rxn_com.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_rxns_mod.f  -o $(DPO)des_rxns_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)des_stl_functions.mod : ./des/des_stl_functions_mod.f \
            $(DPO)param1.mod \
            $(DPO)funits.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cutcell.mod \
            $(DPO)stl.mod \
            $(DPO)indices.mod \
            $(DPO)geometry.mod \
            $(DPO)bc.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)param.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)polygon.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_stl_functions_mod.f  -o $(DPO)des_stl_functions_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)des_thermo.mod : ./des/des_thermo_mod.f \
            $(DPO)param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_mod.f  -o $(DPO)des_thermo_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)discretelement.mod : ./des/discretelement_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/discretelement_mod.f  -o $(DPO)discretelement_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)interpolation.mod : ./des/interpolation_mod.f \
            $(DPO)constant.mod \
            $(DPO)discretelement.mod \
            $(DPO)geometry.mod \
            $(DPO)param1.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)param.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/interpolation_mod.f  -o $(DPO)interpolation_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)mfix_pic.mod : ./des/mfix_pic_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mfix_pic_mod.f  -o $(DPO)mfix_pic_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)mppic_wallbc.mod : ./des/mppic_wallbc_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)randomno.mod \
            $(DPO)cutcell.mod \
            $(DPO)fldvar.mod \
            $(DPO)mfix_pic.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mppic_wallbc_mod.f  -o $(DPO)mppic_wallbc_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)randomno.mod : ./des/randomno_mod.f \
            $(DPO)constant.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/randomno_mod.f  -o $(DPO)randomno_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)sendrecvnode.mod : ./des/sendrecvnode_mod.f \
            $(DPO)parallel_mpi.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)discretelement.mod \
            $(DPO)compar.mod \
            $(DPO)physprop.mod \
            $(DPO)sendrecv.mod \
            $(DPO)desmpi_wrapper.mod \
            $(DPO)desgrid.mod \
            function.inc                                                 \
            des/desgrid_functions.inc                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/sendrecvnode_mod.f  -o $(DPO)sendrecvnode_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)softspring_funcs_cutcell.mod : ./des/softspring_funcs_cutcell_mod.f \
            $(DPO)run.mod \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)cutcell.mod \
            $(DPO)funits.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)parallel.mod \
            $(DPO)stl.mod \
            $(DPO)des_stl_functions.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/softspring_funcs_cutcell_mod.f  -o $(DPO)softspring_funcs_cutcell_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)compar.mod : ./dmp_modules/compar_mod.f \
            $(DPO)mpi.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/compar_mod.f  -o $(DPO)compar_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)dbg_util.mod : ./dmp_modules/dbg_util_mod.f \
            $(DPO)compar.mod \
            $(DPO)geometry.mod \
            $(DPO)parallel_mpi.mod \
            $(DPO)indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/dbg_util_mod.f  -o $(DPO)dbg_util_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)debug.mod : ./dmp_modules/debug_mod.f \
            $(DPO)funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/debug_mod.f  -o $(DPO)debug_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)gridmap.mod : ./dmp_modules/gridmap_mod.f \
            $(DPO)mpi_utility.mod \
            $(DPO)parallel_mpi.mod \
            $(DPO)geometry.mod \
            $(DPO)sendrecv.mod \
            $(DPO)compar.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)error_manager.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/gridmap_mod.f  -o $(DPO)gridmap_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)mpi.mod : ./dmp_modules/mpi_mod.f \
            mpif.h                                                      
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_mod.f  -o $(DPO)mpi_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)mpi_utility.mod : ./dmp_modules/mpi_utility_mod.f \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)parallel_mpi.mod \
            $(DPO)debug.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_utility_mod.f  -o $(DPO)mpi_utility_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)parallel_mpi.mod : ./dmp_modules/parallel_mpi_mod.f \
            $(DPO)geometry.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/parallel_mpi_mod.f  -o $(DPO)parallel_mpi_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)sendrecv3.mod : ./dmp_modules/sendrecv3_mod.f \
            $(DPO)parallel_mpi.mod \
            $(DPO)debug.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)mpi.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv3_mod.f  -o $(DPO)sendrecv3_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)sendrecv.mod : ./dmp_modules/sendrecv_mod.f \
            $(DPO)parallel_mpi.mod \
            $(DPO)debug.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)mpi.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv_mod.f  -o $(DPO)sendrecv_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)ghdtheory.mod : ./GhdTheory/ghdtheory_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ghdtheory_mod.f  -o $(DPO)ghdtheory_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_bc.mod : ./qmomk/qmomk_bc_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)bc.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)qmomk_quadrature.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_bc_mod.f  -o $(DPO)qmomk_bc_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_collision.mod : ./qmomk/qmomk_collision_mod.f \
            $(DPO)constant.mod \
            $(DPO)qmomk_parameters.mod \
            $(DPO)qmomk_quadrature.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_collision_mod.f  -o $(DPO)qmomk_collision_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_fluxes.mod : ./qmomk/qmomk_fluxes_mod.f \
            $(DPO)qmomk_parameters.mod \
            $(DPO)qmomk_collision.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_fluxes_mod.f  -o $(DPO)qmomk_fluxes_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)qmom_kinetic_equation.mod : ./qmomk/qmom_kinetic_equation_mod.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)qmomk_parameters.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmom_kinetic_equation_mod.f  -o $(DPO)qmom_kinetic_equation_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_parameters.mod : ./qmomk/qmomk_parameters_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_parameters_mod.f  -o $(DPO)qmomk_parameters_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_quadrature.mod : ./qmomk/qmomk_quadrature_mod.f \
            $(DPO)qmomk_tools.mod \
            $(DPO)qmomk_parameters.mod \
            $(DPO)constant.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_quadrature_mod.f  -o $(DPO)qmomk_quadrature_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_tools.mod : ./qmomk/qmomk_tools_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_tools_mod.f  -o $(DPO)qmomk_tools_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)accum_resid.$(OBJ_EXT) : accum_resid.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)residual.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) accum_resid.f  -o $(DPO)accum_resid.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_u_g.$(OBJ_EXT) : adjust_a_u_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)usr.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_u_g.f  -o $(DPO)adjust_a_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_u_s.$(OBJ_EXT) : adjust_a_u_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_u_s.f  -o $(DPO)adjust_a_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_v_g.$(OBJ_EXT) : adjust_a_v_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_v_g.f  -o $(DPO)adjust_a_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_v_s.$(OBJ_EXT) : adjust_a_v_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_v_s.f  -o $(DPO)adjust_a_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_w_g.$(OBJ_EXT) : adjust_a_w_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_w_g.f  -o $(DPO)adjust_a_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_w_s.$(OBJ_EXT) : adjust_a_w_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_w_s.f  -o $(DPO)adjust_a_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_dt.$(OBJ_EXT) : adjust_dt.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_dt.f  -o $(DPO)adjust_dt.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_eps.$(OBJ_EXT) : adjust_eps.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_eps.f  -o $(DPO)adjust_eps.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_leq.$(OBJ_EXT) : adjust_leq.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)leqsol.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_leq.f  -o $(DPO)adjust_leq.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_rop.$(OBJ_EXT) : adjust_rop.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_rop.f  -o $(DPO)adjust_rop.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_theta.$(OBJ_EXT) : adjust_theta.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_theta.f  -o $(DPO)adjust_theta.$(OBJ_EXT) -module $(DPO)
$(DPO)allocate_arrays.$(OBJ_EXT) : allocate_arrays.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)ambm.mod \
            $(DPO)cont.mod \
            $(DPO)drag.mod \
            $(DPO)energy.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)pgcor.mod \
            $(DPO)physprop.mod \
            $(DPO)pscor.mod \
            $(DPO)residual.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)scalars.mod \
            $(DPO)turb.mod \
            $(DPO)tau_g.mod \
            $(DPO)tau_s.mod \
            $(DPO)tmp_array.mod \
            $(DPO)tmp_array1.mod \
            $(DPO)trace.mod \
            $(DPO)visc_g.mod \
            $(DPO)visc_s.mod \
            $(DPO)xsi_array.mod \
            $(DPO)vshear.mod \
            $(DPO)mflux.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)kintheory.mod \
            $(DPO)cdist.mod \
            $(DPO)des_rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) allocate_arrays.f  -o $(DPO)allocate_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)bc_phi.$(OBJ_EXT) : bc_phi.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)cutcell.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) bc_phi.f  -o $(DPO)bc_phi.$(OBJ_EXT) -module $(DPO)
$(DPO)bc_theta.$(OBJ_EXT) : bc_theta.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)geometry.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)turb.mod \
            $(DPO)rxns.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) bc_theta.f  -o $(DPO)bc_theta.$(OBJ_EXT) -module $(DPO)
$(DPO)b_m_p_star.$(OBJ_EXT) : b_m_p_star.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)rxns.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) b_m_p_star.f  -o $(DPO)b_m_p_star.$(OBJ_EXT) -module $(DPO)
$(DPO)bound_x.$(OBJ_EXT) : bound_x.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) bound_x.f  -o $(DPO)bound_x.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_cell.$(OBJ_EXT) : calc_cell.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_cell.f  -o $(DPO)calc_cell.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_coeff.$(OBJ_EXT) : calc_coeff.f \
            $(DPO)param1.mod \
            $(DPO)ur_facs.mod \
            $(DPO)rxns.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)des_rxns.mod \
            $(DPO)visc_g.mod \
            $(DPO)visc_s.mod \
            $(DPO)tau_g.mod \
            $(DPO)tau_s.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_coeff.f  -o $(DPO)calc_coeff.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_d.$(OBJ_EXT) : calc_d.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)scales.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)cutcell.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)discretelement.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_d.f  -o $(DPO)calc_d.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_dif_g.$(OBJ_EXT) : calc_dif_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)constant.mod \
            $(DPO)scales.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)run.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_dif_g.f  -o $(DPO)calc_dif_g.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_dif_s.$(OBJ_EXT) : calc_dif_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)run.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_dif_s.f  -o $(DPO)calc_dif_s.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_drag.$(OBJ_EXT) : calc_drag.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)drag.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)qmom_kinetic_equation.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_drag.f  -o $(DPO)calc_drag.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_e.$(OBJ_EXT) : calc_e.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_e.f  -o $(DPO)calc_e.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_gama.$(OBJ_EXT) : calc_gama.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)energy.mod \
            $(DPO)rxns.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_gama.f  -o $(DPO)calc_gama.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_grbdry.$(OBJ_EXT) : calc_grbdry.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)turb.mod \
            $(DPO)visc_s.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)toleranc.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_grbdry.f  -o $(DPO)calc_grbdry.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_h.$(OBJ_EXT) : calc_h.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)constant.mod \
            $(DPO)discretelement.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_h.f  -o $(DPO)calc_h.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_k_cp.$(OBJ_EXT) : calc_k_cp.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)pscor.mod \
            $(DPO)geometry.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)visc_s.mod \
            $(DPO)trace.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_k_cp.f  -o $(DPO)calc_k_cp.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_k_g.$(OBJ_EXT) : calc_k_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod \
            $(DPO)run.mod \
            $(DPO)sendrecv.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_k_g.f  -o $(DPO)calc_k_g.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_k_s.$(OBJ_EXT) : calc_k_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)run.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_k_s.f  -o $(DPO)calc_k_s.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_mflux.$(OBJ_EXT) : calc_mflux.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)mflux.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_mflux.f  -o $(DPO)calc_mflux.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_mu_g.$(OBJ_EXT) : calc_mu_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)visc_s.mod \
            $(DPO)indices.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)compar.mod \
            $(DPO)drag.mod \
            $(DPO)run.mod \
            $(DPO)turb.mod \
            $(DPO)sendrecv.mod \
            $(DPO)mms.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_mu_g.f  -o $(DPO)calc_mu_g.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_mu_s.$(OBJ_EXT) : calc_mu_s.f \
            $(DPO)run.mod \
            $(DPO)vshear.mod \
            $(DPO)visc_s.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)geometry.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)mms.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)trace.mod \
            $(DPO)toleranc.mod \
            $(DPO)turb.mod \
            $(DPO)drag.mod \
            $(DPO)kintheory.mod \
            $(DPO)ur_facs.mod \
            $(DPO)cutcell.mod \
            $(DPO)parallel.mod \
            $(DPO)visc_g.mod \
            $(DPO)is.mod \
            $(DPO)sendrecv.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            s_pr1.inc                                                    \
            s_pr2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_mu_s.f  -o $(DPO)calc_mu_s.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_mw.$(OBJ_EXT) : calc_mw.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_mw.f  -o $(DPO)calc_mw.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_outflow.$(OBJ_EXT) : calc_outflow.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)bc.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_outflow.f  -o $(DPO)calc_outflow.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_p_star.$(OBJ_EXT) : calc_p_star.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)ur_facs.mod \
            $(DPO)residual.mod \
            $(DPO)compar.mod \
            $(DPO)run.mod \
            $(DPO)visc_s.mod \
            $(DPO)fldvar.mod \
            $(DPO)toleranc.mod \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_p_star.f  -o $(DPO)calc_p_star.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_resid.$(OBJ_EXT) : calc_resid.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)run.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)toleranc.mod \
            $(DPO)bc.mod \
            $(DPO)constant.mod \
            $(DPO)residual.mod \
            $(DPO)rxns.mod \
            $(DPO)mflux.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_resid.f  -o $(DPO)calc_resid.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_s_ddot_s.$(OBJ_EXT) : calc_s_ddot_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_s_ddot_s.f  -o $(DPO)calc_s_ddot_s.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_trd_g.$(OBJ_EXT) : calc_trd_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)bc.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_trd_g.f  -o $(DPO)calc_trd_g.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_trd_s.$(OBJ_EXT) : calc_trd_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)bc.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_trd_s.f  -o $(DPO)calc_trd_s.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_u_friction.$(OBJ_EXT) : calc_u_friction.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)bc.mod \
            $(DPO)run.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_u_friction.f  -o $(DPO)calc_u_friction.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_vol_fr.$(OBJ_EXT) : calc_vol_fr.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)visc_s.mod \
            $(DPO)constant.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)fldvar.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_vol_fr.f  -o $(DPO)calc_vol_fr.$(OBJ_EXT) -module $(DPO)
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_vol_fr.f  -o $(DPO)calc_vol_fr.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_xsi.$(OBJ_EXT) : calc_xsi.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)vshear.mod \
            $(DPO)chischeme.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            xsi1.inc                                                     \
            function.inc                                                 \
            xsi2.inc                                                    
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_xsi.f  -o $(DPO)calc_xsi.$(OBJ_EXT) -module $(DPO)
$(DPO)cal_d.$(OBJ_EXT) : cal_d.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_s.mod \
            $(DPO)bc.mod \
            $(DPO)vshear.mod \
            $(DPO)compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) cal_d.f  -o $(DPO)cal_d.$(OBJ_EXT) -module $(DPO)
$(DPO)check_ab_m.$(OBJ_EXT) : check_ab_m.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) check_ab_m.f  -o $(DPO)check_ab_m.$(OBJ_EXT) -module $(DPO)
$(DPO)check_convergence.$(OBJ_EXT) : check_convergence.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)residual.mod \
            $(DPO)toleranc.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)scalars.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_convergence.f  -o $(DPO)check_convergence.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_03.$(OBJ_EXT) : check_data_03.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)bc.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_data_03.f  -o $(DPO)check_data_03.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_06.$(OBJ_EXT) : check_data_06.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)ic.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)scalars.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)rxns.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) check_data_06.f  -o $(DPO)check_data_06.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_07.$(OBJ_EXT) : check_data_07.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)bc.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)scalars.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cutcell.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) check_data_07.f  -o $(DPO)check_data_07.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_08.$(OBJ_EXT) : check_data_08.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_data_08.f  -o $(DPO)check_data_08.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_09.$(OBJ_EXT) : check_data_09.f \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)fldvar.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parse.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_data_09.f  -o $(DPO)check_data_09.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_10.$(OBJ_EXT) : check_data_10.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_data_10.f  -o $(DPO)check_data_10.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_20.$(OBJ_EXT) : check_data_20.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)scalars.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) check_data_20.f  -o $(DPO)check_data_20.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_30.$(OBJ_EXT) : check_data_30.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)fldvar.mod \
            $(DPO)rxns.mod \
            $(DPO)visc_s.mod \
            $(DPO)visc_g.mod \
            $(DPO)geometry.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)discretelement.mod \
            $(DPO)mms.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) check_data_30.f  -o $(DPO)check_data_30.$(OBJ_EXT) -module $(DPO)
$(DPO)check_mass_balance.$(OBJ_EXT) : check_mass_balance.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)fldvar.mod \
            $(DPO)rxns.mod \
            $(DPO)geometry.mod \
            $(DPO)run.mod \
            $(DPO)bc.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)output.mod \
            $(DPO)check.mod \
            $(DPO)mflux.mod \
            $(DPO)xsi_array.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) check_mass_balance.f  -o $(DPO)check_mass_balance.$(OBJ_EXT) -module $(DPO)
$(DPO)check_one_axis.$(OBJ_EXT) : check_one_axis.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_one_axis.f  -o $(DPO)check_one_axis.$(OBJ_EXT) -module $(DPO)
$(DPO)check_plane.$(OBJ_EXT) : check_plane.f \
            $(DPO)funits.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_plane.f  -o $(DPO)check_plane.$(OBJ_EXT) -module $(DPO)
$(DPO)cn_extrapol.$(OBJ_EXT) : cn_extrapol.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)scalars.mod \
            $(DPO)trace.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) cn_extrapol.f  -o $(DPO)cn_extrapol.$(OBJ_EXT) -module $(DPO)
$(DPO)compare.$(OBJ_EXT) : compare.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) compare.f  -o $(DPO)compare.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_phi.$(OBJ_EXT) : conv_dif_phi.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)xsi_array.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)indices.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)toleranc.mod \
            $(DPO)cutcell.mod \
            $(DPO)sendrecv3.mod \
            $(DPO)tmp_array.mod \
            $(DPO)vshear.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)output.mod \
            $(DPO)is.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            function3.inc                                                \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_phi.f  -o $(DPO)conv_dif_phi.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_u_g.$(OBJ_EXT) : conv_dif_u_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)visc_g.mod \
            $(DPO)compar.mod \
            $(DPO)toleranc.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)mflux.mod \
            $(DPO)cutcell.mod \
            $(DPO)vshear.mod \
            $(DPO)xsi_array.mod \
            $(DPO)tmp_array.mod \
            $(DPO)sendrecv.mod \
            $(DPO)sendrecv3.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_u_g.f  -o $(DPO)conv_dif_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_u_s.$(OBJ_EXT) : conv_dif_u_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)visc_s.mod \
            $(DPO)compar.mod \
            $(DPO)toleranc.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)mflux.mod \
            $(DPO)cutcell.mod \
            $(DPO)xsi_array.mod \
            $(DPO)tmp_array.mod \
            $(DPO)sendrecv.mod \
            $(DPO)sendrecv3.mod \
            $(DPO)vshear.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_u_s.f  -o $(DPO)conv_dif_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_v_g.$(OBJ_EXT) : conv_dif_v_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)visc_g.mod \
            $(DPO)compar.mod \
            $(DPO)toleranc.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)mflux.mod \
            $(DPO)cutcell.mod \
            $(DPO)xsi_array.mod \
            $(DPO)vshear.mod \
            $(DPO)tmp_array.mod \
            $(DPO)sendrecv.mod \
            $(DPO)sendrecv3.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_v_g.f  -o $(DPO)conv_dif_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_v_s.$(OBJ_EXT) : conv_dif_v_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)visc_s.mod \
            $(DPO)compar.mod \
            $(DPO)toleranc.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)mflux.mod \
            $(DPO)cutcell.mod \
            $(DPO)xsi_array.mod \
            $(DPO)tmp_array.mod \
            $(DPO)sendrecv.mod \
            $(DPO)sendrecv3.mod \
            $(DPO)vshear.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_v_s.f  -o $(DPO)conv_dif_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_w_g.$(OBJ_EXT) : conv_dif_w_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)visc_g.mod \
            $(DPO)compar.mod \
            $(DPO)toleranc.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)mflux.mod \
            $(DPO)cutcell.mod \
            $(DPO)xsi_array.mod \
            $(DPO)tmp_array.mod \
            $(DPO)sendrecv.mod \
            $(DPO)sendrecv3.mod \
            $(DPO)vshear.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_w_g.f  -o $(DPO)conv_dif_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_w_s.$(OBJ_EXT) : conv_dif_w_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)visc_s.mod \
            $(DPO)compar.mod \
            $(DPO)toleranc.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)mflux.mod \
            $(DPO)cutcell.mod \
            $(DPO)xsi_array.mod \
            $(DPO)tmp_array.mod \
            $(DPO)sendrecv.mod \
            $(DPO)sendrecv3.mod \
            $(DPO)vshear.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_w_s.f  -o $(DPO)conv_dif_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_pp_g.$(OBJ_EXT) : conv_pp_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)pgcor.mod \
            $(DPO)compar.mod \
            $(DPO)mflux.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_pp_g.f  -o $(DPO)conv_pp_g.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_rop.$(OBJ_EXT) : conv_rop.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)mflux.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)xsi_array.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_rop.f  -o $(DPO)conv_rop.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_rop_g.$(OBJ_EXT) : conv_rop_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)pgcor.mod \
            $(DPO)xsi_array.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_rop_g.f  -o $(DPO)conv_rop_g.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_rop_s.$(OBJ_EXT) : conv_rop_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)xsi_array.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_rop_s.f  -o $(DPO)conv_rop_s.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_source_epp.$(OBJ_EXT) : conv_source_epp.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)xsi_array.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)rxns.mod \
            $(DPO)indices.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)vshear.mod \
            $(DPO)ps.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_source_epp.f  -o $(DPO)conv_source_epp.$(OBJ_EXT) -module $(DPO)
$(DPO)copy_a.$(OBJ_EXT) : copy_a.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)physprop.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) copy_a.f  -o $(DPO)copy_a.$(OBJ_EXT) -module $(DPO)
$(DPO)corner.$(OBJ_EXT) : corner.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)matrix.mod \
            $(DPO)corner.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) corner.f  -o $(DPO)corner.$(OBJ_EXT) -module $(DPO)
$(DPO)correct_0.$(OBJ_EXT) : correct_0.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)pgcor.mod \
            $(DPO)ur_facs.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)compar.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) correct_0.f  -o $(DPO)correct_0.$(OBJ_EXT) -module $(DPO)
$(DPO)correct_1.$(OBJ_EXT) : correct_1.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)geometry.mod \
            $(DPO)pscor.mod \
            $(DPO)ur_facs.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) correct_1.f  -o $(DPO)correct_1.$(OBJ_EXT) -module $(DPO)
$(DPO)dgtsl.$(OBJ_EXT) : dgtsl.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) dgtsl.f  -o $(DPO)dgtsl.$(OBJ_EXT) -module $(DPO)
$(DPO)dif_u_is.$(OBJ_EXT) : dif_u_is.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) dif_u_is.f  -o $(DPO)dif_u_is.$(OBJ_EXT) -module $(DPO)
$(DPO)dif_v_is.$(OBJ_EXT) : dif_v_is.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) dif_v_is.f  -o $(DPO)dif_v_is.$(OBJ_EXT) -module $(DPO)
$(DPO)dif_w_is.$(OBJ_EXT) : dif_w_is.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) dif_w_is.f  -o $(DPO)dif_w_is.$(OBJ_EXT) -module $(DPO)
$(DPO)discretize.$(OBJ_EXT) : discretize.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) discretize.f  -o $(DPO)discretize.$(OBJ_EXT) -module $(DPO)
$(DPO)display_resid.$(OBJ_EXT) : display_resid.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)residual.mod \
            $(DPO)fldvar.mod \
            $(DPO)compar.mod \
            $(DPO)geometry.mod \
            $(DPO)scalars.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) display_resid.f  -o $(DPO)display_resid.$(OBJ_EXT) -module $(DPO)
$(DPO)drag_gs.$(OBJ_EXT) : drag_gs.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod \
            $(DPO)drag.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)ur_facs.mod \
            $(DPO)funits.mod \
            $(DPO)mms.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) drag_gs.f  -o $(DPO)drag_gs.$(OBJ_EXT) -module $(DPO)
$(DPO)drag_ss.$(OBJ_EXT) : drag_ss.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)drag.mod \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) drag_ss.f  -o $(DPO)drag_ss.$(OBJ_EXT) -module $(DPO)
$(DPO)eosg.$(OBJ_EXT) : eosg.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)scales.mod \
            sc_p_g1.inc                                                  \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) eosg.f  -o $(DPO)eosg.$(OBJ_EXT) -module $(DPO)
$(DPO)eoss.$(OBJ_EXT) : eoss.f \
            $(DPO)physprop.mod \
            $(DPO)compar.mod \
            $(DPO)funits.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) eoss.f  -o $(DPO)eoss.$(OBJ_EXT) -module $(DPO)
$(DPO)equal.$(OBJ_EXT) : equal.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) equal.f  -o $(DPO)equal.$(OBJ_EXT) -module $(DPO)
$(DPO)error_routine.$(OBJ_EXT) : error_routine.f \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) error_routine.f  -o $(DPO)error_routine.$(OBJ_EXT) -module $(DPO)
$(DPO)exchange.$(OBJ_EXT) : exchange.f \
            $(DPO)coeff.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) exchange.f  -o $(DPO)exchange.$(OBJ_EXT) -module $(DPO)
$(DPO)exit.$(OBJ_EXT) : exit.f \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) exit.f  -o $(DPO)exit.$(OBJ_EXT) -module $(DPO)
$(DPO)flow_to_vel.$(OBJ_EXT) : flow_to_vel.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)discretelement.mod \
            $(DPO)bc.mod \
            $(DPO)error_manager.mod \
            $(DPO)scales.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)mfix_pic.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) flow_to_vel.f  -o $(DPO)flow_to_vel.$(OBJ_EXT) -module $(DPO)
$(DPO)g_0.$(OBJ_EXT) : g_0.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)visc_s.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) g_0.f  -o $(DPO)g_0.$(OBJ_EXT) -module $(DPO)
$(DPO)get_bc_area.$(OBJ_EXT) : get_bc_area.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)parallel.mod \
            $(DPO)indices.mod \
            $(DPO)sendrecv.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)cutcell.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) get_bc_area.f  -o $(DPO)get_bc_area.$(OBJ_EXT) -module $(DPO)
$(DPO)get_data.$(OBJ_EXT) : get_data.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)gridmap.mod \
            $(DPO)discretelement.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod \
            $(DPO)leqsol.mod \
            $(DPO)parallel.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cutcell.mod \
            $(DPO)dashboard.mod \
            $(DPO)visc_g.mod \
            $(DPO)constant.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_data.f  -o $(DPO)get_data.$(OBJ_EXT) -module $(DPO)
$(DPO)get_eq.$(OBJ_EXT) : get_eq.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_eq.f  -o $(DPO)get_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)get_flow_bc.$(OBJ_EXT) : get_flow_bc.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_flow_bc.f  -o $(DPO)get_flow_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)get_hloss.$(OBJ_EXT) : get_hloss.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)bc.mod \
            $(DPO)indices.mod \
            $(DPO)energy.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_hloss.f  -o $(DPO)get_hloss.$(OBJ_EXT) -module $(DPO)
$(DPO)get_is.$(OBJ_EXT) : get_is.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)is.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_is.f  -o $(DPO)get_is.$(OBJ_EXT) -module $(DPO)
$(DPO)get_philoss.$(OBJ_EXT) : get_philoss.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)bc.mod \
            $(DPO)indices.mod \
            $(DPO)energy.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) get_philoss.f  -o $(DPO)get_philoss.$(OBJ_EXT) -module $(DPO)
$(DPO)get_ps.$(OBJ_EXT) : get_ps.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)ps.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_ps.f  -o $(DPO)get_ps.$(OBJ_EXT) -module $(DPO)
$(DPO)get_smass.$(OBJ_EXT) : get_smass.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) get_smass.f  -o $(DPO)get_smass.$(OBJ_EXT) -module $(DPO)
$(DPO)get_stats.$(OBJ_EXT) : get_stats.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)residual.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) get_stats.f  -o $(DPO)get_stats.$(OBJ_EXT) -module $(DPO)
$(DPO)get_walls_bc.$(OBJ_EXT) : get_walls_bc.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_walls_bc.f  -o $(DPO)get_walls_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)in_bin_512.$(OBJ_EXT) : in_bin_512.f \
            $(DPO)machine.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) in_bin_512.f  -o $(DPO)in_bin_512.$(OBJ_EXT) -module $(DPO)
$(DPO)in_bin_512i.$(OBJ_EXT) : in_bin_512i.f \
            $(DPO)machine.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) in_bin_512i.f  -o $(DPO)in_bin_512i.$(OBJ_EXT) -module $(DPO)
$(DPO)init_ab_m.$(OBJ_EXT) : init_ab_m.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)parallel.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) init_ab_m.f  -o $(DPO)init_ab_m.$(OBJ_EXT) -module $(DPO)
$(DPO)init_fvars.$(OBJ_EXT) : init_fvars.f \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)rxns.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) init_fvars.f  -o $(DPO)init_fvars.$(OBJ_EXT) -module $(DPO)
$(DPO)init_namelist.$(OBJ_EXT) : init_namelist.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)ic.mod \
            $(DPO)bc.mod \
            $(DPO)ps.mod \
            $(DPO)fldvar.mod \
            $(DPO)constant.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)toleranc.mod \
            $(DPO)scales.mod \
            $(DPO)ur_facs.mod \
            $(DPO)leqsol.mod \
            $(DPO)residual.mod \
            $(DPO)rxns.mod \
            $(DPO)scalars.mod \
            $(DPO)compar.mod \
            $(DPO)parallel.mod \
            $(DPO)cdist.mod \
            $(DPO)stiff_chem.mod \
            namelist.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) init_namelist.f  -o $(DPO)init_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)init_resid.$(OBJ_EXT) : init_resid.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)residual.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) init_resid.f  -o $(DPO)init_resid.$(OBJ_EXT) -module $(DPO)
$(DPO)iterate.$(OBJ_EXT) : iterate.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)time_cpu.mod \
            $(DPO)pscor.mod \
            $(DPO)leqsol.mod \
            $(DPO)visc_g.mod \
            $(DPO)pgcor.mod \
            $(DPO)cont.mod \
            $(DPO)scalars.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)discretelement.mod \
            $(DPO)residual.mod \
            $(DPO)cutcell.mod \
            $(DPO)vtk.mod \
            $(DPO)dashboard.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)stiff_chem.mod \
            $(DPO)rxns.mod \
            $(DPO)mms.mod \
            $(DPO)bc.mod \
            $(DPO)constant.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) iterate.f  -o $(DPO)iterate.$(OBJ_EXT) -module $(DPO)
$(DPO)k_epsilon_prop.$(OBJ_EXT) : k_epsilon_prop.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)drag.mod \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)visc_s.mod \
            $(DPO)trace.mod \
            $(DPO)indices.mod \
            $(DPO)constant.mod \
            $(DPO)vshear.mod \
            $(DPO)turb.mod \
            $(DPO)toleranc.mod \
            $(DPO)compar.mod \
            $(DPO)tau_g.mod \
            $(DPO)sendrecv.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) k_epsilon_prop.f  -o $(DPO)k_epsilon_prop.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory_drag_ss.$(OBJ_EXT) : kintheory_drag_ss.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)drag.mod \
            $(DPO)kintheory.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_drag_ss.f  -o $(DPO)kintheory_drag_ss.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory_energy_dissipation_ss.$(OBJ_EXT) : kintheory_energy_dissipation_ss.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)kintheory.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_energy_dissipation_ss.f  -o $(DPO)kintheory_energy_dissipation_ss.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory_u_s.$(OBJ_EXT) : kintheory_u_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)kintheory.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_u_s.f  -o $(DPO)kintheory_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory_v_s.$(OBJ_EXT) : kintheory_v_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)kintheory.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_v_s.f  -o $(DPO)kintheory_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory_w_s.$(OBJ_EXT) : kintheory_w_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)kintheory.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_w_s.f  -o $(DPO)kintheory_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)leq_bicgs.$(OBJ_EXT) : leq_bicgs.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)leqsol.mod \
            $(DPO)funits.mod \
            $(DPO)parallel.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)cutcell.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) leq_bicgs.f  -o $(DPO)leq_bicgs.$(OBJ_EXT) -module $(DPO)
$(DPO)leq_bicgst.$(OBJ_EXT) : leq_bicgst.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)leqsol.mod \
            $(DPO)funits.mod \
            $(DPO)parallel.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) leq_bicgst.f  -o $(DPO)leq_bicgst.$(OBJ_EXT) -module $(DPO)
$(DPO)leq_cg.$(OBJ_EXT) : leq_cg.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)leqsol.mod \
            $(DPO)funits.mod \
            $(DPO)parallel.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) leq_cg.f  -o $(DPO)leq_cg.$(OBJ_EXT) -module $(DPO)
$(DPO)leq_gmres.$(OBJ_EXT) : leq_gmres.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)debug.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)parallel.mod \
            $(DPO)funits.mod \
            $(DPO)gridmap.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) leq_gmres.f  -o $(DPO)leq_gmres.$(OBJ_EXT) -module $(DPO)
$(DPO)leq_sor.$(OBJ_EXT) : leq_sor.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)leqsol.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) leq_sor.f  -o $(DPO)leq_sor.$(OBJ_EXT) -module $(DPO)
$(DPO)line_too_big.$(OBJ_EXT) : line_too_big.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) line_too_big.f  -o $(DPO)line_too_big.$(OBJ_EXT) -module $(DPO)
$(DPO)location_check.$(OBJ_EXT) : location_check.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) location_check.f  -o $(DPO)location_check.$(OBJ_EXT) -module $(DPO)
$(DPO)location.$(OBJ_EXT) : location.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) location.f  -o $(DPO)location.$(OBJ_EXT) -module $(DPO)
$(DPO)machine.$(OBJ_EXT) : machine.f \
            $(DPO)machine.mod \
            $(DPO)param.mod \
            $(DPO)run.mod \
            $(DPO)funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) machine.f  -o $(DPO)machine.$(OBJ_EXT) -module $(DPO)
$(DPO)make_upper_case.$(OBJ_EXT) : make_upper_case.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) make_upper_case.f  -o $(DPO)make_upper_case.$(OBJ_EXT) -module $(DPO)
$(DPO)mark_phase_4_cor.$(OBJ_EXT) : mark_phase_4_cor.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod \
            $(DPO)visc_s.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) mark_phase_4_cor.f  -o $(DPO)mark_phase_4_cor.$(OBJ_EXT) -module $(DPO)
	$(FORTRAN_CMD) $(FORT_FLAGS) mark_phase_4_cor.f  -o $(DPO)mark_phase_4_cor.$(OBJ_EXT) -module $(DPO)
$(DPO)mfix.$(OBJ_EXT) : mfix.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)time_cpu.mod \
            $(DPO)funits.mod \
            $(DPO)output.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)parallel_mpi.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cdist.mod \
            $(DPO)mfix_netcdf.mod \
            $(DPO)fldvar.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)dashboard.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)sendrecv.mod \
            $(DPO)sendrecv3.mod \
            $(DPO)indices.mod \
            $(DPO)leqsol.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) mfix.f  -o $(DPO)mfix.$(OBJ_EXT) -module $(DPO)
$(DPO)mod_bc_i.$(OBJ_EXT) : mod_bc_i.f \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)error_manager.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) mod_bc_i.f  -o $(DPO)mod_bc_i.$(OBJ_EXT) -module $(DPO)
$(DPO)mod_bc_j.$(OBJ_EXT) : mod_bc_j.f \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)error_manager.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) mod_bc_j.f  -o $(DPO)mod_bc_j.$(OBJ_EXT) -module $(DPO)
$(DPO)mod_bc_k.$(OBJ_EXT) : mod_bc_k.f \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)error_manager.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) mod_bc_k.f  -o $(DPO)mod_bc_k.$(OBJ_EXT) -module $(DPO)
$(DPO)open_file.$(OBJ_EXT) : open_file.f \
            $(DPO)cdist.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) open_file.f  -o $(DPO)open_file.$(OBJ_EXT) -module $(DPO)
$(DPO)open_files.$(OBJ_EXT) : open_files.f \
            $(DPO)machine.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)cdist.mod \
            $(DPO)error_manager.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) open_files.f  -o $(DPO)open_files.$(OBJ_EXT) -module $(DPO)
$(DPO)out_array_c.$(OBJ_EXT) : out_array_c.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) out_array_c.f  -o $(DPO)out_array_c.$(OBJ_EXT) -module $(DPO)
$(DPO)out_array.$(OBJ_EXT) : out_array.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) out_array.f  -o $(DPO)out_array.$(OBJ_EXT) -module $(DPO)
$(DPO)out_array_kc.$(OBJ_EXT) : out_array_kc.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) out_array_kc.f  -o $(DPO)out_array_kc.$(OBJ_EXT) -module $(DPO)
$(DPO)out_array_k.$(OBJ_EXT) : out_array_k.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) out_array_k.f  -o $(DPO)out_array_k.$(OBJ_EXT) -module $(DPO)
$(DPO)out_bin_512.$(OBJ_EXT) : out_bin_512.f \
            $(DPO)machine.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) out_bin_512.f  -o $(DPO)out_bin_512.$(OBJ_EXT) -module $(DPO)
$(DPO)out_bin_512i.$(OBJ_EXT) : out_bin_512i.f \
            $(DPO)machine.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) out_bin_512i.f  -o $(DPO)out_bin_512i.$(OBJ_EXT) -module $(DPO)
$(DPO)out_bin_512r.$(OBJ_EXT) : out_bin_512r.f \
            $(DPO)machine.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) out_bin_512r.f  -o $(DPO)out_bin_512r.$(OBJ_EXT) -module $(DPO)
$(DPO)out_bin_r.$(OBJ_EXT) : out_bin_r.f \
            $(DPO)param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) out_bin_r.f  -o $(DPO)out_bin_r.$(OBJ_EXT) -module $(DPO)
$(DPO)parse_line.$(OBJ_EXT) : parse_line.f \
            $(DPO)compar.mod \
            $(DPO)des_rxns.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parse.mod \
            $(DPO)rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_line.f  -o $(DPO)parse_line.$(OBJ_EXT) -module $(DPO)
$(DPO)parse_resid_string.$(OBJ_EXT) : parse_resid_string.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)residual.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_resid_string.f  -o $(DPO)parse_resid_string.$(OBJ_EXT) -module $(DPO)
$(DPO)parse_rxn.$(OBJ_EXT) : parse_rxn.f \
            $(DPO)compar.mod \
            $(DPO)funits.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parse.mod \
            $(DPO)rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_rxn.f  -o $(DPO)parse_rxn.$(OBJ_EXT) -module $(DPO)
$(DPO)partial_elim.$(OBJ_EXT) : partial_elim.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)matrix.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)drag.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) partial_elim.f  -o $(DPO)partial_elim.$(OBJ_EXT) -module $(DPO)
$(DPO)physical_prop.$(OBJ_EXT) : physical_prop.f \
            $(DPO)compar.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)coeff.mod \
            $(DPO)fldvar.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)scalars.mod \
            $(DPO)cutcell.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) physical_prop.f  -o $(DPO)physical_prop.$(OBJ_EXT) -module $(DPO)
$(DPO)read_database.$(OBJ_EXT) : read_database.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod \
            $(DPO)rxns.mod \
            $(DPO)funits.mod \
            $(DPO)discretelement.mod \
            $(DPO)des_rxns.mod \
            $(DPO)error_manager.mod \
            mfix_directory_path.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) read_database.f  -o $(DPO)read_database.$(OBJ_EXT) -module $(DPO)
$(DPO)read_namelist.$(OBJ_EXT) : read_namelist.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)ic.mod \
            $(DPO)is.mod \
            $(DPO)bc.mod \
            $(DPO)ps.mod \
            $(DPO)fldvar.mod \
            $(DPO)constant.mod \
            $(DPO)indices.mod \
            $(DPO)toleranc.mod \
            $(DPO)funits.mod \
            $(DPO)scales.mod \
            $(DPO)ur_facs.mod \
            $(DPO)leqsol.mod \
            $(DPO)residual.mod \
            $(DPO)rxns.mod \
            $(DPO)scalars.mod \
            $(DPO)compar.mod \
            $(DPO)parallel.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)usr.mod \
            $(DPO)des_bc.mod \
            $(DPO)des_ic.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod \
            $(DPO)stiff_chem.mod \
            $(DPO)cdist.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            $(DPO)vtk.mod \
            $(DPO)polygon.mod \
            $(DPO)dashboard.mod \
            $(DPO)stl.mod \
            $(DPO)qmom_kinetic_equation.mod \
            usrnlst.inc                                                  \
            namelist.inc                                                 \
            des/desnamelist.inc                                          \
            cartesian_grid/cartesian_grid_namelist.inc                   \
            qmomk/qmomknamelist.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) read_namelist.f  -o $(DPO)read_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)read_res0.$(OBJ_EXT) : read_res0.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)ic.mod \
            $(DPO)bc.mod \
            $(DPO)is.mod \
            $(DPO)constant.mod \
            $(DPO)funits.mod \
            $(DPO)output.mod \
            $(DPO)scales.mod \
            $(DPO)ur_facs.mod \
            $(DPO)toleranc.mod \
            $(DPO)leqsol.mod \
            $(DPO)scalars.mod \
            $(DPO)rxns.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)fldvar.mod \
            $(DPO)stiff_chem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) read_res0.f  -o $(DPO)read_res0.$(OBJ_EXT) -module $(DPO)
$(DPO)read_res1.$(OBJ_EXT) : read_res1.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)rxns.mod \
            $(DPO)scalars.mod \
            $(DPO)funits.mod \
            $(DPO)energy.mod \
            $(DPO)compar.mod \
            $(DPO)cdist.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)mfix_netcdf.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) read_res1.f  -o $(DPO)read_res1.$(OBJ_EXT) -module $(DPO)
$(DPO)remove_comment.$(OBJ_EXT) : remove_comment.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) remove_comment.f  -o $(DPO)remove_comment.$(OBJ_EXT) -module $(DPO)
$(DPO)reset_new.$(OBJ_EXT) : reset_new.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)trace.mod \
            $(DPO)run.mod \
            $(DPO)scalars.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) reset_new.f  -o $(DPO)reset_new.$(OBJ_EXT) -module $(DPO)
$(DPO)rrates0.$(OBJ_EXT) : rrates0.f \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)rxns.mod \
            $(DPO)energy.mod \
            $(DPO)param.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)discretelement.mod \
            $(DPO)toleranc.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) rrates0.f  -o $(DPO)rrates0.$(OBJ_EXT) -module $(DPO)
$(DPO)rrates.$(OBJ_EXT) : rrates.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)rxns.mod \
            $(DPO)energy.mod \
            $(DPO)geometry.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) rrates.f  -o $(DPO)rrates.$(OBJ_EXT) -module $(DPO)
$(DPO)rrates_init.$(OBJ_EXT) : rrates_init.f \
            $(DPO)energy.mod \
            $(DPO)param1.mod \
            $(DPO)rxns.mod \
            $(DPO)stiff_chem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) rrates_init.f  -o $(DPO)rrates_init.$(OBJ_EXT) -module $(DPO)
$(DPO)scalar_prop.$(OBJ_EXT) : scalar_prop.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)scalars.mod \
            $(DPO)toleranc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) scalar_prop.f  -o $(DPO)scalar_prop.$(OBJ_EXT) -module $(DPO)
$(DPO)seek_comment.$(OBJ_EXT) : seek_comment.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) seek_comment.f  -o $(DPO)seek_comment.$(OBJ_EXT) -module $(DPO)
$(DPO)seek_end.$(OBJ_EXT) : seek_end.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) seek_end.f  -o $(DPO)seek_end.$(OBJ_EXT) -module $(DPO)
$(DPO)set_bc0.$(OBJ_EXT) : set_bc0.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)bc.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)funits.mod \
            $(DPO)scales.mod \
            $(DPO)scalars.mod \
            $(DPO)boundfunijk.mod \
            $(DPO)toleranc.mod \
            $(DPO)sendrecv.mod \
            $(DPO)mms.mod \
            sc_p_g1.inc                                                  \
            function.inc                                                 \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_bc0.f  -o $(DPO)set_bc0.$(OBJ_EXT) -module $(DPO)
$(DPO)set_bc1.$(OBJ_EXT) : set_bc1.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)bc.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_bc1.f  -o $(DPO)set_bc1.$(OBJ_EXT) -module $(DPO)
$(DPO)set_bc_flow.$(OBJ_EXT) : set_bc_flow.f \
            $(DPO)physprop.mod \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            $(DPO)bc.mod \
            $(DPO)param1.mod \
            $(DPO)param.mod \
            $(DPO)error_manager.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)scalars.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_bc_flow.f  -o $(DPO)set_bc_flow.$(OBJ_EXT) -module $(DPO)
$(DPO)set_constants.$(OBJ_EXT) : set_constants.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)energy.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)funits.mod \
            $(DPO)drag.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_constants.f  -o $(DPO)set_constants.$(OBJ_EXT) -module $(DPO)
$(DPO)set_constprop.$(OBJ_EXT) : set_constprop.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)visc_g.mod \
            $(DPO)energy.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)funits.mod \
            $(DPO)drag.mod \
            $(DPO)compar.mod \
            $(DPO)kintheory.mod \
            $(DPO)mms.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_constprop.f  -o $(DPO)set_constprop.$(OBJ_EXT) -module $(DPO)
$(DPO)set_flags.$(OBJ_EXT) : set_flags.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)bc.mod \
            $(DPO)is.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)sendrecv3.mod \
            $(DPO)boundfunijk.mod \
            $(DPO)mpi_utility.mod \
            function.inc                                                 \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) set_flags.f  -o $(DPO)set_flags.$(OBJ_EXT) -module $(DPO)
$(DPO)set_fluidbed_p.$(OBJ_EXT) : set_fluidbed_p.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)bc.mod \
            $(DPO)ic.mod \
            $(DPO)fldvar.mod \
            $(DPO)constant.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)scales.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            sc_p_g1.inc                                                  \
            b_force1.inc                                                 \
            function.inc                                                 \
            b_force2.inc                                                 \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_fluidbed_p.f  -o $(DPO)set_fluidbed_p.$(OBJ_EXT) -module $(DPO)
$(DPO)set_geometry1.$(OBJ_EXT) : set_geometry1.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_geometry1.f  -o $(DPO)set_geometry1.$(OBJ_EXT) -module $(DPO)
$(DPO)set_geometry.$(OBJ_EXT) : set_geometry.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)bc.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_geometry.f  -o $(DPO)set_geometry.$(OBJ_EXT) -module $(DPO)
$(DPO)set_icbc_flags.$(OBJ_EXT) : set_icbc_flags.f \
            $(DPO)run.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)error_manager.mod \
            $(DPO)ic.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)bc.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_icbc_flags.f  -o $(DPO)set_icbc_flags.$(OBJ_EXT) -module $(DPO)
$(DPO)set_ic.$(OBJ_EXT) : set_ic.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)ic.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)indices.mod \
            $(DPO)scales.mod \
            $(DPO)energy.mod \
            $(DPO)scalars.mod \
            $(DPO)compar.mod \
            $(DPO)run.mod \
            $(DPO)sendrecv.mod \
            sc_p_g1.inc                                                  \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_ic.f  -o $(DPO)set_ic.$(OBJ_EXT) -module $(DPO)
$(DPO)set_increments3.$(OBJ_EXT) : set_increments3.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)indices.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)funits.mod \
            function.inc                                                 \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) set_increments3.f  -o $(DPO)set_increments3.$(OBJ_EXT) -module $(DPO)
$(DPO)set_increments.$(OBJ_EXT) : set_increments.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)indices.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)funits.mod \
            $(DPO)scalars.mod \
            $(DPO)run.mod \
            $(DPO)visc_g.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)cutcell.mod \
            $(DPO)stl.mod \
            $(DPO)sendrecv.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)parallel.mod \
            $(DPO)bc.mod \
            $(DPO)discretelement.mod \
            $(DPO)cdist.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_increments.f  -o $(DPO)set_increments.$(OBJ_EXT) -module $(DPO)
$(DPO)set_index1a3.$(OBJ_EXT) : set_index1a3.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)boundfunijk3.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_index1a3.f  -o $(DPO)set_index1a3.$(OBJ_EXT) -module $(DPO)
$(DPO)set_index1a.$(OBJ_EXT) : set_index1a.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)boundfunijk.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_index1a.f  -o $(DPO)set_index1a.$(OBJ_EXT) -module $(DPO)
$(DPO)set_index1.$(OBJ_EXT) : set_index1.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)constant.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_index1.f  -o $(DPO)set_index1.$(OBJ_EXT) -module $(DPO)
$(DPO)set_l_scale.$(OBJ_EXT) : set_l_scale.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)visc_g.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_l_scale.f  -o $(DPO)set_l_scale.$(OBJ_EXT) -module $(DPO)
$(DPO)set_max2.$(OBJ_EXT) : set_max2.f \
            $(DPO)compar.mod \
            $(DPO)geometry.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_max2.f  -o $(DPO)set_max2.$(OBJ_EXT) -module $(DPO)
$(DPO)set_mw_mix_g.$(OBJ_EXT) : set_mw_mix_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)constant.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_mw_mix_g.f  -o $(DPO)set_mw_mix_g.$(OBJ_EXT) -module $(DPO)
$(DPO)set_outflow.$(OBJ_EXT) : set_outflow.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)bc.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)scalars.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)mflux.mod \
            $(DPO)discretelement.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) set_outflow.f  -o $(DPO)set_outflow.$(OBJ_EXT) -module $(DPO)
$(DPO)set_ps.$(OBJ_EXT) : set_ps.f \
            $(DPO)param.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)ps.mod \
            $(DPO)compar.mod \
            $(DPO)geometry.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)constant.mod \
            $(DPO)bc.mod \
            $(DPO)cutcell.mod \
            $(DPO)fldvar.mod \
            $(DPO)ic.mod \
            $(DPO)indices.mod \
            $(DPO)mflux.mod \
            $(DPO)parallel.mod \
            $(DPO)param1.mod \
            $(DPO)sendrecv.mod \
            $(DPO)toleranc.mod \
            $(DPO)usr.mod \
            $(DPO)rxns.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_ps.f  -o $(DPO)set_ps.$(OBJ_EXT) -module $(DPO)
$(DPO)set_ro_g.$(OBJ_EXT) : set_ro_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)constant.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_ro_g.f  -o $(DPO)set_ro_g.$(OBJ_EXT) -module $(DPO)
$(DPO)set_ro_s.$(OBJ_EXT) : set_ro_s.f \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_ro_s.f  -o $(DPO)set_ro_s.$(OBJ_EXT) -module $(DPO)
$(DPO)set_wall_bc.$(OBJ_EXT) : set_wall_bc.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)bc.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_wall_bc.f  -o $(DPO)set_wall_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)shift_dxyz.$(OBJ_EXT) : shift_dxyz.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) shift_dxyz.f  -o $(DPO)shift_dxyz.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_continuity.$(OBJ_EXT) : solve_continuity.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)residual.mod \
            $(DPO)cont.mod \
            $(DPO)leqsol.mod \
            $(DPO)ambm.mod \
            $(DPO)ur_facs.mod \
            $(DPO)run.mod \
            $(DPO)ps.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_continuity.f  -o $(DPO)solve_continuity.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_energy_eq.$(OBJ_EXT) : solve_energy_eq.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)drag.mod \
            $(DPO)residual.mod \
            $(DPO)ur_facs.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)leqsol.mod \
            $(DPO)bc.mod \
            $(DPO)energy.mod \
            $(DPO)rxns.mod \
            $(DPO)ambm.mod \
            $(DPO)tmp_array.mod \
            $(DPO)tmp_array1.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)des_thermo.mod \
            $(DPO)mflux.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)ps.mod \
            $(DPO)mms.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_energy_eq.f  -o $(DPO)solve_energy_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_epp.$(OBJ_EXT) : solve_epp.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)pscor.mod \
            $(DPO)residual.mod \
            $(DPO)leqsol.mod \
            $(DPO)physprop.mod \
            $(DPO)ambm.mod \
            $(DPO)tmp_array1.mod \
            $(DPO)ps.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_epp.f  -o $(DPO)solve_epp.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_granular_energy.$(OBJ_EXT) : solve_granular_energy.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)visc_s.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)constant.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)drag.mod \
            $(DPO)residual.mod \
            $(DPO)ur_facs.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)leqsol.mod \
            $(DPO)bc.mod \
            $(DPO)energy.mod \
            $(DPO)rxns.mod \
            $(DPO)ambm.mod \
            $(DPO)tmp_array.mod \
            $(DPO)compar.mod \
            $(DPO)mflux.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)mms.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_granular_energy.f  -o $(DPO)solve_granular_energy.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_k_epsilon_eq.$(OBJ_EXT) : solve_k_epsilon_eq.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)constant.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)drag.mod \
            $(DPO)residual.mod \
            $(DPO)ur_facs.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)leqsol.mod \
            $(DPO)bc.mod \
            $(DPO)energy.mod \
            $(DPO)rxns.mod \
            $(DPO)turb.mod \
            $(DPO)usr.mod \
            $(DPO)ambm.mod \
            $(DPO)tmp_array.mod \
            $(DPO)compar.mod \
            $(DPO)mflux.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_k_epsilon_eq.f  -o $(DPO)solve_k_epsilon_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_lin_eq.$(OBJ_EXT) : solve_lin_eq.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)residual.mod \
            $(DPO)toleranc.mod \
            $(DPO)leqsol.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_lin_eq.f  -o $(DPO)solve_lin_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_pp_g.$(OBJ_EXT) : solve_pp_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)pgcor.mod \
            $(DPO)residual.mod \
            $(DPO)leqsol.mod \
            $(DPO)run.mod \
            $(DPO)ambm.mod \
            $(DPO)tmp_array1.mod \
            $(DPO)ps.mod \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_pp_g.f  -o $(DPO)solve_pp_g.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_scalar_eq.$(OBJ_EXT) : solve_scalar_eq.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)drag.mod \
            $(DPO)residual.mod \
            $(DPO)ur_facs.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)leqsol.mod \
            $(DPO)bc.mod \
            $(DPO)energy.mod \
            $(DPO)rxns.mod \
            $(DPO)scalars.mod \
            $(DPO)ambm.mod \
            $(DPO)tmp_array.mod \
            $(DPO)compar.mod \
            $(DPO)mflux.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_scalar_eq.f  -o $(DPO)solve_scalar_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_species_eq.$(OBJ_EXT) : solve_species_eq.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)drag.mod \
            $(DPO)residual.mod \
            $(DPO)ur_facs.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)leqsol.mod \
            $(DPO)bc.mod \
            $(DPO)energy.mod \
            $(DPO)rxns.mod \
            $(DPO)ambm.mod \
            $(DPO)matrix.mod \
            $(DPO)chischeme.mod \
            $(DPO)tmp_array.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)mflux.mod \
            $(DPO)ps.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_species_eq.f  -o $(DPO)solve_species_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_vel_star.$(OBJ_EXT) : solve_vel_star.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)drag.mod \
            $(DPO)residual.mod \
            $(DPO)ur_facs.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)leqsol.mod \
            $(DPO)ambm.mod \
            $(DPO)tmp_array1.mod \
            $(DPO)tmp_array.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)ps.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_vel_star.f  -o $(DPO)solve_vel_star.$(OBJ_EXT) -module $(DPO)
$(DPO)source_granular_energy.$(OBJ_EXT) : source_granular_energy.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)drag.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)visc_s.mod \
            $(DPO)trace.mod \
            $(DPO)turb.mod \
            $(DPO)indices.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)compar.mod \
            $(DPO)kintheory.mod \
            $(DPO)mms.mod \
            $(DPO)residual.mod \
            s_pr1.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            s_pr2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) source_granular_energy.f  -o $(DPO)source_granular_energy.$(OBJ_EXT) -module $(DPO)
$(DPO)source_phi.$(OBJ_EXT) : source_phi.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_s.mod \
            $(DPO)compar.mod \
            $(DPO)ps.mod \
            $(DPO)bc.mod \
            $(DPO)usr.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) source_phi.f  -o $(DPO)source_phi.$(OBJ_EXT) -module $(DPO)
$(DPO)source_pp_g.$(OBJ_EXT) : source_pp_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)pgcor.mod \
            $(DPO)bc.mod \
            $(DPO)vshear.mod \
            $(DPO)xsi_array.mod \
            $(DPO)compar.mod \
            $(DPO)ur_facs.mod \
            $(DPO)constant.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_pp_g.f  -o $(DPO)source_pp_g.$(OBJ_EXT) -module $(DPO)
$(DPO)source_rop_g.$(OBJ_EXT) : source_rop_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)fldvar.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)pgcor.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_rop_g.f  -o $(DPO)source_rop_g.$(OBJ_EXT) -module $(DPO)
$(DPO)source_rop_s.$(OBJ_EXT) : source_rop_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)fldvar.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)ps.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_rop_s.f  -o $(DPO)source_rop_s.$(OBJ_EXT) -module $(DPO)
$(DPO)source_u_g.$(OBJ_EXT) : source_u_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_g.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)mms.mod \
            $(DPO)output.mod \
            $(DPO)turb.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)ps.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_u_g.f  -o $(DPO)source_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)source_u_s.$(OBJ_EXT) : source_u_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_s.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)kintheory.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)mms.mod \
            $(DPO)output.mod \
            $(DPO)ps.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_u_s.f  -o $(DPO)source_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)source_v_g.$(OBJ_EXT) : source_v_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_g.mod \
            $(DPO)bc.mod \
            $(DPO)vshear.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)mms.mod \
            $(DPO)output.mod \
            $(DPO)ps.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_v_g.f  -o $(DPO)source_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)source_v_s.$(OBJ_EXT) : source_v_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_s.mod \
            $(DPO)bc.mod \
            $(DPO)vshear.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)kintheory.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)mms.mod \
            $(DPO)output.mod \
            $(DPO)ps.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_v_s.f  -o $(DPO)source_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)source_w_g.$(OBJ_EXT) : source_w_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_g.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)mms.mod \
            $(DPO)output.mod \
            $(DPO)ps.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_w_g.f  -o $(DPO)source_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)source_w_s.$(OBJ_EXT) : source_w_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_s.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)kintheory.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)mms.mod \
            $(DPO)output.mod \
            $(DPO)ps.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_w_s.f  -o $(DPO)source_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_u_g.$(OBJ_EXT) : tau_u_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)bc.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_u_g.f  -o $(DPO)tau_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_u_s.$(OBJ_EXT) : tau_u_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)vshear.mod \
            $(DPO)sendrecv.mod \
            $(DPO)compar.mod \
            $(DPO)bc.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_u_s.f  -o $(DPO)tau_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_v_g.$(OBJ_EXT) : tau_v_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)sendrecv.mod \
            $(DPO)compar.mod \
            $(DPO)bc.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_v_g.f  -o $(DPO)tau_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_v_s.$(OBJ_EXT) : tau_v_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)sendrecv.mod \
            $(DPO)compar.mod \
            $(DPO)bc.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_v_s.f  -o $(DPO)tau_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_w_g.$(OBJ_EXT) : tau_w_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)sendrecv.mod \
            $(DPO)compar.mod \
            $(DPO)bc.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_w_g.f  -o $(DPO)tau_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_w_s.$(OBJ_EXT) : tau_w_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)sendrecv.mod \
            $(DPO)compar.mod \
            $(DPO)bc.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_w_s.f  -o $(DPO)tau_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)test_lin_eq.$(OBJ_EXT) : test_lin_eq.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) test_lin_eq.f  -o $(DPO)test_lin_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)time_march.$(OBJ_EXT) : time_march.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)cont.mod \
            $(DPO)tau_g.mod \
            $(DPO)tau_s.mod \
            $(DPO)visc_g.mod \
            $(DPO)visc_s.mod \
            $(DPO)funits.mod \
            $(DPO)vshear.mod \
            $(DPO)scalars.mod \
            $(DPO)toleranc.mod \
            $(DPO)drag.mod \
            $(DPO)rxns.mod \
            $(DPO)compar.mod \
            $(DPO)time_cpu.mod \
            $(DPO)discretelement.mod \
            $(DPO)leqsol.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)cdist.mod \
            $(DPO)mfix_netcdf.mod \
            $(DPO)cutcell.mod \
            $(DPO)vtk.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)dashboard.mod \
            $(DPO)indices.mod \
            $(DPO)bc.mod \
            $(DPO)coeff.mod \
            $(DPO)stiff_chem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) time_march.f  -o $(DPO)time_march.$(OBJ_EXT) -module $(DPO)
$(DPO)transfer.$(OBJ_EXT) : transfer.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) transfer.f  -o $(DPO)transfer.$(OBJ_EXT) -module $(DPO)
$(DPO)transport_prop.$(OBJ_EXT) : transport_prop.f \
            $(DPO)physprop.mod \
            $(DPO)coeff.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) transport_prop.f  -o $(DPO)transport_prop.$(OBJ_EXT) -module $(DPO)
$(DPO)undef_2_0.$(OBJ_EXT) : undef_2_0.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) undef_2_0.f  -o $(DPO)undef_2_0.$(OBJ_EXT) -module $(DPO)
$(DPO)under_relax.$(OBJ_EXT) : under_relax.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) under_relax.f  -o $(DPO)under_relax.$(OBJ_EXT) -module $(DPO)
$(DPO)update_old.$(OBJ_EXT) : update_old.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)trace.mod \
            $(DPO)visc_s.mod \
            $(DPO)scalars.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) update_old.f  -o $(DPO)update_old.$(OBJ_EXT) -module $(DPO)
$(DPO)usr0.$(OBJ_EXT) : usr0.f \
            $(DPO)usr.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr0.f  -o $(DPO)usr0.$(OBJ_EXT) -module $(DPO)
$(DPO)usr1.$(OBJ_EXT) : usr1.f \
            $(DPO)usr.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr1.f  -o $(DPO)usr1.$(OBJ_EXT) -module $(DPO)
$(DPO)usr2.$(OBJ_EXT) : usr2.f \
            $(DPO)usr.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr2.f  -o $(DPO)usr2.$(OBJ_EXT) -module $(DPO)
$(DPO)usr3.$(OBJ_EXT) : usr3.f \
            $(DPO)usr.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr3.f  -o $(DPO)usr3.$(OBJ_EXT) -module $(DPO)
$(DPO)usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f \
            $(DPO)usr.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_init_namelist.f  -o $(DPO)usr_init_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)usr_rates.$(OBJ_EXT) : usr_rates.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)rxns.mod \
            $(DPO)energy.mod \
            $(DPO)geometry.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)toleranc.mod \
            $(DPO)usr.mod \
            species.inc                                                  \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_rates.f  -o $(DPO)usr_rates.$(OBJ_EXT) -module $(DPO)
$(DPO)usr_write_out0.$(OBJ_EXT) : usr_write_out0.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_write_out0.f  -o $(DPO)usr_write_out0.$(OBJ_EXT) -module $(DPO)
$(DPO)usr_write_out1.$(OBJ_EXT) : usr_write_out1.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_write_out1.f  -o $(DPO)usr_write_out1.$(OBJ_EXT) -module $(DPO)
$(DPO)utilities.$(OBJ_EXT) : utilities.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)bc.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)toleranc.mod \
            $(DPO)mpi_utility.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) utilities.f  -o $(DPO)utilities.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_u_g.$(OBJ_EXT) : vavg_u_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)mflux.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_u_g.f  -o $(DPO)vavg_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_u_s.$(OBJ_EXT) : vavg_u_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_u_s.f  -o $(DPO)vavg_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_v_g.$(OBJ_EXT) : vavg_v_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)mflux.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_v_g.f  -o $(DPO)vavg_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_v_s.$(OBJ_EXT) : vavg_v_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_v_s.f  -o $(DPO)vavg_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_w_g.$(OBJ_EXT) : vavg_w_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)mflux.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_w_g.f  -o $(DPO)vavg_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_w_s.$(OBJ_EXT) : vavg_w_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_w_s.f  -o $(DPO)vavg_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)vf_gs_x.$(OBJ_EXT) : vf_gs_x.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)compar.mod \
            $(DPO)drag.mod \
            $(DPO)discretelement.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vf_gs_x.f  -o $(DPO)vf_gs_x.$(OBJ_EXT) -module $(DPO)
$(DPO)vf_gs_y.$(OBJ_EXT) : vf_gs_y.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)compar.mod \
            $(DPO)drag.mod \
            $(DPO)discretelement.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vf_gs_y.f  -o $(DPO)vf_gs_y.$(OBJ_EXT) -module $(DPO)
$(DPO)vf_gs_z.$(OBJ_EXT) : vf_gs_z.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)compar.mod \
            $(DPO)drag.mod \
            $(DPO)discretelement.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vf_gs_z.f  -o $(DPO)vf_gs_z.$(OBJ_EXT) -module $(DPO)
$(DPO)vtc_scalar.$(OBJ_EXT) : vtc_scalar.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)compar.mod \
            $(DPO)kintheory.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vtc_scalar.f  -o $(DPO)vtc_scalar.$(OBJ_EXT) -module $(DPO)
$(DPO)write_ab_m.$(OBJ_EXT) : write_ab_m.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) write_ab_m.f  -o $(DPO)write_ab_m.$(OBJ_EXT) -module $(DPO)
$(DPO)write_ab_m_var.$(OBJ_EXT) : write_ab_m_var.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) write_ab_m_var.f  -o $(DPO)write_ab_m_var.$(OBJ_EXT) -module $(DPO)
$(DPO)write_error.$(OBJ_EXT) : write_error.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_error.f  -o $(DPO)write_error.$(OBJ_EXT) -module $(DPO)
$(DPO)write_header.$(OBJ_EXT) : write_header.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_header.f  -o $(DPO)write_header.$(OBJ_EXT) -module $(DPO)
$(DPO)write_out0.$(OBJ_EXT) : write_out0.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)ic.mod \
            $(DPO)bc.mod \
            $(DPO)is.mod \
            $(DPO)fldvar.mod \
            $(DPO)constant.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)toleranc.mod \
            $(DPO)scales.mod \
            $(DPO)scalars.mod \
            $(DPO)ur_facs.mod \
            $(DPO)leqsol.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)rxns.mod \
            $(DPO)mfix_pic.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) write_out0.f  -o $(DPO)write_out0.$(OBJ_EXT) -module $(DPO)
$(DPO)write_out1.$(OBJ_EXT) : write_out1.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)scalars.mod \
            $(DPO)funits.mod \
            $(DPO)rxns.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_out1.f  -o $(DPO)write_out1.$(OBJ_EXT) -module $(DPO)
$(DPO)write_out3.$(OBJ_EXT) : write_out3.f \
            $(DPO)funits.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_out3.f  -o $(DPO)write_out3.$(OBJ_EXT) -module $(DPO)
$(DPO)write_res0.$(OBJ_EXT) : write_res0.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)ic.mod \
            $(DPO)is.mod \
            $(DPO)bc.mod \
            $(DPO)constant.mod \
            $(DPO)funits.mod \
            $(DPO)output.mod \
            $(DPO)scales.mod \
            $(DPO)scalars.mod \
            $(DPO)rxns.mod \
            $(DPO)ur_facs.mod \
            $(DPO)leqsol.mod \
            $(DPO)toleranc.mod \
            $(DPO)cdist.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)stiff_chem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_res0.f  -o $(DPO)write_res0.$(OBJ_EXT) -module $(DPO)
$(DPO)write_res1.$(OBJ_EXT) : write_res1.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)scalars.mod \
            $(DPO)rxns.mod \
            $(DPO)funits.mod \
            $(DPO)output.mod \
            $(DPO)energy.mod \
            $(DPO)cdist.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)cutcell.mod \
            $(DPO)mfix_netcdf.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_res1.f  -o $(DPO)write_res1.$(OBJ_EXT) -module $(DPO)
$(DPO)write_spx0.$(OBJ_EXT) : write_spx0.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)funits.mod \
            $(DPO)cdist.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_spx0.f  -o $(DPO)write_spx0.$(OBJ_EXT) -module $(DPO)
$(DPO)write_spx1.$(OBJ_EXT) : write_spx1.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)funits.mod \
            $(DPO)scalars.mod \
            $(DPO)output.mod \
            $(DPO)rxns.mod \
            $(DPO)cdist.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)cutcell.mod \
            $(DPO)sendrecv.mod \
            $(DPO)mfix_netcdf.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_spx1.f  -o $(DPO)write_spx1.$(OBJ_EXT) -module $(DPO)
$(DPO)write_table.$(OBJ_EXT) : write_table.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_table.f  -o $(DPO)write_table.$(OBJ_EXT) -module $(DPO)
$(DPO)write_usr0.$(OBJ_EXT) : write_usr0.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_usr0.f  -o $(DPO)write_usr0.$(OBJ_EXT) -module $(DPO)
$(DPO)write_usr1.$(OBJ_EXT) : write_usr1.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_usr1.f  -o $(DPO)write_usr1.$(OBJ_EXT) -module $(DPO)
$(DPO)xerbla.$(OBJ_EXT) : xerbla.f \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) xerbla.f  -o $(DPO)xerbla.$(OBJ_EXT) -module $(DPO)
$(DPO)zero_array.$(OBJ_EXT) : zero_array.f \
            $(DPO)param.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) zero_array.f  -o $(DPO)zero_array.$(OBJ_EXT) -module $(DPO)
$(DPO)zero_norm_vel.$(OBJ_EXT) : zero_norm_vel.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) zero_norm_vel.f  -o $(DPO)zero_norm_vel.$(OBJ_EXT) -module $(DPO)
$(DPO)allocate_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/allocate_cut_cell_arrays.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)indices.mod \
            $(DPO)cutcell.mod \
            $(DPO)stl.mod \
            $(DPO)discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/allocate_cut_cell_arrays.f  -o $(DPO)allocate_cut_cell_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)allocate_dummy_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/allocate_dummy_cut_cell_arrays.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)indices.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/allocate_dummy_cut_cell_arrays.f  -o $(DPO)allocate_dummy_cut_cell_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_vort_out.$(OBJ_EXT) : ./cartesian_grid/calc_vort_out.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)fldvar.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/calc_vort_out.f  -o $(DPO)calc_vort_out.$(OBJ_EXT) -module $(DPO)
$(DPO)cartesian_grid_init_namelist.$(OBJ_EXT) : ./cartesian_grid/cartesian_grid_init_namelist.f \
            $(DPO)param1.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            $(DPO)polygon.mod \
            $(DPO)vtk.mod \
            $(DPO)progress_bar.mod \
            $(DPO)dashboard.mod \
            $(DPO)stl.mod \
            cartesian_grid/cartesian_grid_namelist.inc                  
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cartesian_grid_init_namelist.f  -o $(DPO)cartesian_grid_init_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_set_bc0.$(OBJ_EXT) : ./cartesian_grid/CG_set_bc0.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)physprop.mod \
            $(DPO)bc.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)funits.mod \
            $(DPO)scales.mod \
            $(DPO)scalars.mod \
            $(DPO)boundfunijk.mod \
            $(DPO)toleranc.mod \
            $(DPO)sendrecv.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            sc_p_g1.inc                                                  \
            function.inc                                                 \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_set_bc0.f  -o $(DPO)CG_set_bc0.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_set_outflow.$(OBJ_EXT) : ./cartesian_grid/CG_set_outflow.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)bc.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)scalars.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)mflux.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_set_outflow.f  -o $(DPO)CG_set_outflow.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_u_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_u_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_g.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_u_g.f  -o $(DPO)CG_source_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_u_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_u_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_s.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)kintheory.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_u_s.f  -o $(DPO)CG_source_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_v_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_v_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_g.mod \
            $(DPO)bc.mod \
            $(DPO)vshear.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_v_g.f  -o $(DPO)CG_source_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_v_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_v_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_s.mod \
            $(DPO)bc.mod \
            $(DPO)vshear.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)kintheory.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_v_s.f  -o $(DPO)CG_source_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_w_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_w_g.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_g.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_w_g.f  -o $(DPO)CG_source_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_w_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_w_s.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_s.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)kintheory.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)drag.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_w_s.f  -o $(DPO)CG_source_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_cartesian.$(OBJ_EXT) : ./cartesian_grid/check_data_cartesian.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)scalars.mod \
            $(DPO)funits.mod \
            $(DPO)leqsol.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)bc.mod \
            $(DPO)discretelement.mod \
            $(DPO)cutcell.mod \
            $(DPO)quadric.mod \
            $(DPO)vtk.mod \
            $(DPO)polygon.mod \
            $(DPO)dashboard.mod \
            $(DPO)stl.mod \
            $(DPO)fldvar.mod \
            $(DPO)scales.mod \
            $(DPO)parallel.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)sendrecv.mod \
            $(DPO)ps.mod \
            $(DPO)gridmap.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/check_data_cartesian.f  -o $(DPO)check_data_cartesian.$(OBJ_EXT) -module $(DPO)
$(DPO)cut_cell_preprocessing.$(OBJ_EXT) : ./cartesian_grid/cut_cell_preprocessing.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            $(DPO)vtk.mod \
            $(DPO)cdist.mod \
            $(DPO)discretelement.mod \
            $(DPO)des_stl_functions.mod \
            $(DPO)fldvar.mod \
            $(DPO)polygon.mod \
            $(DPO)stl.mod \
            $(DPO)stl.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cut_cell_preprocessing.f  -o $(DPO)cut_cell_preprocessing.$(OBJ_EXT) -module $(DPO)
$(DPO)deallocate_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/deallocate_cut_cell_arrays.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)indices.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/deallocate_cut_cell_arrays.f  -o $(DPO)deallocate_cut_cell_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)define_quadrics.$(OBJ_EXT) : ./cartesian_grid/define_quadrics.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            $(DPO)vtk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/define_quadrics.f  -o $(DPO)define_quadrics.$(OBJ_EXT) -module $(DPO)
$(DPO)dmp_cartesian.$(OBJ_EXT) : ./cartesian_grid/dmp_cartesian.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/dmp_cartesian.f  -o $(DPO)dmp_cartesian.$(OBJ_EXT) -module $(DPO)
$(DPO)eval_usr_fct.$(OBJ_EXT) : ./cartesian_grid/eval_usr_fct.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)fldvar.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/eval_usr_fct.f  -o $(DPO)eval_usr_fct.$(OBJ_EXT) -module $(DPO)
$(DPO)get_alpha.$(OBJ_EXT) : ./cartesian_grid/get_alpha.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)bc.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_alpha.f  -o $(DPO)get_alpha.$(OBJ_EXT) -module $(DPO)
$(DPO)get_connectivity.$(OBJ_EXT) : ./cartesian_grid/get_connectivity.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            $(DPO)polygon.mod \
            $(DPO)stl.mod \
            $(DPO)fldvar.mod \
            $(DPO)vtk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_connectivity.f  -o $(DPO)get_connectivity.$(OBJ_EXT) -module $(DPO)
$(DPO)get_cut_cell_flags.$(OBJ_EXT) : ./cartesian_grid/get_cut_cell_flags.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            $(DPO)vtk.mod \
            $(DPO)polygon.mod \
            $(DPO)stl.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)scalars.mod \
            $(DPO)funits.mod \
            $(DPO)rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_cut_cell_flags.f  -o $(DPO)get_cut_cell_flags.$(OBJ_EXT) -module $(DPO)
$(DPO)get_cut_cell_volume_area.$(OBJ_EXT) : ./cartesian_grid/get_cut_cell_volume_area.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            $(DPO)polygon.mod \
            $(DPO)stl.mod \
            $(DPO)bc.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_cut_cell_volume_area.f  -o $(DPO)get_cut_cell_volume_area.$(OBJ_EXT) -module $(DPO)
$(DPO)get_delh.$(OBJ_EXT) : ./cartesian_grid/get_delh.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_delh.f  -o $(DPO)get_delh.$(OBJ_EXT) -module $(DPO)
$(DPO)get_master.$(OBJ_EXT) : ./cartesian_grid/get_master.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)bc.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_master.f  -o $(DPO)get_master.$(OBJ_EXT) -module $(DPO)
$(DPO)get_poly_data.$(OBJ_EXT) : ./cartesian_grid/get_poly_data.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)scalars.mod \
            $(DPO)funits.mod \
            $(DPO)rxns.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)progress_bar.mod \
            $(DPO)polygon.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_poly_data.f  -o $(DPO)get_poly_data.$(OBJ_EXT) -module $(DPO)
$(DPO)get_stl_data.$(OBJ_EXT) : ./cartesian_grid/get_stl_data.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)scalars.mod \
            $(DPO)funits.mod \
            $(DPO)rxns.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)progress_bar.mod \
            $(DPO)stl.mod \
            $(DPO)vtk.mod \
            $(DPO)quadric.mod \
            $(DPO)constant.mod \
            $(DPO)bc.mod \
            $(DPO)cutcell.mod \
            $(DPO)parallel.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)sendrecv.mod \
            $(DPO)stl.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_stl_data.f  -o $(DPO)get_stl_data.$(OBJ_EXT) -module $(DPO)
$(DPO)set_Odxyz.$(OBJ_EXT) : ./cartesian_grid/set_Odxyz.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            $(DPO)vtk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/set_Odxyz.f  -o $(DPO)set_Odxyz.$(OBJ_EXT) -module $(DPO)
$(DPO)update_dashboard.$(OBJ_EXT) : ./cartesian_grid/update_dashboard.f \
            $(DPO)compar.mod \
            $(DPO)parallel.mod \
            $(DPO)sendrecv.mod \
            $(DPO)run.mod \
            $(DPO)leqsol.mod \
            $(DPO)time_cpu.mod \
            $(DPO)residual.mod \
            $(DPO)dashboard.mod \
            $(DPO)vtk.mod \
            $(DPO)constant.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/update_dashboard.f  -o $(DPO)update_dashboard.$(OBJ_EXT) -module $(DPO)
$(DPO)vtk_out.$(OBJ_EXT) : ./cartesian_grid/vtk_out.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)quadric.mod \
            $(DPO)cutcell.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)physprop.mod \
            $(DPO)pgcor.mod \
            $(DPO)vtk.mod \
            $(DPO)rxns.mod \
            $(DPO)output.mod \
            $(DPO)scalars.mod \
            $(DPO)stl.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)parallel_mpi.mod \
            $(DPO)pscor.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cdist.mod \
            $(DPO)polygon.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/vtk_out.f  -o $(DPO)vtk_out.$(OBJ_EXT) -module $(DPO)
$(DPO)write_progress_bar.$(OBJ_EXT) : ./cartesian_grid/write_progress_bar.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)scalars.mod \
            $(DPO)funits.mod \
            $(DPO)rxns.mod \
            $(DPO)compar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)progress_bar.mod \
            $(DPO)parallel.mod \
            $(DPO)sendrecv.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/write_progress_bar.f  -o $(DPO)write_progress_bar.$(OBJ_EXT) -module $(DPO)
$(DPO)check_axis.$(OBJ_EXT) : ./check_data/check_axis.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_axis.f  -o $(DPO)check_axis.$(OBJ_EXT) -module $(DPO)
$(DPO)check_bc_geometry.$(OBJ_EXT) : ./check_data/check_bc_geometry.f \
            $(DPO)bc.mod \
            $(DPO)geometry.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)error_manager.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_bc_geometry.f  -o $(DPO)check_bc_geometry.$(OBJ_EXT) -module $(DPO)
$(DPO)check_bc_inflow.$(OBJ_EXT) : ./check_data/check_bc_inflow.f \
            $(DPO)run.mod \
            $(DPO)scalars.mod \
            $(DPO)physprop.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)bc.mod \
            $(DPO)error_manager.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_bc_inflow.f  -o $(DPO)check_bc_inflow.$(OBJ_EXT) -module $(DPO)
$(DPO)check_bc_outflow.$(OBJ_EXT) : ./check_data/check_bc_outflow.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)discretelement.mod \
            $(DPO)bc.mod \
            $(DPO)error_manager.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)scalars.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_bc_outflow.f  -o $(DPO)check_bc_outflow.$(OBJ_EXT) -module $(DPO)
$(DPO)check_bc_walls.$(OBJ_EXT) : ./check_data/check_bc_walls.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)bc.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)scalars.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cutcell.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_bc_walls.f  -o $(DPO)check_bc_walls.$(OBJ_EXT) -module $(DPO)
$(DPO)check_boundary_conditions.$(OBJ_EXT) : ./check_data/check_boundary_conditions.f \
            $(DPO)physprop.mod \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            $(DPO)bc.mod \
            $(DPO)param1.mod \
            $(DPO)param.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_boundary_conditions.f  -o $(DPO)check_boundary_conditions.$(OBJ_EXT) -module $(DPO)
$(DPO)check_dmp_prereqs.$(OBJ_EXT) : ./check_data/check_dmp_prereqs.f \
            $(DPO)compar.mod \
            $(DPO)geometry.mod \
            $(DPO)param1.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_dmp_prereqs.f  -o $(DPO)check_dmp_prereqs.$(OBJ_EXT) -module $(DPO)
$(DPO)check_gas_phase.$(OBJ_EXT) : ./check_data/check_gas_phase.f \
            $(DPO)run.mod \
            $(DPO)rxns.mod \
            $(DPO)physprop.mod \
            $(DPO)param1.mod \
            $(DPO)error_manager.mod \
            $(DPO)param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_gas_phase.f  -o $(DPO)check_gas_phase.$(OBJ_EXT) -module $(DPO)
$(DPO)check_geometry.$(OBJ_EXT) : ./check_data/check_geometry.f \
            $(DPO)geometry.mod \
            $(DPO)bc.mod \
            $(DPO)cutcell.mod \
            $(DPO)param1.mod \
            $(DPO)param.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_geometry.f  -o $(DPO)check_geometry.$(OBJ_EXT) -module $(DPO)
$(DPO)check_geometry_prereqs.$(OBJ_EXT) : ./check_data/check_geometry_prereqs.f \
            $(DPO)geometry.mod \
            $(DPO)bc.mod \
            $(DPO)cutcell.mod \
            $(DPO)param1.mod \
            $(DPO)param.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_geometry_prereqs.f  -o $(DPO)check_geometry_prereqs.$(OBJ_EXT) -module $(DPO)
$(DPO)check_ic_common_discrete.$(OBJ_EXT) : ./check_data/check_ic_common_discrete.f \
            $(DPO)discretelement.mod \
            $(DPO)ic.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)param.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_ic_common_discrete.f  -o $(DPO)check_ic_common_discrete.$(OBJ_EXT) -module $(DPO)
$(DPO)check_initial_conditions_dem.$(OBJ_EXT) : ./check_data/check_initial_conditions_dem.f \
            $(DPO)discretelement.mod \
            $(DPO)ic.mod \
            $(DPO)param1.mod \
            $(DPO)param.mod \
            $(DPO)geometry.mod \
            $(DPO)constant.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_initial_conditions_dem.f  -o $(DPO)check_initial_conditions_dem.$(OBJ_EXT) -module $(DPO)
$(DPO)check_initial_conditions.$(OBJ_EXT) : ./check_data/check_initial_conditions.f \
            $(DPO)ic.mod \
            $(DPO)run.mod \
            $(DPO)param.mod \
            $(DPO)error_manager.mod \
            $(DPO)geometry.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)scalars.mod \
            $(DPO)discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_initial_conditions.f  -o $(DPO)check_initial_conditions.$(OBJ_EXT) -module $(DPO)
$(DPO)check_initial_conditions_mppic.$(OBJ_EXT) : ./check_data/check_initial_conditions_mppic.f \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)ic.mod \
            $(DPO)param1.mod \
            $(DPO)param.mod \
            $(DPO)geometry.mod \
            $(DPO)constant.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)error_manager.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_initial_conditions_mppic.f  -o $(DPO)check_initial_conditions_mppic.$(OBJ_EXT) -module $(DPO)
$(DPO)check_internal_surfaces.$(OBJ_EXT) : ./check_data/check_internal_surfaces.f \
            $(DPO)is.mod \
            $(DPO)param.mod \
            $(DPO)error_manager.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_internal_surfaces.f  -o $(DPO)check_internal_surfaces.$(OBJ_EXT) -module $(DPO)
$(DPO)check_numerics.$(OBJ_EXT) : ./check_data/check_numerics.f \
            $(DPO)run.mod \
            $(DPO)leqsol.mod \
            $(DPO)parallel.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_numerics.f  -o $(DPO)check_numerics.$(OBJ_EXT) -module $(DPO)
$(DPO)check_output_control.$(OBJ_EXT) : ./check_data/check_output_control.f \
            $(DPO)output.mod \
            $(DPO)run.mod \
            $(DPO)rxns.mod \
            $(DPO)param1.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_output_control.f  -o $(DPO)check_output_control.$(OBJ_EXT) -module $(DPO)
$(DPO)check_point_sources.$(OBJ_EXT) : ./check_data/check_point_sources.f \
            $(DPO)ps.mod \
            $(DPO)param.mod \
            $(DPO)error_manager.mod \
            $(DPO)geometry.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_point_sources.f  -o $(DPO)check_point_sources.$(OBJ_EXT) -module $(DPO)
$(DPO)check_run_control.$(OBJ_EXT) : ./check_data/check_run_control.f \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)scalars.mod \
            $(DPO)discretelement.mod \
            $(DPO)param1.mod \
            $(DPO)error_manager.mod \
            $(DPO)cutcell.mod \
            $(DPO)mfix_pic.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_run_control.f  -o $(DPO)check_run_control.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_common_all.$(OBJ_EXT) : ./check_data/check_solids_common_all.f \
            $(DPO)rxns.mod \
            $(DPO)physprop.mod \
            $(DPO)discretelement.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)error_manager.mod \
            $(DPO)run.mod \
            $(DPO)drag.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_common_all.f  -o $(DPO)check_solids_common_all.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_common_discrete.$(OBJ_EXT) : ./check_data/check_solids_common_discrete.f \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cutcell.mod \
            $(DPO)physprop.mod \
            $(DPO)desgrid.mod \
            $(DPO)run.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)error_manager.mod \
            $(DPO)compar.mod \
            $(DPO)des_thermo.mod \
            $(DPO)param1.mod \
            $(DPO)des_rxns.mod \
            $(DPO)geometry.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_common_discrete.f  -o $(DPO)check_solids_common_discrete.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_continuum.$(OBJ_EXT) : ./check_data/check_solids_continuum.f \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)scalars.mod \
            $(DPO)funits.mod \
            $(DPO)rxns.mod \
            $(DPO)cutcell.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_continuum.f  -o $(DPO)check_solids_continuum.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_des.$(OBJ_EXT) : ./check_data/check_solids_des.f \
            $(DPO)discretelement.mod \
            $(DPO)param1.mod \
            $(DPO)error_manager.mod \
            $(DPO)des_thermo.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_des.f  -o $(DPO)check_solids_des.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_model_prereqs.$(OBJ_EXT) : ./check_data/check_solids_model_prereqs.f \
            $(DPO)run.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)physprop.mod \
            $(DPO)param.mod \
            $(DPO)error_manager.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_model_prereqs.f  -o $(DPO)check_solids_model_prereqs.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_mppic.$(OBJ_EXT) : ./check_data/check_solids_mppic.f \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)funits.mod \
            $(DPO)discretelement.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)toleranc.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cutcell.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)error_manager.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_mppic.f  -o $(DPO)check_solids_mppic.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_phases.$(OBJ_EXT) : ./check_data/check_solids_phases.f \
            $(DPO)run.mod \
            $(DPO)error_manager.mod \
            $(DPO)discretelement.mod \
            $(DPO)physprop.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_phases.f  -o $(DPO)check_solids_phases.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_odepack.$(OBJ_EXT) : ./chem/check_data_odepack.f \
            $(DPO)funits.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)rxns.mod \
            $(DPO)stiff_chem.mod \
            $(DPO)compar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)stiff_chem_dbg.mod \
            $(DPO)stiff_chem_stats.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/check_data_odepack.f  -o $(DPO)check_data_odepack.$(OBJ_EXT) -module $(DPO)
$(DPO)stiff_chem_rrates.$(OBJ_EXT) : ./chem/stiff_chem_rrates.f \
            $(DPO)stiff_chem_maps.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)stiff_chem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/stiff_chem_rrates.f  -o $(DPO)stiff_chem_rrates.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_force_des_cutcell.$(OBJ_EXT) : ./des/calc_force_des_cutcell.f \
            $(DPO)run.mod \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)cutcell.mod \
            $(DPO)funits.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)parallel.mod \
            $(DPO)softspring_funcs_cutcell.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_force_des_cutcell.f  -o $(DPO)calc_force_des_cutcell.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_force_des.$(OBJ_EXT) : ./des/calc_force_des.f \
            $(DPO)run.mod \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_force_des.f  -o $(DPO)calc_force_des.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_rrate_des.$(OBJ_EXT) : ./des/calc_rrate_des.f \
            $(DPO)compar.mod \
            $(DPO)des_rxns.mod \
            $(DPO)discretelement.mod \
            $(DPO)energy.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)interpolation.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)rxns.mod \
            $(DPO)usr.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_rrate_des.f  -o $(DPO)calc_rrate_des.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_thermo_des.$(OBJ_EXT) : ./des/calc_thermo_des.f \
            $(DPO)compar.mod \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)interpolation.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_thermo_des.f  -o $(DPO)calc_thermo_des.$(OBJ_EXT) -module $(DPO)
$(DPO)cfassign.$(OBJ_EXT) : ./des/cfassign.f \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod \
            $(DPO)parallel.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)error_manager.mod \
            $(DPO)param.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)indices.mod \
            $(DPO)bc.mod \
            $(DPO)cutcell.mod \
            b_force1.inc                                                 \
            b_force2.inc                                                 \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfassign.f  -o $(DPO)cfassign.$(OBJ_EXT) -module $(DPO)
$(DPO)cffctowall.$(OBJ_EXT) : ./des/cffctowall.f \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffctowall.f  -o $(DPO)cffctowall.$(OBJ_EXT) -module $(DPO)
$(DPO)cffctow.$(OBJ_EXT) : ./des/cffctow.f \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffctow.f  -o $(DPO)cffctow.$(OBJ_EXT) -module $(DPO)
$(DPO)cfnewvalues.$(OBJ_EXT) : ./des/cfnewvalues.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)discretelement.mod \
            $(DPO)des_bc.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)mppic_wallbc.mod \
            $(DPO)cutcell.mod \
            $(DPO)randomno.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfnewvalues.f  -o $(DPO)cfnewvalues.$(OBJ_EXT) -module $(DPO)
$(DPO)cfrelvel.$(OBJ_EXT) : ./des/cfrelvel.f \
            $(DPO)discretelement.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfrelvel.f  -o $(DPO)cfrelvel.$(OBJ_EXT) -module $(DPO)
$(DPO)cfslide.$(OBJ_EXT) : ./des/cfslide.f \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslide.f  -o $(DPO)cfslide.$(OBJ_EXT) -module $(DPO)
$(DPO)cfslidewall.$(OBJ_EXT) : ./des/cfslidewall.f \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslidewall.f  -o $(DPO)cfslidewall.$(OBJ_EXT) -module $(DPO)
$(DPO)cfupdateold.$(OBJ_EXT) : ./des/cfupdateold.f \
            $(DPO)discretelement.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfupdateold.f  -o $(DPO)cfupdateold.$(OBJ_EXT) -module $(DPO)
$(DPO)cfwallcontact.$(OBJ_EXT) : ./des/cfwallcontact.f \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)parallel.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)des_bc.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallcontact.f  -o $(DPO)cfwallcontact.$(OBJ_EXT) -module $(DPO)
$(DPO)cfwallposvel.$(OBJ_EXT) : ./des/cfwallposvel.f \
            $(DPO)discretelement.mod \
            $(DPO)des_bc.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)matrix.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)drag.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallposvel.f  -o $(DPO)cfwallposvel.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_bc.$(OBJ_EXT) : ./des/check_des_bc.f \
            $(DPO)constant.mod \
            $(DPO)des_bc.mod \
            $(DPO)discretelement.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_bc.f  -o $(DPO)check_des_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_cohesion.$(OBJ_EXT) : ./des/check_des_cohesion.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_cohesion.f  -o $(DPO)check_des_cohesion.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_collision.$(OBJ_EXT) : ./des/check_des_collision.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_collision.f  -o $(DPO)check_des_collision.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_coupling.$(OBJ_EXT) : ./des/check_des_coupling.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_coupling.f  -o $(DPO)check_des_coupling.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_data.$(OBJ_EXT) : ./des/check_des_data.f \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)desgrid.mod \
            $(DPO)funits.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_data.f  -o $(DPO)check_des_data.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_energy.$(OBJ_EXT) : ./des/check_des_energy.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_energy.f  -o $(DPO)check_des_energy.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_geometry.$(OBJ_EXT) : ./des/check_des_geometry.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_geometry.f  -o $(DPO)check_des_geometry.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_hybrid.$(OBJ_EXT) : ./des/check_des_hybrid.f \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)funits.mod \
            $(DPO)param1.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_hybrid.f  -o $(DPO)check_des_hybrid.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_mppic.$(OBJ_EXT) : ./des/check_des_mppic.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_mppic.f  -o $(DPO)check_des_mppic.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_rxns.$(OBJ_EXT) : ./des/check_des_rxns.f \
            $(DPO)compar.mod \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)funits.mod \
            $(DPO)run.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parse.mod \
            $(DPO)physprop.mod \
            $(DPO)rxns.mod \
            $(DPO)stiff_chem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_rxns.f  -o $(DPO)check_des_rxns.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_thermo.$(OBJ_EXT) : ./des/check_des_thermo.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_thermo.f  -o $(DPO)check_des_thermo.$(OBJ_EXT) -module $(DPO)
$(DPO)des_allocate_arrays.$(OBJ_EXT) : ./des/des_allocate_arrays.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)discretelement.mod \
            $(DPO)indices.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)physprop.mod \
            $(DPO)des_bc.mod \
            $(DPO)funits.mod \
            $(DPO)desgrid.mod \
            $(DPO)desmpi.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod \
            $(DPO)cutcell.mod \
            $(DPO)error_manager.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_allocate_arrays.f  -o $(DPO)des_allocate_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)des_check_particle.$(OBJ_EXT) : ./des/des_check_particle.f \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)des_bc.mod \
            $(DPO)discretelement.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_check_particle.f  -o $(DPO)des_check_particle.$(OBJ_EXT) -module $(DPO)
$(DPO)des_cluster_identification.$(OBJ_EXT) : ./des/des_cluster_identification.f \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            $(DPO)des_cluster.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_cluster_identification.f  -o $(DPO)des_cluster_identification.$(OBJ_EXT) -module $(DPO)
$(DPO)des_functions.$(OBJ_EXT) : ./des/des_functions.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod \
            $(DPO)compar.mod \
            $(DPO)funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_functions.f  -o $(DPO)des_functions.$(OBJ_EXT) -module $(DPO)
$(DPO)des_granular_temperature.$(OBJ_EXT) : ./des/des_granular_temperature.f \
            $(DPO)discretelement.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)physprop.mod \
            $(DPO)des_bc.mod \
            $(DPO)fldvar.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_granular_temperature.f  -o $(DPO)des_granular_temperature.$(OBJ_EXT) -module $(DPO)
$(DPO)des_init_arrays.$(OBJ_EXT) : ./des/des_init_arrays.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod \
            $(DPO)indices.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)physprop.mod \
            $(DPO)des_bc.mod \
            $(DPO)run.mod \
            $(DPO)desgrid.mod \
            $(DPO)desmpi.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_arrays.f  -o $(DPO)des_init_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)des_init_bc.$(OBJ_EXT) : ./des/des_init_bc.f \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)des_bc.mod \
            $(DPO)discretelement.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_bc.f  -o $(DPO)des_init_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)des_init_namelist.$(OBJ_EXT) : ./des/des_init_namelist.f \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)des_bc.mod \
            $(DPO)des_ic.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod \
            des/desnamelist.inc                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_namelist.f  -o $(DPO)des_init_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)des_mass_inlet.$(OBJ_EXT) : ./des/des_mass_inlet.f \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)des_bc.mod \
            $(DPO)discretelement.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)desgrid.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod \
            function.inc                                                 \
            des/desgrid_functions.inc                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_mass_inlet.f  -o $(DPO)des_mass_inlet.$(OBJ_EXT) -module $(DPO)
$(DPO)des_physical_prop.$(OBJ_EXT) : ./des/des_physical_prop.f \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)funits.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_physical_prop.f  -o $(DPO)des_physical_prop.$(OBJ_EXT) -module $(DPO)
$(DPO)des_reaction_model.$(OBJ_EXT) : ./des/des_reaction_model.f \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)param1.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_reaction_model.f  -o $(DPO)des_reaction_model.$(OBJ_EXT) -module $(DPO)
$(DPO)des_rrates0.$(OBJ_EXT) : ./des/des_rrates0.f \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod \
            $(DPO)discretelement.mod \
            $(DPO)energy.mod \
            $(DPO)fldvar.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)parallel.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parse.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)rxns.mod \
            $(DPO)sendrecv.mod \
            $(DPO)usr.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_rrates0.f  -o $(DPO)des_rrates0.$(OBJ_EXT) -module $(DPO)
$(DPO)des_set_ic.$(OBJ_EXT) : ./des/des_set_ic.f \
            $(DPO)compar.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)des_ic.mod \
            $(DPO)des_rxns.mod \
            $(DPO)funits.mod \
            $(DPO)physprop.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_set_ic.f  -o $(DPO)des_set_ic.$(OBJ_EXT) -module $(DPO)
$(DPO)des_thermo_cond.$(OBJ_EXT) : ./des/des_thermo_cond.f \
            $(DPO)constant.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)funits.mod \
            $(DPO)physprop.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_cond.f  -o $(DPO)des_thermo_cond.$(OBJ_EXT) -module $(DPO)
$(DPO)des_thermo_conv.$(OBJ_EXT) : ./des/des_thermo_conv.f \
            $(DPO)constant.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)fldvar.mod \
            $(DPO)interpolation.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)usr.mod \
            $(DPO)compar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_conv.f  -o $(DPO)des_thermo_conv.$(OBJ_EXT) -module $(DPO)
$(DPO)des_thermo_newvalues.$(OBJ_EXT) : ./des/des_thermo_newvalues.f \
            $(DPO)compar.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod \
            $(DPO)discretelement.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_newvalues.f  -o $(DPO)des_thermo_newvalues.$(OBJ_EXT) -module $(DPO)
$(DPO)des_thermo_rad.$(OBJ_EXT) : ./des/des_thermo_rad.f \
            $(DPO)constant.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)fldvar.mod \
            $(DPO)param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_rad.f  -o $(DPO)des_thermo_rad.$(OBJ_EXT) -module $(DPO)
$(DPO)des_time_march.$(OBJ_EXT) : ./des/des_time_march.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)cont.mod \
            $(DPO)tau_g.mod \
            $(DPO)tau_s.mod \
            $(DPO)visc_g.mod \
            $(DPO)visc_s.mod \
            $(DPO)funits.mod \
            $(DPO)vshear.mod \
            $(DPO)scalars.mod \
            $(DPO)drag.mod \
            $(DPO)rxns.mod \
            $(DPO)compar.mod \
            $(DPO)time_cpu.mod \
            $(DPO)discretelement.mod \
            $(DPO)constant.mod \
            $(DPO)sendrecv.mod \
            $(DPO)des_bc.mod \
            $(DPO)cutcell.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod \
            $(DPO)interpolation.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_time_march.f  -o $(DPO)des_time_march.$(OBJ_EXT) -module $(DPO)
$(DPO)des_wallbc_preprocessing.$(OBJ_EXT) : ./des/des_wallbc_preprocessing.f \
            $(DPO)param1.mod \
            $(DPO)funits.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)cutcell.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)bc.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)sendrecv.mod \
            $(DPO)error_manager.mod \
            $(DPO)softspring_funcs_cutcell.mod \
            $(DPO)desmpi.mod \
            $(DPO)cdist.mod \
            $(DPO)discretelement.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_wallbc_preprocessing.f  -o $(DPO)des_wallbc_preprocessing.$(OBJ_EXT) -module $(DPO)
$(DPO)drag_fgs.$(OBJ_EXT) : ./des/drag_fgs.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)cutcell.mod \
            $(DPO)interpolation.mod \
            $(DPO)desmpi.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)constant.mod \
            $(DPO)drag.mod \
            $(DPO)ur_facs.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/drag_fgs.f  -o $(DPO)drag_fgs.$(OBJ_EXT) -module $(DPO)
$(DPO)gas_drag.$(OBJ_EXT) : ./des/gas_drag.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_g.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)drag.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/gas_drag.f  -o $(DPO)gas_drag.$(OBJ_EXT) -module $(DPO)
$(DPO)generate_particle_config.$(OBJ_EXT) : ./des/generate_particle_config.f \
            $(DPO)mfix_pic.mod \
            $(DPO)des_linked_list_data.mod \
            $(DPO)des_linked_list_funcs.mod \
            $(DPO)cutcell.mod \
            $(DPO)discretelement.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)error_manager.mod \
            $(DPO)des_linked_list_data.mod \
            $(DPO)softspring_funcs_cutcell.mod \
            $(DPO)indices.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)ic.mod \
            $(DPO)param1.mod \
            $(DPO)param.mod \
            $(DPO)randomno.mod \
            $(DPO)funits.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)parallel.mod \
            $(DPO)desmpi.mod \
            $(DPO)cdist.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/generate_particle_config.f  -o $(DPO)generate_particle_config.$(OBJ_EXT) -module $(DPO)
$(DPO)make_arrays_des.$(OBJ_EXT) : ./des/make_arrays_des.f \
            $(DPO)param1.mod \
            $(DPO)funits.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)cutcell.mod \
            $(DPO)desmpi.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)geometry.mod \
            $(DPO)des_ic.mod \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_stl_functions.mod \
            $(DPO)cdist.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/make_arrays_des.f  -o $(DPO)make_arrays_des.$(OBJ_EXT) -module $(DPO)
$(DPO)mppic_routines.$(OBJ_EXT) : ./des/mppic_routines.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)pgcor.mod \
            $(DPO)pscor.mod \
            $(DPO)cont.mod \
            $(DPO)tau_g.mod \
            $(DPO)tau_s.mod \
            $(DPO)visc_g.mod \
            $(DPO)visc_s.mod \
            $(DPO)funits.mod \
            $(DPO)vshear.mod \
            $(DPO)scalars.mod \
            $(DPO)drag.mod \
            $(DPO)rxns.mod \
            $(DPO)compar.mod \
            $(DPO)time_cpu.mod \
            $(DPO)discretelement.mod \
            $(DPO)constant.mod \
            $(DPO)sendrecv.mod \
            $(DPO)des_bc.mod \
            $(DPO)cutcell.mod \
            $(DPO)mppic_wallbc.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod \
            $(DPO)interpolation.mod \
            $(DPO)parallel.mod \
            $(DPO)indices.mod \
            $(DPO)desmpi.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)bc.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)toleranc.mod \
            $(DPO)is.mod \
            $(DPO)quadric.mod \
            $(DPO)vtk.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mppic_routines.f  -o $(DPO)mppic_routines.$(OBJ_EXT) -module $(DPO)
$(DPO)neighbour.$(OBJ_EXT) : ./des/neighbour.f \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod \
            $(DPO)desgrid.mod \
            $(DPO)des_thermo.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/neighbour.f  -o $(DPO)neighbour.$(OBJ_EXT) -module $(DPO)
$(DPO)nsquare.$(OBJ_EXT) : ./des/nsquare.f \
            $(DPO)param1.mod \
            $(DPO)discretelement.mod \
            $(DPO)geometry.mod \
            $(DPO)des_bc.mod \
            $(DPO)des_thermo.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/nsquare.f  -o $(DPO)nsquare.$(OBJ_EXT) -module $(DPO)
$(DPO)particles_in_cell.$(OBJ_EXT) : ./des/particles_in_cell.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)compar.mod \
            $(DPO)parallel.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)desgrid.mod \
            $(DPO)desmpi.mod \
            $(DPO)cutcell.mod \
            $(DPO)mfix_pic.mod \
            $(DPO)des_rxns.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)bc.mod \
            $(DPO)drag.mod \
            $(DPO)interpolation.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)funits.mod \
            $(DPO)desmpi_wrapper.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/particles_in_cell.f  -o $(DPO)particles_in_cell.$(OBJ_EXT) -module $(DPO)
$(DPO)read_des_restart.$(OBJ_EXT) : ./des/read_des_restart.f \
            $(DPO)param1.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            $(DPO)des_bc.mod \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)desmpi.mod \
            $(DPO)machine.mod \
            $(DPO)cdist.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)constant.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/read_des_restart.f  -o $(DPO)read_des_restart.$(OBJ_EXT) -module $(DPO)
$(DPO)solid_drag.$(OBJ_EXT) : ./des/solid_drag.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)matrix.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/solid_drag.f  -o $(DPO)solid_drag.$(OBJ_EXT) -module $(DPO)
$(DPO)usr0_des.$(OBJ_EXT) : ./des/usr0_des.f \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            $(DPO)usr.mod \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr0_des.f  -o $(DPO)usr0_des.$(OBJ_EXT) -module $(DPO)
$(DPO)usr1_des.$(OBJ_EXT) : ./des/usr1_des.f \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            $(DPO)usr.mod \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr1_des.f  -o $(DPO)usr1_des.$(OBJ_EXT) -module $(DPO)
$(DPO)usr2_des.$(OBJ_EXT) : ./des/usr2_des.f \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            $(DPO)usr.mod \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr2_des.f  -o $(DPO)usr2_des.$(OBJ_EXT) -module $(DPO)
$(DPO)usr3_des.$(OBJ_EXT) : ./des/usr3_des.f \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            $(DPO)usr.mod \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr3_des.f  -o $(DPO)usr3_des.$(OBJ_EXT) -module $(DPO)
$(DPO)usr4_des.$(OBJ_EXT) : ./des/usr4_des.f \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)discretelement.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)usr.mod \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr4_des.f  -o $(DPO)usr4_des.$(OBJ_EXT) -module $(DPO)
$(DPO)usr_rates_des.$(OBJ_EXT) : ./des/usr_rates_des.f \
            $(DPO)compar.mod \
            $(DPO)constant.mod \
            $(DPO)des_thermo.mod \
            $(DPO)des_rxns.mod \
            $(DPO)discretelement.mod \
            $(DPO)energy.mod \
            $(DPO)fldvar.mod \
            $(DPO)funits.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)physprop.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)usr.mod \
            species.inc                                                  \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr_rates_des.f  -o $(DPO)usr_rates_des.$(OBJ_EXT) -module $(DPO)
$(DPO)walledgecontact.$(OBJ_EXT) : ./des/walledgecontact.f \
            $(DPO)discretelement.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)matrix.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)drag.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/walledgecontact.f  -o $(DPO)walledgecontact.$(OBJ_EXT) -module $(DPO)
$(DPO)wallfacecontact.$(OBJ_EXT) : ./des/wallfacecontact.f \
            $(DPO)discretelement.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)matrix.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)drag.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/wallfacecontact.f  -o $(DPO)wallfacecontact.$(OBJ_EXT) -module $(DPO)
$(DPO)wallnodecontact.$(OBJ_EXT) : ./des/wallnodecontact.f \
            $(DPO)discretelement.mod \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)matrix.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)drag.mod \
            $(DPO)constant.mod \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/wallnodecontact.f  -o $(DPO)wallnodecontact.$(OBJ_EXT) -module $(DPO)
$(DPO)write_des_data.$(OBJ_EXT) : ./des/write_des_data.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)sendrecv.mod \
            $(DPO)des_bc.mod \
            $(DPO)mpi_utility.mod \
            $(DPO)compar.mod \
            $(DPO)desmpi.mod \
            $(DPO)cdist.mod \
            $(DPO)des_thermo.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/write_des_data.f  -o $(DPO)write_des_data.$(OBJ_EXT) -module $(DPO)
$(DPO)write_des_restart.$(OBJ_EXT) : ./des/write_des_restart.f \
            $(DPO)param1.mod \
            $(DPO)compar.mod \
            $(DPO)discretelement.mod \
            $(DPO)run.mod \
            $(DPO)des_bc.mod \
            $(DPO)des_rxns.mod \
            $(DPO)des_thermo.mod \
            $(DPO)desmpi.mod \
            $(DPO)machine.mod \
            $(DPO)cdist.mod \
            $(DPO)mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/write_des_restart.f  -o $(DPO)write_des_restart.$(OBJ_EXT) -module $(DPO)
$(DPO)gaussj.$(OBJ_EXT) : ./dqmom/gaussj.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/gaussj.f  -o $(DPO)gaussj.$(OBJ_EXT) -module $(DPO)
$(DPO)odeint.$(OBJ_EXT) : ./dqmom/odeint.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/odeint.f  -o $(DPO)odeint.$(OBJ_EXT) -module $(DPO)
$(DPO)rkck.$(OBJ_EXT) : ./dqmom/rkck.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/rkck.f  -o $(DPO)rkck.$(OBJ_EXT) -module $(DPO)
$(DPO)rkqs.$(OBJ_EXT) : ./dqmom/rkqs.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/rkqs.f  -o $(DPO)rkqs.$(OBJ_EXT) -module $(DPO)
$(DPO)source_population_eq.$(OBJ_EXT) : ./dqmom/source_population_eq.f \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)scalars.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/source_population_eq.f  -o $(DPO)source_population_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)usr_dqmom.$(OBJ_EXT) : ./dqmom/usr_dqmom.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)run.mod \
            $(DPO)physprop.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)output.mod \
            $(DPO)indices.mod \
            $(DPO)rxns.mod \
            $(DPO)constant.mod \
            $(DPO)ambm.mod \
            $(DPO)compar.mod \
            $(DPO)scalars.mod \
            $(DPO)usr.mod \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/usr_dqmom.f  -o $(DPO)usr_dqmom.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_eps_ghd.$(OBJ_EXT) : ./GhdTheory/adjust_eps_ghd.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)toleranc.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)ghdtheory.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/adjust_eps_ghd.f  -o $(DPO)adjust_eps_ghd.$(OBJ_EXT) -module $(DPO)
$(DPO)bulk_viscosity.$(OBJ_EXT) : ./GhdTheory/bulk_viscosity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/bulk_viscosity.f  -o $(DPO)bulk_viscosity.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_d_ghd.$(OBJ_EXT) : ./GhdTheory/calc_d_ghd.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)scales.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_d_ghd.f  -o $(DPO)calc_d_ghd.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_external_forces.$(OBJ_EXT) : ./GhdTheory/calc_external_forces.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)drag.mod \
            $(DPO)bc.mod \
            $(DPO)scales.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            b_force1.inc                                                 \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_external_forces.f  -o $(DPO)calc_external_forces.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_nflux.$(OBJ_EXT) : ./GhdTheory/calc_nflux.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)constant.mod \
            $(DPO)indices.mod \
            $(DPO)mflux.mod \
            $(DPO)compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_nflux.f  -o $(DPO)calc_nflux.$(OBJ_EXT) -module $(DPO)
$(DPO)chi_ij_GHD.$(OBJ_EXT) : ./GhdTheory/chi_ij_GHD.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/chi_ij_GHD.f  -o $(DPO)chi_ij_GHD.$(OBJ_EXT) -module $(DPO)
$(DPO)cooling_rate.$(OBJ_EXT) : ./GhdTheory/cooling_rate.f \
            $(DPO)compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/cooling_rate.f  -o $(DPO)cooling_rate.$(OBJ_EXT) -module $(DPO)
$(DPO)cooling_rate_tc.$(OBJ_EXT) : ./GhdTheory/cooling_rate_tc.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/cooling_rate_tc.f  -o $(DPO)cooling_rate_tc.$(OBJ_EXT) -module $(DPO)
$(DPO)dufour_coeff.$(OBJ_EXT) : ./GhdTheory/dufour_coeff.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/dufour_coeff.f  -o $(DPO)dufour_coeff.$(OBJ_EXT) -module $(DPO)
$(DPO)ghd.$(OBJ_EXT) : ./GhdTheory/ghd.f \
            $(DPO)drag.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ghd.f  -o $(DPO)ghd.$(OBJ_EXT) -module $(DPO)
$(DPO)ghdmassflux.$(OBJ_EXT) : ./GhdTheory/ghdmassflux.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)visc_s.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)drag.mod \
            $(DPO)is.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ghdmassflux.f  -o $(DPO)ghdmassflux.$(OBJ_EXT) -module $(DPO)
$(DPO)mass_mobility.$(OBJ_EXT) : ./GhdTheory/mass_mobility.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/mass_mobility.f  -o $(DPO)mass_mobility.$(OBJ_EXT) -module $(DPO)
$(DPO)ordinary_diff.$(OBJ_EXT) : ./GhdTheory/ordinary_diff.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ordinary_diff.f  -o $(DPO)ordinary_diff.$(OBJ_EXT) -module $(DPO)
$(DPO)partial_elim_ghd.$(OBJ_EXT) : ./GhdTheory/partial_elim_ghd.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)geometry.mod \
            $(DPO)matrix.mod \
            $(DPO)physprop.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)drag.mod \
            $(DPO)fldvar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/partial_elim_ghd.f  -o $(DPO)partial_elim_ghd.$(OBJ_EXT) -module $(DPO)
$(DPO)pressure.$(OBJ_EXT) : ./GhdTheory/pressure.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/pressure.f  -o $(DPO)pressure.$(OBJ_EXT) -module $(DPO)
$(DPO)shear_viscosity.$(OBJ_EXT) : ./GhdTheory/shear_viscosity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/shear_viscosity.f  -o $(DPO)shear_viscosity.$(OBJ_EXT) -module $(DPO)
$(DPO)source_ghd_granular_energy.$(OBJ_EXT) : ./GhdTheory/source_ghd_granular_energy.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)drag.mod \
            $(DPO)geometry.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_s.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)trace.mod \
            $(DPO)indices.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            $(DPO)compar.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            b_force1.inc                                                 \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/source_ghd_granular_energy.f  -o $(DPO)source_ghd_granular_energy.$(OBJ_EXT) -module $(DPO)
$(DPO)thermal_conductivity.$(OBJ_EXT) : ./GhdTheory/thermal_conductivity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/thermal_conductivity.f  -o $(DPO)thermal_conductivity.$(OBJ_EXT) -module $(DPO)
$(DPO)thermal_diffusivity.$(OBJ_EXT) : ./GhdTheory/thermal_diffusivity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/thermal_diffusivity.f  -o $(DPO)thermal_diffusivity.$(OBJ_EXT) -module $(DPO)
$(DPO)thermal_mobility.$(OBJ_EXT) : ./GhdTheory/thermal_mobility.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/thermal_mobility.f  -o $(DPO)thermal_mobility.$(OBJ_EXT) -module $(DPO)
$(DPO)transport_coeff_ghd.$(OBJ_EXT) : ./GhdTheory/transport_coeff_ghd.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)fldvar.mod \
            $(DPO)indices.mod \
            $(DPO)visc_s.mod \
            $(DPO)ghdtheory.mod \
            $(DPO)physprop.mod \
            $(DPO)run.mod \
            $(DPO)constant.mod \
            $(DPO)toleranc.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/transport_coeff_ghd.f  -o $(DPO)transport_coeff_ghd.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_allocate_arrays.$(OBJ_EXT) : ./qmomk/qmomk_allocate_arrays.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)indices.mod \
            $(DPO)geometry.mod \
            $(DPO)physprop.mod \
            $(DPO)qmom_kinetic_equation.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_allocate_arrays.f  -o $(DPO)qmomk_allocate_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_gas_drag.$(OBJ_EXT) : ./qmomk/qmomk_gas_drag.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)parallel.mod \
            $(DPO)matrix.mod \
            $(DPO)scales.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)visc_g.mod \
            $(DPO)rxns.mod \
            $(DPO)run.mod \
            $(DPO)toleranc.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)is.mod \
            $(DPO)tau_g.mod \
            $(DPO)bc.mod \
            $(DPO)compar.mod \
            $(DPO)sendrecv.mod \
            $(DPO)discretelement.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)drag.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_gas_drag.f  -o $(DPO)qmomk_gas_drag.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_init_bc.$(OBJ_EXT) : ./qmomk/qmomk_init_bc.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)bc.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)qmomk_quadrature.mod \
            $(DPO)qmomk_bc.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_init_bc.f  -o $(DPO)qmomk_init_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_initial_conditions.$(OBJ_EXT) : ./qmomk/qmomk_initial_conditions.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)qmomk_quadrature.mod \
            $(DPO)qmomk_parameters.mod \
            $(DPO)ic.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_initial_conditions.f  -o $(DPO)qmomk_initial_conditions.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_init_namelist.$(OBJ_EXT) : ./qmomk/qmomk_init_namelist.f \
            $(DPO)param1.mod \
            $(DPO)qmom_kinetic_equation.mod \
            qmomk/qmomknamelist.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_init_namelist.f  -o $(DPO)qmomk_init_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_make_arrays.$(OBJ_EXT) : ./qmomk/qmomk_make_arrays.f \
            $(DPO)param1.mod \
            $(DPO)geometry.mod \
            $(DPO)funits.mod \
            $(DPO)compar.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_make_arrays.f  -o $(DPO)qmomk_make_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_read_restart.$(OBJ_EXT) : ./qmomk/qmomk_read_restart.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)fldvar.mod \
            $(DPO)cont.mod \
            $(DPO)geometry.mod \
            $(DPO)indices.mod \
            $(DPO)run.mod \
            $(DPO)compar.mod \
            $(DPO)physprop.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)qmomk_quadrature.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_read_restart.f  -o $(DPO)qmomk_read_restart.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_set_bc.$(OBJ_EXT) : ./qmomk/qmomk_set_bc.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)compar.mod \
            $(DPO)indices.mod \
            $(DPO)bc.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)qmomk_quadrature.mod \
            $(DPO)qmomk_bc.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_set_bc.f  -o $(DPO)qmomk_set_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_time_march.$(OBJ_EXT) : ./qmomk/qmomk_time_march.f \
            $(DPO)param.mod \
            $(DPO)param1.mod \
            $(DPO)constant.mod \
            $(DPO)run.mod \
            $(DPO)output.mod \
            $(DPO)physprop.mod \
            $(DPO)fldvar.mod \
            $(DPO)geometry.mod \
            $(DPO)cont.mod \
            $(DPO)tau_g.mod \
            $(DPO)tau_s.mod \
            $(DPO)visc_g.mod \
            $(DPO)visc_s.mod \
            $(DPO)funits.mod \
            $(DPO)vshear.mod \
            $(DPO)scalars.mod \
            $(DPO)drag.mod \
            $(DPO)rxns.mod \
            $(DPO)compar.mod \
            $(DPO)time_cpu.mod \
            $(DPO)is.mod \
            $(DPO)indices.mod \
            $(DPO)matrix.mod \
            $(DPO)sendrecv.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)qmomk_fluxes.mod \
            $(DPO)qmomk_quadrature.mod \
            $(DPO)qmomk_collision.mod \
            $(DPO)qmomk_parameters.mod \
            $(DPO)ur_facs.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_time_march.f  -o $(DPO)qmomk_time_march.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_write_restart.$(OBJ_EXT) : ./qmomk/qmomk_write_restart.f \
            $(DPO)param1.mod \
            $(DPO)qmom_kinetic_equation.mod \
            $(DPO)run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_write_restart.f  -o $(DPO)qmomk_write_restart.$(OBJ_EXT) -module $(DPO)
$(DPO)get_values.$(OBJ_EXT) : ./thermochemical/get_values.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./thermochemical/get_values.f  -o $(DPO)get_values.$(OBJ_EXT) -module $(DPO)
$(DPO)readTherm.$(OBJ_EXT) : ./thermochemical/readTherm.f \
            $(DPO)physprop.mod \
            $(DPO)des_rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./thermochemical/readTherm.f  -o $(DPO)readTherm.$(OBJ_EXT) -module $(DPO)
