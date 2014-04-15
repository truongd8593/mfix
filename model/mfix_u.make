.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
$(EXEC_FILE) : \
    $(DPO)AMBM.mod \
    $(DPO)BC.mod \
    $(DPO)BOUNDFUNIJK3.mod \
    $(DPO)BOUNDFUNIJK.mod \
    $(DPO)CDIST.mod \
    $(DPO)CHECK.mod \
    $(DPO)CHISCHEME.mod \
    $(DPO)COEFF.mod \
    $(DPO)CONSTANT.mod \
    $(DPO)CONT.mod \
    $(DPO)CORNER.mod \
    $(DPO)DBG.mod \
    $(DPO)DRAG.mod \
    $(DPO)ENERGY.mod \
    $(DPO)ERROR_MANAGER.mod \
    $(DPO)FLDVAR.mod \
    $(DPO)FUNCTION.mod \
    $(DPO)FUNITS.mod \
    $(DPO)GEOMETRY.mod \
    $(DPO)IC.mod \
    $(DPO)INDICES.mod \
    $(DPO)IS.mod \
    $(DPO)KINTHEORY.mod \
    $(DPO)LEQSOL.mod \
    $(DPO)MACHINE.mod \
    $(DPO)MATRIX.mod \
    $(DPO)MFIX_NETCDF.mod \
    $(DPO)MFLUX.mod \
    $(DPO)MMS.mod \
    $(DPO)OUTPUT.mod \
    $(DPO)PARALLEL.mod \
    $(DPO)PARAM1.mod \
    $(DPO)PARAM.mod \
    $(DPO)PARSE.mod \
    $(DPO)PGCOR.mod \
    $(DPO)PHYSPROP.mod \
    $(DPO)PSCOR.mod \
    $(DPO)PS.mod \
    $(DPO)RESIDUAL.mod \
    $(DPO)RUN.mod \
    $(DPO)RXN_COM.mod \
    $(DPO)RXNS.mod \
    $(DPO)SCALARS.mod \
    $(DPO)SCALES.mod \
    $(DPO)TAU_G.mod \
    $(DPO)TAU_S.mod \
    $(DPO)TIME_CPU.mod \
    $(DPO)TMP_ARRAY1.mod \
    $(DPO)TMP_ARRAY.mod \
    $(DPO)TOLERANC.mod \
    $(DPO)TRACE.mod \
    $(DPO)TURB.mod \
    $(DPO)UR_FACS.mod \
    $(DPO)USR.mod \
    $(DPO)VISC_G.mod \
    $(DPO)VISC_S.mod \
    $(DPO)VSHEAR.mod \
    $(DPO)XSI_ARRAY.mod \
    $(DPO)CUTCELL.mod \
    $(DPO)DASHBOARD.mod \
    $(DPO)POLYGON.mod \
    $(DPO)PROGRESS_BAR.mod \
    $(DPO)QUADRIC.mod \
    $(DPO)STL.mod \
    $(DPO)VTK.mod \
    $(DPO)STIFF_CHEM_DBG.mod \
    $(DPO)STIFF_CHEM_MAPS.mod \
    $(DPO)STIFF_CHEM.mod \
    $(DPO)STIFF_CHEM_STATS.mod \
    $(DPO)DES_BC.mod \
    $(DPO)DES_CLUSTER.mod \
    $(DPO)DESGRID.mod \
    $(DPO)DES_IC.mod \
    $(DPO)DES_LINKED_LIST_DATA.mod \
    $(DPO)DES_LINKED_LIST_FUNCS.mod \
    $(DPO)DESMPI.mod \
    $(DPO)DESMPI_WRAPPER.mod \
    $(DPO)DES_RXNS.mod \
    $(DPO)DES_STL_FUNCTIONS.mod \
    $(DPO)DES_THERMO.mod \
    $(DPO)DISCRETELEMENT.mod \
    $(DPO)INTERPOLATION.mod \
    $(DPO)MFIX_PIC.mod \
    $(DPO)MPPIC_WALLBC.mod \
    $(DPO)RANDOMNO.mod \
    $(DPO)SENDRECVNODE.mod \
    $(DPO)SOFTSPRING_FUNCS_CUTCELL.mod \
    $(DPO)COMPAR.mod \
    $(DPO)DBG_UTIL.mod \
    $(DPO)DEBUG.mod \
    $(DPO)GRIDMAP.mod \
    $(DPO)MPI.mod \
    $(DPO)MPI_UTILITY.mod \
    $(DPO)PARALLEL_MPI.mod \
    $(DPO)SENDRECV3.mod \
    $(DPO)SENDRECV.mod \
    $(DPO)GHDTHEORY.mod \
    $(DPO)QMOMK_BC.mod \
    $(DPO)QMOMK_COLLISION.mod \
    $(DPO)QMOMK_FLUXES.mod \
    $(DPO)QMOM_KINETIC_EQUATION.mod \
    $(DPO)QMOMK_PARAMETERS.mod \
    $(DPO)QMOMK_QUADRATURE.mod \
    $(DPO)QMOMK_TOOLS.mod \
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
    $(DPO)check_chemical_rxns.$(OBJ_EXT) \
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
    $(DPO)calc_force_dem.$(OBJ_EXT) \
    $(DPO)calc_force_dem_stl.$(OBJ_EXT) \
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
    $(DPO)check_des_data.$(OBJ_EXT) \
    $(DPO)check_des_hybrid.$(OBJ_EXT) \
    $(DPO)check_des_rxns.$(OBJ_EXT) \
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
    $(DPO)check_chemical_rxns.$(OBJ_EXT) \
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
    $(DPO)calc_force_dem.$(OBJ_EXT) \
    $(DPO)calc_force_dem_stl.$(OBJ_EXT) \
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
    $(DPO)check_des_data.$(OBJ_EXT) \
    $(DPO)check_des_hybrid.$(OBJ_EXT) \
    $(DPO)check_des_rxns.$(OBJ_EXT) \
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
$(DPO)AMBM.mod : ambm_mod.f \
            $(DPO)COMPAR.mod \
            $(DPO)FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ambm_mod.f  -o $(DPO)ambm_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)BC.mod : bc_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) bc_mod.f  -o $(DPO)bc_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)BOUNDFUNIJK3.mod : boundfunijk3_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk3_mod.f  -o $(DPO)boundfunijk3_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)BOUNDFUNIJK.mod : boundfunijk_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk_mod.f  -o $(DPO)boundfunijk_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)CDIST.mod : cdist_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) cdist_mod.f  -o $(DPO)cdist_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)CHECK.mod : check_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_mod.f  -o $(DPO)check_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)CHISCHEME.mod : chischeme_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) chischeme_mod.f  -o $(DPO)chischeme_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)COEFF.mod : coeff_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)VISC_G.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MMS.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) coeff_mod.f  -o $(DPO)coeff_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)CONSTANT.mod : constant_mod.f \
            $(DPO)PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) constant_mod.f  -o $(DPO)constant_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)CONT.mod : cont_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) cont_mod.f  -o $(DPO)cont_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)CORNER.mod : corner_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) corner_mod.f  -o $(DPO)corner_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DBG.mod : dbg_mod.f \
            $(DPO)PARAM1.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)MFLUX.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PARAM.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)SENDRECV.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) dbg_mod.f  -o $(DPO)dbg_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DRAG.mod : drag_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) drag_mod.f  -o $(DPO)drag_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)ENERGY.mod : energy_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) energy_mod.f  -o $(DPO)energy_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)ERROR_MANAGER.mod : error_manager_mod.f \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) error_manager_mod.f  -o $(DPO)error_manager_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)FLDVAR.mod : fldvar_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) fldvar_mod.f  -o $(DPO)fldvar_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)FUNCTION.mod : function_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) function_mod.f  -o $(DPO)function_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)FUNITS.mod : funits_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) funits_mod.f  -o $(DPO)funits_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)GEOMETRY.mod : geometry_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) geometry_mod.f  -o $(DPO)geometry_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)IC.mod : ic_mod.f \
            $(DPO)PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ic_mod.f  -o $(DPO)ic_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)INDICES.mod : indices_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) indices_mod.f  -o $(DPO)indices_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)IS.mod : is_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) is_mod.f  -o $(DPO)is_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)KINTHEORY.mod : kintheory_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_mod.f  -o $(DPO)kintheory_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)LEQSOL.mod : leqsol_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) leqsol_mod.f  -o $(DPO)leqsol_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)MACHINE.mod : machine_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) machine_mod.f  -o $(DPO)machine_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)MATRIX.mod : matrix_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) matrix_mod.f  -o $(DPO)matrix_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)MFIX_NETCDF.mod : mfix_netcdf_mod.f \
            MFIX_netcdf_constants.fi                                     \
            MFIX_netcdf_overloads.fi                                     \
            MFIX_netcdf_variables.fi                                     \
            MFIX_netcdf_misc.fi                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) mfix_netcdf_mod.f  -o $(DPO)mfix_netcdf_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)MFLUX.mod : mflux_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) mflux_mod.f  -o $(DPO)mflux_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)MMS.mod : mms_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)FLDVAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) mms_mod.f  -o $(DPO)mms_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)OUTPUT.mod : output_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) output_mod.f  -o $(DPO)output_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)PARALLEL.mod : parallel_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parallel_mod.f  -o $(DPO)parallel_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)PARAM1.mod : param1_mod.f \
            $(DPO)PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) param1_mod.f  -o $(DPO)param1_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)PARAM.mod : param_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) param_mod.f  -o $(DPO)param_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)PARSE.mod : parse_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)RXN_COM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_mod.f  -o $(DPO)parse_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)PGCOR.mod : pgcor_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pgcor_mod.f  -o $(DPO)pgcor_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)PHYSPROP.mod : physprop_mod.f \
            $(DPO)PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) physprop_mod.f  -o $(DPO)physprop_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)PSCOR.mod : pscor_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pscor_mod.f  -o $(DPO)pscor_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)PS.mod : ps_mod.f \
            $(DPO)PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ps_mod.f  -o $(DPO)ps_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)RESIDUAL.mod : residual_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) residual_mod.f  -o $(DPO)residual_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)RUN.mod : run_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) run_mod.f  -o $(DPO)run_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)RXN_COM.mod : rxn_com_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FUNITS.mod \
            $(DPO)ERROR_MANAGER.mod \
            mfix_directory_path.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) rxn_com_mod.f  -o $(DPO)rxn_com_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)RXNS.mod : rxns_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RXN_COM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) rxns_mod.f  -o $(DPO)rxns_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)SCALARS.mod : scalars_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scalars_mod.f  -o $(DPO)scalars_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)SCALES.mod : scales_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scales_mod.f  -o $(DPO)scales_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)TAU_G.mod : tau_g_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_g_mod.f  -o $(DPO)tau_g_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)TAU_S.mod : tau_s_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_s_mod.f  -o $(DPO)tau_s_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)TIME_CPU.mod : time_cpu_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) time_cpu_mod.f  -o $(DPO)time_cpu_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)TMP_ARRAY1.mod : tmp_array1_mod.f \
            $(DPO)COMPAR.mod \
            $(DPO)FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array1_mod.f  -o $(DPO)tmp_array1_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)TMP_ARRAY.mod : tmp_array_mod.f \
            $(DPO)COMPAR.mod \
            $(DPO)FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array_mod.f  -o $(DPO)tmp_array_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)TOLERANC.mod : toleranc_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) toleranc_mod.f  -o $(DPO)toleranc_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)TRACE.mod : trace_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) trace_mod.f  -o $(DPO)trace_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)TURB.mod : turb_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) turb_mod.f  -o $(DPO)turb_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)UR_FACS.mod : ur_facs_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ur_facs_mod.f  -o $(DPO)ur_facs_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)USR.mod : usr_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_mod.f  -o $(DPO)usr_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)VISC_G.mod : visc_g_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_g_mod.f  -o $(DPO)visc_g_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)VISC_S.mod : visc_s_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_s_mod.f  -o $(DPO)visc_s_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)VSHEAR.mod : vshear_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) vshear_mod.f  -o $(DPO)vshear_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)XSI_ARRAY.mod : xsi_array_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) xsi_array_mod.f  -o $(DPO)xsi_array_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)CUTCELL.mod : ./cartesian_grid/cutcell_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PROGRESS_BAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cutcell_mod.f  -o $(DPO)cutcell_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DASHBOARD.mod : ./cartesian_grid/dashboard_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/dashboard_mod.f  -o $(DPO)dashboard_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)POLYGON.mod : ./cartesian_grid/polygon_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/polygon_mod.f  -o $(DPO)polygon_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)PROGRESS_BAR.mod : ./cartesian_grid/progress_bar_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/progress_bar_mod.f  -o $(DPO)progress_bar_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)QUADRIC.mod : ./cartesian_grid/quadric_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/quadric_mod.f  -o $(DPO)quadric_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)STL.mod : ./cartesian_grid/stl_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/stl_mod.f  -o $(DPO)stl_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)VTK.mod : ./cartesian_grid/vtk_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/vtk_mod.f  -o $(DPO)vtk_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)STIFF_CHEM_DBG.mod : ./chem/stiff_chem_dbg_mod.f \
            $(DPO)COMPAR.mod \
            $(DPO)STIFF_CHEM_STATS.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)RXNS.mod \
            $(DPO)INDICES.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/stiff_chem_dbg_mod.f  -o $(DPO)stiff_chem_dbg_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)STIFF_CHEM_MAPS.mod : ./chem/stiff_chem_maps_mod.f \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)RUN.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/stiff_chem_maps_mod.f  -o $(DPO)stiff_chem_maps_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)STIFF_CHEM.mod : ./chem/stiff_chem_mod.f \
            $(DPO)STIFF_CHEM_MAPS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)RUN.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)STIFF_CHEM_DBG.mod \
            $(DPO)STIFF_CHEM_STATS.mod \
            $(DPO)RXNS.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/stiff_chem_mod.f  -o $(DPO)stiff_chem_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)STIFF_CHEM_STATS.mod : ./chem/stiff_chem_stats_mod.f \
            $(DPO)COMPAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/stiff_chem_stats_mod.f  -o $(DPO)stiff_chem_stats_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DES_BC.mod : ./des/des_bc_mod.f \
            $(DPO)PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_bc_mod.f  -o $(DPO)des_bc_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DES_CLUSTER.mod : ./des/des_cluster_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)DESMPI_WRAPPER.mod \
            $(DPO)DESMPI.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_cluster_mod.f  -o $(DPO)des_cluster_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DESGRID.mod : ./des/desgrid_mod.f \
            $(DPO)PARAM1.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DESMPI_WRAPPER.mod \
            $(DPO)DES_THERMO.mod \
            des/desgrid_functions.inc                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desgrid_mod.f  -o $(DPO)desgrid_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DES_IC.mod : ./des/des_ic_mod.f \
            $(DPO)PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_ic_mod.f  -o $(DPO)des_ic_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DES_LINKED_LIST_DATA.mod : ./des/des_linked_list_data_mod.f \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_linked_list_data_mod.f  -o $(DPO)des_linked_list_data_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DES_LINKED_LIST_FUNCS.mod : ./des/des_linked_list_funcs_mod.f \
            $(DPO)DES_LINKED_LIST_DATA.mod \
            $(DPO)DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_linked_list_funcs_mod.f  -o $(DPO)des_linked_list_funcs_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DESMPI.mod : ./des/desmpi_mod.f \
            $(DPO)PARALLEL_MPI.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DESGRID.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DES_BC.mod \
            $(DPO)DESMPI_WRAPPER.mod \
            $(DPO)SENDRECVNODE.mod \
            $(DPO)MFIX_PIC.mod \
            des/desgrid_functions.inc                                    \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desmpi_mod.f  -o $(DPO)desmpi_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DESMPI_WRAPPER.mod : ./des/desmpi_wrapper_mod.f \
            $(DPO)PARALLEL_MPI.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desmpi_wrapper_mod.f  -o $(DPO)desmpi_wrapper_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DES_RXNS.mod : ./des/des_rxns_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)RXN_COM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_rxns_mod.f  -o $(DPO)des_rxns_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DES_STL_FUNCTIONS.mod : ./des/des_stl_functions_mod.f \
            $(DPO)STL.mod \
            $(DPO)PARAM.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)MPI_UTILITY.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_stl_functions_mod.f  -o $(DPO)des_stl_functions_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DES_THERMO.mod : ./des/des_thermo_mod.f \
            $(DPO)PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_mod.f  -o $(DPO)des_thermo_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DISCRETELEMENT.mod : ./des/discretelement_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/discretelement_mod.f  -o $(DPO)discretelement_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)INTERPOLATION.mod : ./des/interpolation_mod.f \
            $(DPO)CONSTANT.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PARAM1.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)PARAM.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/interpolation_mod.f  -o $(DPO)interpolation_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)MFIX_PIC.mod : ./des/mfix_pic_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mfix_pic_mod.f  -o $(DPO)mfix_pic_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)MPPIC_WALLBC.mod : ./des/mppic_wallbc_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RANDOMNO.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)RUN.mod \
            $(DPO)STL.mod \
            $(DPO)DES_STL_FUNCTIONS.mod \
            $(DPO)PARALLEL.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mppic_wallbc_mod.f  -o $(DPO)mppic_wallbc_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)RANDOMNO.mod : ./des/randomno_mod.f \
            $(DPO)CONSTANT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/randomno_mod.f  -o $(DPO)randomno_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)SENDRECVNODE.mod : ./des/sendrecvnode_mod.f \
            $(DPO)PARALLEL_MPI.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DESMPI_WRAPPER.mod \
            $(DPO)DESGRID.mod \
            function.inc                                                 \
            des/desgrid_functions.inc                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/sendrecvnode_mod.f  -o $(DPO)sendrecvnode_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)SOFTSPRING_FUNCS_CUTCELL.mod : ./des/softspring_funcs_cutcell_mod.f \
            $(DPO)RUN.mod \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)FUNITS.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)STL.mod \
            $(DPO)DES_STL_FUNCTIONS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/softspring_funcs_cutcell_mod.f  -o $(DPO)softspring_funcs_cutcell_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)COMPAR.mod : ./dmp_modules/compar_mod.f \
            $(DPO)MPI.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/compar_mod.f  -o $(DPO)compar_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DBG_UTIL.mod : ./dmp_modules/dbg_util_mod.f \
            $(DPO)COMPAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PARALLEL_MPI.mod \
            $(DPO)INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/dbg_util_mod.f  -o $(DPO)dbg_util_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)DEBUG.mod : ./dmp_modules/debug_mod.f \
            $(DPO)FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/debug_mod.f  -o $(DPO)debug_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)GRIDMAP.mod : ./dmp_modules/gridmap_mod.f \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PARALLEL_MPI.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)COMPAR.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)ERROR_MANAGER.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/gridmap_mod.f  -o $(DPO)gridmap_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)MPI.mod : ./dmp_modules/mpi_mod.f \
            mpif.h                                                      
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_mod.f  -o $(DPO)mpi_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)MPI_UTILITY.mod : ./dmp_modules/mpi_utility_mod.f \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PARALLEL_MPI.mod \
            $(DPO)DEBUG.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_utility_mod.f  -o $(DPO)mpi_utility_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)PARALLEL_MPI.mod : ./dmp_modules/parallel_mpi_mod.f \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/parallel_mpi_mod.f  -o $(DPO)parallel_mpi_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)SENDRECV3.mod : ./dmp_modules/sendrecv3_mod.f \
            $(DPO)PARALLEL_MPI.mod \
            $(DPO)DEBUG.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)MPI.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv3_mod.f  -o $(DPO)sendrecv3_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)SENDRECV.mod : ./dmp_modules/sendrecv_mod.f \
            $(DPO)PARALLEL_MPI.mod \
            $(DPO)DEBUG.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)MPI.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv_mod.f  -o $(DPO)sendrecv_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)GHDTHEORY.mod : ./GhdTheory/ghdtheory_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ghdtheory_mod.f  -o $(DPO)ghdtheory_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)QMOMK_BC.mod : ./qmomk/qmomk_bc_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)BC.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)QMOMK_QUADRATURE.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_bc_mod.f  -o $(DPO)qmomk_bc_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)QMOMK_COLLISION.mod : ./qmomk/qmomk_collision_mod.f \
            $(DPO)QMOMK_PARAMETERS.mod \
            $(DPO)QMOMK_QUADRATURE.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_collision_mod.f  -o $(DPO)qmomk_collision_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)QMOMK_FLUXES.mod : ./qmomk/qmomk_fluxes_mod.f \
            $(DPO)QMOMK_PARAMETERS.mod \
            $(DPO)QMOMK_COLLISION.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_fluxes_mod.f  -o $(DPO)qmomk_fluxes_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)QMOM_KINETIC_EQUATION.mod : ./qmomk/qmom_kinetic_equation_mod.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)QMOMK_PARAMETERS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmom_kinetic_equation_mod.f  -o $(DPO)qmom_kinetic_equation_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)QMOMK_PARAMETERS.mod : ./qmomk/qmomk_parameters_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_parameters_mod.f  -o $(DPO)qmomk_parameters_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)QMOMK_QUADRATURE.mod : ./qmomk/qmomk_quadrature_mod.f \
            $(DPO)QMOMK_TOOLS.mod \
            $(DPO)QMOMK_PARAMETERS.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_quadrature_mod.f  -o $(DPO)qmomk_quadrature_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)QMOMK_TOOLS.mod : ./qmomk/qmomk_tools_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_tools_mod.f  -o $(DPO)qmomk_tools_mod.$(OBJ_EXT) -module $(DPO)
$(DPO)accum_resid.$(OBJ_EXT) : accum_resid.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) accum_resid.f  -o $(DPO)accum_resid.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_u_g.$(OBJ_EXT) : adjust_a_u_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)USR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_u_g.f  -o $(DPO)adjust_a_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_u_s.$(OBJ_EXT) : adjust_a_u_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_u_s.f  -o $(DPO)adjust_a_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_v_g.$(OBJ_EXT) : adjust_a_v_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_v_g.f  -o $(DPO)adjust_a_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_v_s.$(OBJ_EXT) : adjust_a_v_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_v_s.f  -o $(DPO)adjust_a_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_w_g.$(OBJ_EXT) : adjust_a_w_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_w_g.f  -o $(DPO)adjust_a_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_a_w_s.$(OBJ_EXT) : adjust_a_w_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_a_w_s.f  -o $(DPO)adjust_a_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_dt.$(OBJ_EXT) : adjust_dt.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_dt.f  -o $(DPO)adjust_dt.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_eps.$(OBJ_EXT) : adjust_eps.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_eps.f  -o $(DPO)adjust_eps.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_leq.$(OBJ_EXT) : adjust_leq.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)LEQSOL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_leq.f  -o $(DPO)adjust_leq.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_rop.$(OBJ_EXT) : adjust_rop.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_rop.f  -o $(DPO)adjust_rop.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_theta.$(OBJ_EXT) : adjust_theta.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) adjust_theta.f  -o $(DPO)adjust_theta.$(OBJ_EXT) -module $(DPO)
$(DPO)allocate_arrays.$(OBJ_EXT) : allocate_arrays.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)AMBM.mod \
            $(DPO)CONT.mod \
            $(DPO)DRAG.mod \
            $(DPO)ENERGY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PSCOR.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)SCALARS.mod \
            $(DPO)TURB.mod \
            $(DPO)TAU_G.mod \
            $(DPO)TAU_S.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)TMP_ARRAY1.mod \
            $(DPO)TRACE.mod \
            $(DPO)VISC_G.mod \
            $(DPO)VISC_S.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)MFLUX.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)CDIST.mod \
            $(DPO)DES_RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) allocate_arrays.f  -o $(DPO)allocate_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)bc_phi.$(OBJ_EXT) : bc_phi.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)CUTCELL.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) bc_phi.f  -o $(DPO)bc_phi.$(OBJ_EXT) -module $(DPO)
$(DPO)bc_theta.$(OBJ_EXT) : bc_theta.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)TURB.mod \
            $(DPO)RXNS.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) bc_theta.f  -o $(DPO)bc_theta.$(OBJ_EXT) -module $(DPO)
$(DPO)b_m_p_star.$(OBJ_EXT) : b_m_p_star.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)RXNS.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) b_m_p_star.f  -o $(DPO)b_m_p_star.$(OBJ_EXT) -module $(DPO)
$(DPO)bound_x.$(OBJ_EXT) : bound_x.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) bound_x.f  -o $(DPO)bound_x.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_cell.$(OBJ_EXT) : calc_cell.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_cell.f  -o $(DPO)calc_cell.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_coeff.$(OBJ_EXT) : calc_coeff.f \
            $(DPO)PARAM1.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)RXNS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)VISC_G.mod \
            $(DPO)VISC_S.mod \
            $(DPO)TAU_G.mod \
            $(DPO)TAU_S.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_coeff.f  -o $(DPO)calc_coeff.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_d.$(OBJ_EXT) : calc_d.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)SCALES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)DISCRETELEMENT.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_d.f  -o $(DPO)calc_d.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_dif_g.$(OBJ_EXT) : calc_dif_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)SCALES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)RUN.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_dif_g.f  -o $(DPO)calc_dif_g.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_dif_s.$(OBJ_EXT) : calc_dif_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)RUN.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_dif_s.f  -o $(DPO)calc_dif_s.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_drag.$(OBJ_EXT) : calc_drag.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)DRAG.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_drag.f  -o $(DPO)calc_drag.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_e.$(OBJ_EXT) : calc_e.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_e.f  -o $(DPO)calc_e.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_gama.$(OBJ_EXT) : calc_gama.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)ENERGY.mod \
            $(DPO)RXNS.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_gama.f  -o $(DPO)calc_gama.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_grbdry.$(OBJ_EXT) : calc_grbdry.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)TURB.mod \
            $(DPO)VISC_S.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_grbdry.f  -o $(DPO)calc_grbdry.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_h.$(OBJ_EXT) : calc_h.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_h.f  -o $(DPO)calc_h.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_k_cp.$(OBJ_EXT) : calc_k_cp.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)PSCOR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)VISC_S.mod \
            $(DPO)TRACE.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_k_cp.f  -o $(DPO)calc_k_cp.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_k_g.$(OBJ_EXT) : calc_k_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)RUN.mod \
            $(DPO)SENDRECV.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_k_g.f  -o $(DPO)calc_k_g.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_k_s.$(OBJ_EXT) : calc_k_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)RUN.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_k_s.f  -o $(DPO)calc_k_s.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_mflux.$(OBJ_EXT) : calc_mflux.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)MFLUX.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_mflux.f  -o $(DPO)calc_mflux.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_mu_g.$(OBJ_EXT) : calc_mu_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)VISC_S.mod \
            $(DPO)INDICES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DRAG.mod \
            $(DPO)RUN.mod \
            $(DPO)TURB.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)MMS.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_mu_g.f  -o $(DPO)calc_mu_g.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_mu_s.$(OBJ_EXT) : calc_mu_s.f \
            $(DPO)RUN.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)MMS.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TRACE.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)TURB.mod \
            $(DPO)DRAG.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)VISC_G.mod \
            $(DPO)IS.mod \
            $(DPO)SENDRECV.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            s_pr1.inc                                                    \
            s_pr2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_mu_s.f  -o $(DPO)calc_mu_s.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_mw.$(OBJ_EXT) : calc_mw.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_mw.f  -o $(DPO)calc_mw.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_outflow.$(OBJ_EXT) : calc_outflow.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)BC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_outflow.f  -o $(DPO)calc_outflow.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_p_star.$(OBJ_EXT) : calc_p_star.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)COMPAR.mod \
            $(DPO)RUN.mod \
            $(DPO)VISC_S.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)TOLERANC.mod \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_p_star.f  -o $(DPO)calc_p_star.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_resid.$(OBJ_EXT) : calc_resid.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)RUN.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)BC.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)RXNS.mod \
            $(DPO)MFLUX.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_resid.f  -o $(DPO)calc_resid.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_s_ddot_s.$(OBJ_EXT) : calc_s_ddot_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_s_ddot_s.f  -o $(DPO)calc_s_ddot_s.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_trd_g.$(OBJ_EXT) : calc_trd_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)BC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_trd_g.f  -o $(DPO)calc_trd_g.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_trd_s.$(OBJ_EXT) : calc_trd_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)BC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_trd_s.f  -o $(DPO)calc_trd_s.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_u_friction.$(OBJ_EXT) : calc_u_friction.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)BC.mod \
            $(DPO)RUN.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_u_friction.f  -o $(DPO)calc_u_friction.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_vol_fr.$(OBJ_EXT) : calc_vol_fr.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)VISC_S.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)FLDVAR.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_vol_fr.f  -o $(DPO)calc_vol_fr.$(OBJ_EXT) -module $(DPO)
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_vol_fr.f  -o $(DPO)calc_vol_fr.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_xsi.$(OBJ_EXT) : calc_xsi.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)CHISCHEME.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            xsi1.inc                                                     \
            function.inc                                                 \
            xsi2.inc                                                    
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_xsi.f  -o $(DPO)calc_xsi.$(OBJ_EXT) -module $(DPO)
$(DPO)cal_d.$(OBJ_EXT) : cal_d.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_S.mod \
            $(DPO)BC.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) cal_d.f  -o $(DPO)cal_d.$(OBJ_EXT) -module $(DPO)
$(DPO)check_ab_m.$(OBJ_EXT) : check_ab_m.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) check_ab_m.f  -o $(DPO)check_ab_m.$(OBJ_EXT) -module $(DPO)
$(DPO)check_convergence.$(OBJ_EXT) : check_convergence.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SCALARS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_convergence.f  -o $(DPO)check_convergence.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_20.$(OBJ_EXT) : check_data_20.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)SCALARS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) check_data_20.f  -o $(DPO)check_data_20.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_30.$(OBJ_EXT) : check_data_30.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RXNS.mod \
            $(DPO)VISC_S.mod \
            $(DPO)VISC_G.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MMS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) check_data_30.f  -o $(DPO)check_data_30.$(OBJ_EXT) -module $(DPO)
$(DPO)check_mass_balance.$(OBJ_EXT) : check_mass_balance.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RXNS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)RUN.mod \
            $(DPO)BC.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)CHECK.mod \
            $(DPO)MFLUX.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) check_mass_balance.f  -o $(DPO)check_mass_balance.$(OBJ_EXT) -module $(DPO)
$(DPO)check_one_axis.$(OBJ_EXT) : check_one_axis.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_one_axis.f  -o $(DPO)check_one_axis.$(OBJ_EXT) -module $(DPO)
$(DPO)check_plane.$(OBJ_EXT) : check_plane.f \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_plane.f  -o $(DPO)check_plane.$(OBJ_EXT) -module $(DPO)
$(DPO)cn_extrapol.$(OBJ_EXT) : cn_extrapol.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)SCALARS.mod \
            $(DPO)TRACE.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) cn_extrapol.f  -o $(DPO)cn_extrapol.$(OBJ_EXT) -module $(DPO)
$(DPO)compare.$(OBJ_EXT) : compare.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) compare.f  -o $(DPO)compare.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_phi.$(OBJ_EXT) : conv_dif_phi.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)SENDRECV3.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)IS.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            function3.inc                                                \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_phi.f  -o $(DPO)conv_dif_phi.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_u_g.$(OBJ_EXT) : conv_dif_u_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)VISC_G.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)MFLUX.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)SENDRECV3.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_u_g.f  -o $(DPO)conv_dif_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_u_s.$(OBJ_EXT) : conv_dif_u_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)VISC_S.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)MFLUX.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)SENDRECV3.mod \
            $(DPO)VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_u_s.f  -o $(DPO)conv_dif_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_v_g.$(OBJ_EXT) : conv_dif_v_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)VISC_G.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)MFLUX.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)SENDRECV3.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_v_g.f  -o $(DPO)conv_dif_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_v_s.$(OBJ_EXT) : conv_dif_v_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)VISC_S.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)MFLUX.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)SENDRECV3.mod \
            $(DPO)VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_v_s.f  -o $(DPO)conv_dif_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_w_g.$(OBJ_EXT) : conv_dif_w_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)VISC_G.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)MFLUX.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)SENDRECV3.mod \
            $(DPO)VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_w_g.f  -o $(DPO)conv_dif_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_dif_w_s.$(OBJ_EXT) : conv_dif_w_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)VISC_S.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)MFLUX.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)SENDRECV3.mod \
            $(DPO)VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_dif_w_s.f  -o $(DPO)conv_dif_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_pp_g.$(OBJ_EXT) : conv_pp_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PGCOR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MFLUX.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_pp_g.f  -o $(DPO)conv_pp_g.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_rop.$(OBJ_EXT) : conv_rop.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)MFLUX.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)XSI_ARRAY.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_rop.f  -o $(DPO)conv_rop.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_rop_g.$(OBJ_EXT) : conv_rop_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PGCOR.mod \
            $(DPO)XSI_ARRAY.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_rop_g.f  -o $(DPO)conv_rop_g.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_rop_s.$(OBJ_EXT) : conv_rop_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)XSI_ARRAY.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_rop_s.f  -o $(DPO)conv_rop_s.$(OBJ_EXT) -module $(DPO)
$(DPO)conv_source_epp.$(OBJ_EXT) : conv_source_epp.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RXNS.mod \
            $(DPO)INDICES.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)PS.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) conv_source_epp.f  -o $(DPO)conv_source_epp.$(OBJ_EXT) -module $(DPO)
$(DPO)copy_a.$(OBJ_EXT) : copy_a.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PHYSPROP.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) copy_a.f  -o $(DPO)copy_a.$(OBJ_EXT) -module $(DPO)
$(DPO)corner.$(OBJ_EXT) : corner.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)MATRIX.mod \
            $(DPO)CORNER.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) corner.f  -o $(DPO)corner.$(OBJ_EXT) -module $(DPO)
$(DPO)correct_0.$(OBJ_EXT) : correct_0.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PGCOR.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)COMPAR.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) correct_0.f  -o $(DPO)correct_0.$(OBJ_EXT) -module $(DPO)
$(DPO)correct_1.$(OBJ_EXT) : correct_1.f \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PSCOR.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)INDICES.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)VISC_S.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) correct_1.f  -o $(DPO)correct_1.$(OBJ_EXT) -module $(DPO)
$(DPO)dgtsl.$(OBJ_EXT) : dgtsl.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) dgtsl.f  -o $(DPO)dgtsl.$(OBJ_EXT) -module $(DPO)
$(DPO)dif_u_is.$(OBJ_EXT) : dif_u_is.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) dif_u_is.f  -o $(DPO)dif_u_is.$(OBJ_EXT) -module $(DPO)
$(DPO)dif_v_is.$(OBJ_EXT) : dif_v_is.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) dif_v_is.f  -o $(DPO)dif_v_is.$(OBJ_EXT) -module $(DPO)
$(DPO)dif_w_is.$(OBJ_EXT) : dif_w_is.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) dif_w_is.f  -o $(DPO)dif_w_is.$(OBJ_EXT) -module $(DPO)
$(DPO)discretize.$(OBJ_EXT) : discretize.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) discretize.f  -o $(DPO)discretize.$(OBJ_EXT) -module $(DPO)
$(DPO)display_resid.$(OBJ_EXT) : display_resid.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)SCALARS.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) display_resid.f  -o $(DPO)display_resid.$(OBJ_EXT) -module $(DPO)
$(DPO)drag_gs.$(OBJ_EXT) : drag_gs.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DRAG.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)MMS.mod \
            $(DPO)CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) drag_gs.f  -o $(DPO)drag_gs.$(OBJ_EXT) -module $(DPO)
$(DPO)drag_ss.$(OBJ_EXT) : drag_ss.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DRAG.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) drag_ss.f  -o $(DPO)drag_ss.$(OBJ_EXT) -module $(DPO)
$(DPO)eosg.$(OBJ_EXT) : eosg.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)SCALES.mod \
            sc_p_g1.inc                                                  \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) eosg.f  -o $(DPO)eosg.$(OBJ_EXT) -module $(DPO)
$(DPO)eoss.$(OBJ_EXT) : eoss.f \
            $(DPO)PHYSPROP.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FUNITS.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) eoss.f  -o $(DPO)eoss.$(OBJ_EXT) -module $(DPO)
$(DPO)equal.$(OBJ_EXT) : equal.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) equal.f  -o $(DPO)equal.$(OBJ_EXT) -module $(DPO)
$(DPO)error_routine.$(OBJ_EXT) : error_routine.f \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) error_routine.f  -o $(DPO)error_routine.$(OBJ_EXT) -module $(DPO)
$(DPO)exchange.$(OBJ_EXT) : exchange.f \
            $(DPO)COEFF.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) exchange.f  -o $(DPO)exchange.$(OBJ_EXT) -module $(DPO)
$(DPO)exit.$(OBJ_EXT) : exit.f \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) exit.f  -o $(DPO)exit.$(OBJ_EXT) -module $(DPO)
$(DPO)flow_to_vel.$(OBJ_EXT) : flow_to_vel.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)BC.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)SCALES.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MFIX_PIC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) flow_to_vel.f  -o $(DPO)flow_to_vel.$(OBJ_EXT) -module $(DPO)
$(DPO)g_0.$(OBJ_EXT) : g_0.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) g_0.f  -o $(DPO)g_0.$(OBJ_EXT) -module $(DPO)
$(DPO)get_bc_area.$(OBJ_EXT) : get_bc_area.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)INDICES.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)CUTCELL.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) get_bc_area.f  -o $(DPO)get_bc_area.$(OBJ_EXT) -module $(DPO)
$(DPO)get_data.$(OBJ_EXT) : get_data.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)GRIDMAP.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)DASHBOARD.mod \
            $(DPO)DES_STL_FUNCTIONS.mod \
            $(DPO)VISC_G.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_data.f  -o $(DPO)get_data.$(OBJ_EXT) -module $(DPO)
$(DPO)get_eq.$(OBJ_EXT) : get_eq.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_eq.f  -o $(DPO)get_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)get_flow_bc.$(OBJ_EXT) : get_flow_bc.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_flow_bc.f  -o $(DPO)get_flow_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)get_hloss.$(OBJ_EXT) : get_hloss.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)BC.mod \
            $(DPO)INDICES.mod \
            $(DPO)ENERGY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_hloss.f  -o $(DPO)get_hloss.$(OBJ_EXT) -module $(DPO)
$(DPO)get_is.$(OBJ_EXT) : get_is.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)IS.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_is.f  -o $(DPO)get_is.$(OBJ_EXT) -module $(DPO)
$(DPO)get_philoss.$(OBJ_EXT) : get_philoss.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)BC.mod \
            $(DPO)INDICES.mod \
            $(DPO)ENERGY.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) get_philoss.f  -o $(DPO)get_philoss.$(OBJ_EXT) -module $(DPO)
$(DPO)get_ps.$(OBJ_EXT) : get_ps.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PS.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_ps.f  -o $(DPO)get_ps.$(OBJ_EXT) -module $(DPO)
$(DPO)get_smass.$(OBJ_EXT) : get_smass.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) get_smass.f  -o $(DPO)get_smass.$(OBJ_EXT) -module $(DPO)
$(DPO)get_stats.$(OBJ_EXT) : get_stats.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) get_stats.f  -o $(DPO)get_stats.$(OBJ_EXT) -module $(DPO)
$(DPO)get_walls_bc.$(OBJ_EXT) : get_walls_bc.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) get_walls_bc.f  -o $(DPO)get_walls_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)in_bin_512.$(OBJ_EXT) : in_bin_512.f \
            $(DPO)MACHINE.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) in_bin_512.f  -o $(DPO)in_bin_512.$(OBJ_EXT) -module $(DPO)
$(DPO)in_bin_512i.$(OBJ_EXT) : in_bin_512i.f \
            $(DPO)MACHINE.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) in_bin_512i.f  -o $(DPO)in_bin_512i.$(OBJ_EXT) -module $(DPO)
$(DPO)init_ab_m.$(OBJ_EXT) : init_ab_m.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) init_ab_m.f  -o $(DPO)init_ab_m.$(OBJ_EXT) -module $(DPO)
$(DPO)init_fvars.$(OBJ_EXT) : init_fvars.f \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RXNS.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) init_fvars.f  -o $(DPO)init_fvars.$(OBJ_EXT) -module $(DPO)
$(DPO)init_namelist.$(OBJ_EXT) : init_namelist.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)IC.mod \
            $(DPO)BC.mod \
            $(DPO)PS.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)SCALES.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)RXNS.mod \
            $(DPO)SCALARS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CDIST.mod \
            $(DPO)STIFF_CHEM.mod \
            namelist.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) init_namelist.f  -o $(DPO)init_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)init_resid.$(OBJ_EXT) : init_resid.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RESIDUAL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) init_resid.f  -o $(DPO)init_resid.$(OBJ_EXT) -module $(DPO)
$(DPO)iterate.$(OBJ_EXT) : iterate.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)TIME_CPU.mod \
            $(DPO)PSCOR.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)VISC_G.mod \
            $(DPO)PGCOR.mod \
            $(DPO)CONT.mod \
            $(DPO)SCALARS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)VTK.mod \
            $(DPO)DASHBOARD.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)STIFF_CHEM.mod \
            $(DPO)RXNS.mod \
            $(DPO)MMS.mod \
            $(DPO)BC.mod \
            $(DPO)CONSTANT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) iterate.f  -o $(DPO)iterate.$(OBJ_EXT) -module $(DPO)
$(DPO)k_epsilon_prop.$(OBJ_EXT) : k_epsilon_prop.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DRAG.mod \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)VISC_S.mod \
            $(DPO)TRACE.mod \
            $(DPO)INDICES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)TURB.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TAU_G.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) k_epsilon_prop.f  -o $(DPO)k_epsilon_prop.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory_drag_ss.$(OBJ_EXT) : kintheory_drag_ss.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DRAG.mod \
            $(DPO)KINTHEORY.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_drag_ss.f  -o $(DPO)kintheory_drag_ss.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory_energy_dissipation_ss.$(OBJ_EXT) : kintheory_energy_dissipation_ss.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)KINTHEORY.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_energy_dissipation_ss.f  -o $(DPO)kintheory_energy_dissipation_ss.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory_u_s.$(OBJ_EXT) : kintheory_u_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_u_s.f  -o $(DPO)kintheory_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory_v_s.$(OBJ_EXT) : kintheory_v_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_v_s.f  -o $(DPO)kintheory_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)kintheory_w_s.$(OBJ_EXT) : kintheory_w_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_w_s.f  -o $(DPO)kintheory_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)leq_bicgs.$(OBJ_EXT) : leq_bicgs.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)FUNITS.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)CUTCELL.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) leq_bicgs.f  -o $(DPO)leq_bicgs.$(OBJ_EXT) -module $(DPO)
$(DPO)leq_bicgst.$(OBJ_EXT) : leq_bicgst.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)FUNITS.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) leq_bicgst.f  -o $(DPO)leq_bicgst.$(OBJ_EXT) -module $(DPO)
$(DPO)leq_cg.$(OBJ_EXT) : leq_cg.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)FUNITS.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) leq_cg.f  -o $(DPO)leq_cg.$(OBJ_EXT) -module $(DPO)
$(DPO)leq_gmres.$(OBJ_EXT) : leq_gmres.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)DEBUG.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GRIDMAP.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) leq_gmres.f  -o $(DPO)leq_gmres.$(OBJ_EXT) -module $(DPO)
$(DPO)leq_sor.$(OBJ_EXT) : leq_sor.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)LEQSOL.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) leq_sor.f  -o $(DPO)leq_sor.$(OBJ_EXT) -module $(DPO)
$(DPO)line_too_big.$(OBJ_EXT) : line_too_big.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) line_too_big.f  -o $(DPO)line_too_big.$(OBJ_EXT) -module $(DPO)
$(DPO)location_check.$(OBJ_EXT) : location_check.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) location_check.f  -o $(DPO)location_check.$(OBJ_EXT) -module $(DPO)
$(DPO)location.$(OBJ_EXT) : location.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) location.f  -o $(DPO)location.$(OBJ_EXT) -module $(DPO)
$(DPO)machine.$(OBJ_EXT) : machine.f \
            $(DPO)MACHINE.mod \
            $(DPO)PARAM.mod \
            $(DPO)RUN.mod \
            $(DPO)FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) machine.f  -o $(DPO)machine.$(OBJ_EXT) -module $(DPO)
$(DPO)make_upper_case.$(OBJ_EXT) : make_upper_case.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) make_upper_case.f  -o $(DPO)make_upper_case.$(OBJ_EXT) -module $(DPO)
$(DPO)mark_phase_4_cor.$(OBJ_EXT) : mark_phase_4_cor.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)VISC_S.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) mark_phase_4_cor.f  -o $(DPO)mark_phase_4_cor.$(OBJ_EXT) -module $(DPO)
	$(FORTRAN_CMD) $(FORT_FLAGS) mark_phase_4_cor.f  -o $(DPO)mark_phase_4_cor.$(OBJ_EXT) -module $(DPO)
$(DPO)mfix.$(OBJ_EXT) : mfix.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)TIME_CPU.mod \
            $(DPO)FUNITS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PARALLEL_MPI.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)CDIST.mod \
            $(DPO)MFIX_NETCDF.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)DASHBOARD.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)SENDRECV3.mod \
            $(DPO)INDICES.mod \
            $(DPO)LEQSOL.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) mfix.f  -o $(DPO)mfix.$(OBJ_EXT) -module $(DPO)
$(DPO)mod_bc_i.$(OBJ_EXT) : mod_bc_i.f \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)ERROR_MANAGER.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) mod_bc_i.f  -o $(DPO)mod_bc_i.$(OBJ_EXT) -module $(DPO)
$(DPO)mod_bc_j.$(OBJ_EXT) : mod_bc_j.f \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)ERROR_MANAGER.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) mod_bc_j.f  -o $(DPO)mod_bc_j.$(OBJ_EXT) -module $(DPO)
$(DPO)mod_bc_k.$(OBJ_EXT) : mod_bc_k.f \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)ERROR_MANAGER.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) mod_bc_k.f  -o $(DPO)mod_bc_k.$(OBJ_EXT) -module $(DPO)
$(DPO)open_file.$(OBJ_EXT) : open_file.f \
            $(DPO)CDIST.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) open_file.f  -o $(DPO)open_file.$(OBJ_EXT) -module $(DPO)
$(DPO)open_files.$(OBJ_EXT) : open_files.f \
            $(DPO)MACHINE.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)CDIST.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) open_files.f  -o $(DPO)open_files.$(OBJ_EXT) -module $(DPO)
$(DPO)out_array_c.$(OBJ_EXT) : out_array_c.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) out_array_c.f  -o $(DPO)out_array_c.$(OBJ_EXT) -module $(DPO)
$(DPO)out_array.$(OBJ_EXT) : out_array.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) out_array.f  -o $(DPO)out_array.$(OBJ_EXT) -module $(DPO)
$(DPO)out_array_kc.$(OBJ_EXT) : out_array_kc.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) out_array_kc.f  -o $(DPO)out_array_kc.$(OBJ_EXT) -module $(DPO)
$(DPO)out_array_k.$(OBJ_EXT) : out_array_k.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) out_array_k.f  -o $(DPO)out_array_k.$(OBJ_EXT) -module $(DPO)
$(DPO)out_bin_512.$(OBJ_EXT) : out_bin_512.f \
            $(DPO)MACHINE.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) out_bin_512.f  -o $(DPO)out_bin_512.$(OBJ_EXT) -module $(DPO)
$(DPO)out_bin_512i.$(OBJ_EXT) : out_bin_512i.f \
            $(DPO)MACHINE.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) out_bin_512i.f  -o $(DPO)out_bin_512i.$(OBJ_EXT) -module $(DPO)
$(DPO)out_bin_512r.$(OBJ_EXT) : out_bin_512r.f \
            $(DPO)MACHINE.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) out_bin_512r.f  -o $(DPO)out_bin_512r.$(OBJ_EXT) -module $(DPO)
$(DPO)out_bin_r.$(OBJ_EXT) : out_bin_r.f \
            $(DPO)PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) out_bin_r.f  -o $(DPO)out_bin_r.$(OBJ_EXT) -module $(DPO)
$(DPO)parse_line.$(OBJ_EXT) : parse_line.f \
            $(DPO)COMPAR.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARSE.mod \
            $(DPO)RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_line.f  -o $(DPO)parse_line.$(OBJ_EXT) -module $(DPO)
$(DPO)parse_resid_string.$(OBJ_EXT) : parse_resid_string.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_resid_string.f  -o $(DPO)parse_resid_string.$(OBJ_EXT) -module $(DPO)
$(DPO)parse_rxn.$(OBJ_EXT) : parse_rxn.f \
            $(DPO)COMPAR.mod \
            $(DPO)FUNITS.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARSE.mod \
            $(DPO)RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_rxn.f  -o $(DPO)parse_rxn.$(OBJ_EXT) -module $(DPO)
$(DPO)partial_elim.$(OBJ_EXT) : partial_elim.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)MATRIX.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DRAG.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) partial_elim.f  -o $(DPO)partial_elim.$(OBJ_EXT) -module $(DPO)
$(DPO)physical_prop.$(OBJ_EXT) : physical_prop.f \
            $(DPO)COMPAR.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)COEFF.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)SCALARS.mod \
            $(DPO)CUTCELL.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) physical_prop.f  -o $(DPO)physical_prop.$(OBJ_EXT) -module $(DPO)
$(DPO)read_database.$(OBJ_EXT) : read_database.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)RXNS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)ERROR_MANAGER.mod \
            mfix_directory_path.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) read_database.f  -o $(DPO)read_database.$(OBJ_EXT) -module $(DPO)
$(DPO)read_namelist.$(OBJ_EXT) : read_namelist.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)IC.mod \
            $(DPO)IS.mod \
            $(DPO)BC.mod \
            $(DPO)PS.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)INDICES.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)FUNITS.mod \
            $(DPO)SCALES.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)RXNS.mod \
            $(DPO)SCALARS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)USR.mod \
            $(DPO)DES_BC.mod \
            $(DPO)DES_IC.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)STIFF_CHEM.mod \
            $(DPO)CDIST.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)VTK.mod \
            $(DPO)POLYGON.mod \
            $(DPO)DASHBOARD.mod \
            $(DPO)STL.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            usrnlst.inc                                                  \
            namelist.inc                                                 \
            des/desnamelist.inc                                          \
            cartesian_grid/cartesian_grid_namelist.inc                   \
            qmomk/qmomknamelist.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) read_namelist.f  -o $(DPO)read_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)read_res0.$(OBJ_EXT) : read_res0.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)IC.mod \
            $(DPO)BC.mod \
            $(DPO)IS.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)SCALES.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)SCALARS.mod \
            $(DPO)RXNS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)STIFF_CHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) read_res0.f  -o $(DPO)read_res0.$(OBJ_EXT) -module $(DPO)
$(DPO)read_res1.$(OBJ_EXT) : read_res1.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)RXNS.mod \
            $(DPO)SCALARS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)ENERGY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)CDIST.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)MFIX_NETCDF.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) read_res1.f  -o $(DPO)read_res1.$(OBJ_EXT) -module $(DPO)
$(DPO)remove_comment.$(OBJ_EXT) : remove_comment.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) remove_comment.f  -o $(DPO)remove_comment.$(OBJ_EXT) -module $(DPO)
$(DPO)reset_new.$(OBJ_EXT) : reset_new.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)TRACE.mod \
            $(DPO)RUN.mod \
            $(DPO)SCALARS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) reset_new.f  -o $(DPO)reset_new.$(OBJ_EXT) -module $(DPO)
$(DPO)rrates0.$(OBJ_EXT) : rrates0.f \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RXNS.mod \
            $(DPO)ENERGY.mod \
            $(DPO)PARAM.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)TOLERANC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) rrates0.f  -o $(DPO)rrates0.$(OBJ_EXT) -module $(DPO)
$(DPO)rrates.$(OBJ_EXT) : rrates.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RXNS.mod \
            $(DPO)ENERGY.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) rrates.f  -o $(DPO)rrates.$(OBJ_EXT) -module $(DPO)
$(DPO)rrates_init.$(OBJ_EXT) : rrates_init.f \
            $(DPO)ENERGY.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RXNS.mod \
            $(DPO)STIFF_CHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) rrates_init.f  -o $(DPO)rrates_init.$(OBJ_EXT) -module $(DPO)
$(DPO)scalar_prop.$(OBJ_EXT) : scalar_prop.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)SCALARS.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) scalar_prop.f  -o $(DPO)scalar_prop.$(OBJ_EXT) -module $(DPO)
$(DPO)seek_comment.$(OBJ_EXT) : seek_comment.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) seek_comment.f  -o $(DPO)seek_comment.$(OBJ_EXT) -module $(DPO)
$(DPO)seek_end.$(OBJ_EXT) : seek_end.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) seek_end.f  -o $(DPO)seek_end.$(OBJ_EXT) -module $(DPO)
$(DPO)set_bc0.$(OBJ_EXT) : set_bc0.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)BC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)FUNITS.mod \
            $(DPO)SCALES.mod \
            $(DPO)SCALARS.mod \
            $(DPO)BOUNDFUNIJK.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)MMS.mod \
            sc_p_g1.inc                                                  \
            function.inc                                                 \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_bc0.f  -o $(DPO)set_bc0.$(OBJ_EXT) -module $(DPO)
$(DPO)set_bc1.$(OBJ_EXT) : set_bc1.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)BC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_bc1.f  -o $(DPO)set_bc1.$(OBJ_EXT) -module $(DPO)
$(DPO)set_bc_flow.$(OBJ_EXT) : set_bc_flow.f \
            $(DPO)PHYSPROP.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            $(DPO)BC.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARAM.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)SCALARS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_bc_flow.f  -o $(DPO)set_bc_flow.$(OBJ_EXT) -module $(DPO)
$(DPO)set_constants.$(OBJ_EXT) : set_constants.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)ENERGY.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)FUNITS.mod \
            $(DPO)DRAG.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_constants.f  -o $(DPO)set_constants.$(OBJ_EXT) -module $(DPO)
$(DPO)set_constprop.$(OBJ_EXT) : set_constprop.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)VISC_G.mod \
            $(DPO)ENERGY.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)DRAG.mod \
            $(DPO)COMPAR.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)MMS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_constprop.f  -o $(DPO)set_constprop.$(OBJ_EXT) -module $(DPO)
$(DPO)set_flags.$(OBJ_EXT) : set_flags.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)BC.mod \
            $(DPO)IS.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)SENDRECV3.mod \
            $(DPO)BOUNDFUNIJK.mod \
            $(DPO)MPI_UTILITY.mod \
            function.inc                                                 \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) set_flags.f  -o $(DPO)set_flags.$(OBJ_EXT) -module $(DPO)
$(DPO)set_fluidbed_p.$(OBJ_EXT) : set_fluidbed_p.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)BC.mod \
            $(DPO)IC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)SCALES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            sc_p_g1.inc                                                  \
            b_force1.inc                                                 \
            function.inc                                                 \
            b_force2.inc                                                 \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_fluidbed_p.f  -o $(DPO)set_fluidbed_p.$(OBJ_EXT) -module $(DPO)
$(DPO)set_geometry1.$(OBJ_EXT) : set_geometry1.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_geometry1.f  -o $(DPO)set_geometry1.$(OBJ_EXT) -module $(DPO)
$(DPO)set_geometry.$(OBJ_EXT) : set_geometry.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)BC.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_geometry.f  -o $(DPO)set_geometry.$(OBJ_EXT) -module $(DPO)
$(DPO)set_icbc_flags.$(OBJ_EXT) : set_icbc_flags.f \
            $(DPO)RUN.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)IC.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)BC.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_icbc_flags.f  -o $(DPO)set_icbc_flags.$(OBJ_EXT) -module $(DPO)
$(DPO)set_ic.$(OBJ_EXT) : set_ic.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)IC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)INDICES.mod \
            $(DPO)SCALES.mod \
            $(DPO)ENERGY.mod \
            $(DPO)SCALARS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)RUN.mod \
            $(DPO)SENDRECV.mod \
            sc_p_g1.inc                                                  \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_ic.f  -o $(DPO)set_ic.$(OBJ_EXT) -module $(DPO)
$(DPO)set_increments3.$(OBJ_EXT) : set_increments3.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)INDICES.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)FUNITS.mod \
            function.inc                                                 \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) set_increments3.f  -o $(DPO)set_increments3.$(OBJ_EXT) -module $(DPO)
$(DPO)set_increments.$(OBJ_EXT) : set_increments.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)INDICES.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)FUNITS.mod \
            $(DPO)SCALARS.mod \
            $(DPO)RUN.mod \
            $(DPO)VISC_G.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)STL.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)BC.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)CDIST.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_increments.f  -o $(DPO)set_increments.$(OBJ_EXT) -module $(DPO)
$(DPO)set_index1a3.$(OBJ_EXT) : set_index1a3.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)BOUNDFUNIJK3.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_index1a3.f  -o $(DPO)set_index1a3.$(OBJ_EXT) -module $(DPO)
$(DPO)set_index1a.$(OBJ_EXT) : set_index1a.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)BOUNDFUNIJK.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_index1a.f  -o $(DPO)set_index1a.$(OBJ_EXT) -module $(DPO)
$(DPO)set_index1.$(OBJ_EXT) : set_index1.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_index1.f  -o $(DPO)set_index1.$(OBJ_EXT) -module $(DPO)
$(DPO)set_l_scale.$(OBJ_EXT) : set_l_scale.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)VISC_G.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_l_scale.f  -o $(DPO)set_l_scale.$(OBJ_EXT) -module $(DPO)
$(DPO)set_max2.$(OBJ_EXT) : set_max2.f \
            $(DPO)COMPAR.mod \
            $(DPO)GEOMETRY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) set_max2.f  -o $(DPO)set_max2.$(OBJ_EXT) -module $(DPO)
$(DPO)set_mw_mix_g.$(OBJ_EXT) : set_mw_mix_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_mw_mix_g.f  -o $(DPO)set_mw_mix_g.$(OBJ_EXT) -module $(DPO)
$(DPO)set_outflow.$(OBJ_EXT) : set_outflow.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)BC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)SCALARS.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MFLUX.mod \
            $(DPO)DISCRETELEMENT.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) set_outflow.f  -o $(DPO)set_outflow.$(OBJ_EXT) -module $(DPO)
$(DPO)set_ps.$(OBJ_EXT) : set_ps.f \
            $(DPO)PARAM.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)BC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)IC.mod \
            $(DPO)INDICES.mod \
            $(DPO)MFLUX.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PARAM1.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)USR.mod \
            $(DPO)RXNS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_ps.f  -o $(DPO)set_ps.$(OBJ_EXT) -module $(DPO)
$(DPO)set_ro_g.$(OBJ_EXT) : set_ro_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_ro_g.f  -o $(DPO)set_ro_g.$(OBJ_EXT) -module $(DPO)
$(DPO)set_ro_s.$(OBJ_EXT) : set_ro_s.f \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_ro_s.f  -o $(DPO)set_ro_s.$(OBJ_EXT) -module $(DPO)
$(DPO)set_wall_bc.$(OBJ_EXT) : set_wall_bc.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)BC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) set_wall_bc.f  -o $(DPO)set_wall_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)shift_dxyz.$(OBJ_EXT) : shift_dxyz.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) shift_dxyz.f  -o $(DPO)shift_dxyz.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_continuity.$(OBJ_EXT) : solve_continuity.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)CONT.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)AMBM.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)RUN.mod \
            $(DPO)PS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_continuity.f  -o $(DPO)solve_continuity.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_energy_eq.$(OBJ_EXT) : solve_energy_eq.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)DRAG.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)BC.mod \
            $(DPO)ENERGY.mod \
            $(DPO)RXNS.mod \
            $(DPO)AMBM.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)TMP_ARRAY1.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)MFLUX.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)PS.mod \
            $(DPO)MMS.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_energy_eq.f  -o $(DPO)solve_energy_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_epp.$(OBJ_EXT) : solve_epp.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PSCOR.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)AMBM.mod \
            $(DPO)TMP_ARRAY1.mod \
            $(DPO)PS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_epp.f  -o $(DPO)solve_epp.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_granular_energy.$(OBJ_EXT) : solve_granular_energy.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)VISC_S.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)DRAG.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)BC.mod \
            $(DPO)ENERGY.mod \
            $(DPO)RXNS.mod \
            $(DPO)AMBM.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MFLUX.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)MMS.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_granular_energy.f  -o $(DPO)solve_granular_energy.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_k_epsilon_eq.$(OBJ_EXT) : solve_k_epsilon_eq.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)DRAG.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)BC.mod \
            $(DPO)ENERGY.mod \
            $(DPO)RXNS.mod \
            $(DPO)TURB.mod \
            $(DPO)USR.mod \
            $(DPO)AMBM.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MFLUX.mod \
            $(DPO)CUTCELL.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_k_epsilon_eq.f  -o $(DPO)solve_k_epsilon_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_lin_eq.$(OBJ_EXT) : solve_lin_eq.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)LEQSOL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_lin_eq.f  -o $(DPO)solve_lin_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_pp_g.$(OBJ_EXT) : solve_pp_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PGCOR.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)RUN.mod \
            $(DPO)AMBM.mod \
            $(DPO)TMP_ARRAY1.mod \
            $(DPO)PS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_pp_g.f  -o $(DPO)solve_pp_g.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_scalar_eq.$(OBJ_EXT) : solve_scalar_eq.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)DRAG.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)BC.mod \
            $(DPO)ENERGY.mod \
            $(DPO)RXNS.mod \
            $(DPO)SCALARS.mod \
            $(DPO)AMBM.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MFLUX.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_scalar_eq.f  -o $(DPO)solve_scalar_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_species_eq.$(OBJ_EXT) : solve_species_eq.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)DRAG.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)BC.mod \
            $(DPO)ENERGY.mod \
            $(DPO)RXNS.mod \
            $(DPO)AMBM.mod \
            $(DPO)MATRIX.mod \
            $(DPO)CHISCHEME.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)MFLUX.mod \
            $(DPO)PS.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_species_eq.f  -o $(DPO)solve_species_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)solve_vel_star.$(OBJ_EXT) : solve_vel_star.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)INDICES.mod \
            $(DPO)DRAG.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)AMBM.mod \
            $(DPO)TMP_ARRAY1.mod \
            $(DPO)TMP_ARRAY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)PS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) solve_vel_star.f  -o $(DPO)solve_vel_star.$(OBJ_EXT) -module $(DPO)
$(DPO)source_granular_energy.$(OBJ_EXT) : source_granular_energy.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)DRAG.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)VISC_S.mod \
            $(DPO)TRACE.mod \
            $(DPO)TURB.mod \
            $(DPO)INDICES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)MMS.mod \
            $(DPO)RESIDUAL.mod \
            s_pr1.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            s_pr2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) source_granular_energy.f  -o $(DPO)source_granular_energy.$(OBJ_EXT) -module $(DPO)
$(DPO)source_phi.$(OBJ_EXT) : source_phi.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_S.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PS.mod \
            $(DPO)BC.mod \
            $(DPO)USR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) source_phi.f  -o $(DPO)source_phi.$(OBJ_EXT) -module $(DPO)
$(DPO)source_pp_g.$(OBJ_EXT) : source_pp_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PGCOR.mod \
            $(DPO)BC.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)XSI_ARRAY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_pp_g.f  -o $(DPO)source_pp_g.$(OBJ_EXT) -module $(DPO)
$(DPO)source_rop_g.$(OBJ_EXT) : source_rop_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PGCOR.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_rop_g.f  -o $(DPO)source_rop_g.$(OBJ_EXT) -module $(DPO)
$(DPO)source_rop_s.$(OBJ_EXT) : source_rop_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_rop_s.f  -o $(DPO)source_rop_s.$(OBJ_EXT) -module $(DPO)
$(DPO)source_u_g.$(OBJ_EXT) : source_u_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_G.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)MMS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)TURB.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PS.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_u_g.f  -o $(DPO)source_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)source_u_s.$(OBJ_EXT) : source_u_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_S.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)MMS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PS.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_u_s.f  -o $(DPO)source_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)source_v_g.$(OBJ_EXT) : source_v_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_G.mod \
            $(DPO)BC.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)MMS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PS.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_v_g.f  -o $(DPO)source_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)source_v_s.$(OBJ_EXT) : source_v_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_S.mod \
            $(DPO)BC.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)MMS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PS.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_v_s.f  -o $(DPO)source_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)source_w_g.$(OBJ_EXT) : source_w_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_G.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)MMS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PS.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_w_g.f  -o $(DPO)source_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)source_w_s.$(OBJ_EXT) : source_w_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_S.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)MMS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PS.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) source_w_s.f  -o $(DPO)source_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_u_g.$(OBJ_EXT) : tau_u_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)BC.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_u_g.f  -o $(DPO)tau_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_u_s.$(OBJ_EXT) : tau_u_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)COMPAR.mod \
            $(DPO)BC.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_u_s.f  -o $(DPO)tau_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_v_g.$(OBJ_EXT) : tau_v_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)COMPAR.mod \
            $(DPO)BC.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_v_g.f  -o $(DPO)tau_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_v_s.$(OBJ_EXT) : tau_v_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)COMPAR.mod \
            $(DPO)BC.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_v_s.f  -o $(DPO)tau_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_w_g.$(OBJ_EXT) : tau_w_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)COMPAR.mod \
            $(DPO)BC.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_w_g.f  -o $(DPO)tau_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)tau_w_s.$(OBJ_EXT) : tau_w_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)COMPAR.mod \
            $(DPO)BC.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_w_s.f  -o $(DPO)tau_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)test_lin_eq.$(OBJ_EXT) : test_lin_eq.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) test_lin_eq.f  -o $(DPO)test_lin_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)time_march.$(OBJ_EXT) : time_march.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)CONT.mod \
            $(DPO)TAU_G.mod \
            $(DPO)TAU_S.mod \
            $(DPO)VISC_G.mod \
            $(DPO)VISC_S.mod \
            $(DPO)FUNITS.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)SCALARS.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)DRAG.mod \
            $(DPO)RXNS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TIME_CPU.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)CDIST.mod \
            $(DPO)MFIX_NETCDF.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)VTK.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)DASHBOARD.mod \
            $(DPO)INDICES.mod \
            $(DPO)BC.mod \
            $(DPO)COEFF.mod \
            $(DPO)STIFF_CHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) time_march.f  -o $(DPO)time_march.$(OBJ_EXT) -module $(DPO)
$(DPO)transfer.$(OBJ_EXT) : transfer.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) transfer.f  -o $(DPO)transfer.$(OBJ_EXT) -module $(DPO)
$(DPO)transport_prop.$(OBJ_EXT) : transport_prop.f \
            $(DPO)PHYSPROP.mod \
            $(DPO)COEFF.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) transport_prop.f  -o $(DPO)transport_prop.$(OBJ_EXT) -module $(DPO)
$(DPO)undef_2_0.$(OBJ_EXT) : undef_2_0.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) undef_2_0.f  -o $(DPO)undef_2_0.$(OBJ_EXT) -module $(DPO)
$(DPO)under_relax.$(OBJ_EXT) : under_relax.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) under_relax.f  -o $(DPO)under_relax.$(OBJ_EXT) -module $(DPO)
$(DPO)update_old.$(OBJ_EXT) : update_old.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)TRACE.mod \
            $(DPO)VISC_S.mod \
            $(DPO)SCALARS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) update_old.f  -o $(DPO)update_old.$(OBJ_EXT) -module $(DPO)
$(DPO)usr0.$(OBJ_EXT) : usr0.f \
            $(DPO)USR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr0.f  -o $(DPO)usr0.$(OBJ_EXT) -module $(DPO)
$(DPO)usr1.$(OBJ_EXT) : usr1.f \
            $(DPO)USR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr1.f  -o $(DPO)usr1.$(OBJ_EXT) -module $(DPO)
$(DPO)usr2.$(OBJ_EXT) : usr2.f \
            $(DPO)USR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr2.f  -o $(DPO)usr2.$(OBJ_EXT) -module $(DPO)
$(DPO)usr3.$(OBJ_EXT) : usr3.f \
            $(DPO)USR.mod \
            $(DPO)FLDVAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr3.f  -o $(DPO)usr3.$(OBJ_EXT) -module $(DPO)
$(DPO)usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f \
            $(DPO)USR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_init_namelist.f  -o $(DPO)usr_init_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)usr_rates.$(OBJ_EXT) : usr_rates.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RXNS.mod \
            $(DPO)ENERGY.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)USR.mod \
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
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)BC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) utilities.f  -o $(DPO)utilities.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_u_g.$(OBJ_EXT) : vavg_u_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)MFLUX.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_u_g.f  -o $(DPO)vavg_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_u_s.$(OBJ_EXT) : vavg_u_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_u_s.f  -o $(DPO)vavg_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_v_g.$(OBJ_EXT) : vavg_v_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)MFLUX.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_v_g.f  -o $(DPO)vavg_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_v_s.$(OBJ_EXT) : vavg_v_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_v_s.f  -o $(DPO)vavg_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_w_g.$(OBJ_EXT) : vavg_w_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)MFLUX.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_w_g.f  -o $(DPO)vavg_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)vavg_w_s.$(OBJ_EXT) : vavg_w_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) vavg_w_s.f  -o $(DPO)vavg_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)vf_gs_x.$(OBJ_EXT) : vf_gs_x.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DRAG.mod \
            $(DPO)DISCRETELEMENT.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vf_gs_x.f  -o $(DPO)vf_gs_x.$(OBJ_EXT) -module $(DPO)
$(DPO)vf_gs_y.$(OBJ_EXT) : vf_gs_y.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DRAG.mod \
            $(DPO)DISCRETELEMENT.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vf_gs_y.f  -o $(DPO)vf_gs_y.$(OBJ_EXT) -module $(DPO)
$(DPO)vf_gs_z.$(OBJ_EXT) : vf_gs_z.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DRAG.mod \
            $(DPO)DISCRETELEMENT.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vf_gs_z.f  -o $(DPO)vf_gs_z.$(OBJ_EXT) -module $(DPO)
$(DPO)vtc_scalar.$(OBJ_EXT) : vtc_scalar.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)COMPAR.mod \
            $(DPO)KINTHEORY.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) vtc_scalar.f  -o $(DPO)vtc_scalar.$(OBJ_EXT) -module $(DPO)
$(DPO)write_ab_m.$(OBJ_EXT) : write_ab_m.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) write_ab_m.f  -o $(DPO)write_ab_m.$(OBJ_EXT) -module $(DPO)
$(DPO)write_ab_m_var.$(OBJ_EXT) : write_ab_m_var.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) write_ab_m_var.f  -o $(DPO)write_ab_m_var.$(OBJ_EXT) -module $(DPO)
$(DPO)write_error.$(OBJ_EXT) : write_error.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_error.f  -o $(DPO)write_error.$(OBJ_EXT) -module $(DPO)
$(DPO)write_header.$(OBJ_EXT) : write_header.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_header.f  -o $(DPO)write_header.$(OBJ_EXT) -module $(DPO)
$(DPO)write_out0.$(OBJ_EXT) : write_out0.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)IC.mod \
            $(DPO)BC.mod \
            $(DPO)IS.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)SCALES.mod \
            $(DPO)SCALARS.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RXNS.mod \
            $(DPO)MFIX_PIC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) write_out0.f  -o $(DPO)write_out0.$(OBJ_EXT) -module $(DPO)
$(DPO)write_out1.$(OBJ_EXT) : write_out1.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)SCALARS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)RXNS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_out1.f  -o $(DPO)write_out1.$(OBJ_EXT) -module $(DPO)
$(DPO)write_out3.$(OBJ_EXT) : write_out3.f \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_out3.f  -o $(DPO)write_out3.$(OBJ_EXT) -module $(DPO)
$(DPO)write_res0.$(OBJ_EXT) : write_res0.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)IC.mod \
            $(DPO)IS.mod \
            $(DPO)BC.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)SCALES.mod \
            $(DPO)SCALARS.mod \
            $(DPO)RXNS.mod \
            $(DPO)UR_FACS.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)CDIST.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)STIFF_CHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_res0.f  -o $(DPO)write_res0.$(OBJ_EXT) -module $(DPO)
$(DPO)write_res1.$(OBJ_EXT) : write_res1.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)SCALARS.mod \
            $(DPO)RXNS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)ENERGY.mod \
            $(DPO)CDIST.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)MFIX_NETCDF.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_res1.f  -o $(DPO)write_res1.$(OBJ_EXT) -module $(DPO)
$(DPO)write_spx0.$(OBJ_EXT) : write_spx0.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)FUNITS.mod \
            $(DPO)CDIST.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_spx0.f  -o $(DPO)write_spx0.$(OBJ_EXT) -module $(DPO)
$(DPO)write_spx1.$(OBJ_EXT) : write_spx1.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)FUNITS.mod \
            $(DPO)SCALARS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)RXNS.mod \
            $(DPO)CDIST.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)MFIX_NETCDF.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_spx1.f  -o $(DPO)write_spx1.$(OBJ_EXT) -module $(DPO)
$(DPO)write_table.$(OBJ_EXT) : write_table.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_table.f  -o $(DPO)write_table.$(OBJ_EXT) -module $(DPO)
$(DPO)write_usr0.$(OBJ_EXT) : write_usr0.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_usr0.f  -o $(DPO)write_usr0.$(OBJ_EXT) -module $(DPO)
$(DPO)write_usr1.$(OBJ_EXT) : write_usr1.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) write_usr1.f  -o $(DPO)write_usr1.$(OBJ_EXT) -module $(DPO)
$(DPO)xerbla.$(OBJ_EXT) : xerbla.f \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) xerbla.f  -o $(DPO)xerbla.$(OBJ_EXT) -module $(DPO)
$(DPO)zero_array.$(OBJ_EXT) : zero_array.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) zero_array.f  -o $(DPO)zero_array.$(OBJ_EXT) -module $(DPO)
$(DPO)zero_norm_vel.$(OBJ_EXT) : zero_norm_vel.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) zero_norm_vel.f  -o $(DPO)zero_norm_vel.$(OBJ_EXT) -module $(DPO)
$(DPO)allocate_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/allocate_cut_cell_arrays.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)INDICES.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)STL.mod \
            $(DPO)DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/allocate_cut_cell_arrays.f  -o $(DPO)allocate_cut_cell_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)allocate_dummy_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/allocate_dummy_cut_cell_arrays.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)INDICES.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/allocate_dummy_cut_cell_arrays.f  -o $(DPO)allocate_dummy_cut_cell_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_vort_out.$(OBJ_EXT) : ./cartesian_grid/calc_vort_out.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/calc_vort_out.f  -o $(DPO)calc_vort_out.$(OBJ_EXT) -module $(DPO)
$(DPO)cartesian_grid_init_namelist.$(OBJ_EXT) : ./cartesian_grid/cartesian_grid_init_namelist.f \
            $(DPO)PARAM1.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)POLYGON.mod \
            $(DPO)VTK.mod \
            $(DPO)PROGRESS_BAR.mod \
            $(DPO)DASHBOARD.mod \
            $(DPO)STL.mod \
            cartesian_grid/cartesian_grid_namelist.inc                  
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cartesian_grid_init_namelist.f  -o $(DPO)cartesian_grid_init_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_set_bc0.$(OBJ_EXT) : ./cartesian_grid/CG_set_bc0.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)BC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)FUNITS.mod \
            $(DPO)SCALES.mod \
            $(DPO)SCALARS.mod \
            $(DPO)BOUNDFUNIJK.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            sc_p_g1.inc                                                  \
            function.inc                                                 \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_set_bc0.f  -o $(DPO)CG_set_bc0.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_set_outflow.$(OBJ_EXT) : ./cartesian_grid/CG_set_outflow.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)BC.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)SCALARS.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MFLUX.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_set_outflow.f  -o $(DPO)CG_set_outflow.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_u_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_u_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_G.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_u_g.f  -o $(DPO)CG_source_u_g.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_u_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_u_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_S.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_u_s.f  -o $(DPO)CG_source_u_s.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_v_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_v_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_G.mod \
            $(DPO)BC.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_v_g.f  -o $(DPO)CG_source_v_g.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_v_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_v_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_S.mod \
            $(DPO)BC.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_v_s.f  -o $(DPO)CG_source_v_s.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_w_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_w_g.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_G.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_w_g.f  -o $(DPO)CG_source_w_g.$(OBJ_EXT) -module $(DPO)
$(DPO)CG_source_w_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_w_s.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_S.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)KINTHEORY.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)DRAG.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_w_s.f  -o $(DPO)CG_source_w_s.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_cartesian.$(OBJ_EXT) : ./cartesian_grid/check_data_cartesian.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)SCALARS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)BC.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)VTK.mod \
            $(DPO)POLYGON.mod \
            $(DPO)DASHBOARD.mod \
            $(DPO)STL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)SCALES.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)PS.mod \
            $(DPO)GRIDMAP.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/check_data_cartesian.f  -o $(DPO)check_data_cartesian.$(OBJ_EXT) -module $(DPO)
$(DPO)cut_cell_preprocessing.$(OBJ_EXT) : ./cartesian_grid/cut_cell_preprocessing.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)VTK.mod \
            $(DPO)CDIST.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)POLYGON.mod \
            $(DPO)STL.mod \
            $(DPO)STL.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cut_cell_preprocessing.f  -o $(DPO)cut_cell_preprocessing.$(OBJ_EXT) -module $(DPO)
$(DPO)deallocate_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/deallocate_cut_cell_arrays.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)INDICES.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/deallocate_cut_cell_arrays.f  -o $(DPO)deallocate_cut_cell_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)define_quadrics.$(OBJ_EXT) : ./cartesian_grid/define_quadrics.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)VTK.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/define_quadrics.f  -o $(DPO)define_quadrics.$(OBJ_EXT) -module $(DPO)
$(DPO)dmp_cartesian.$(OBJ_EXT) : ./cartesian_grid/dmp_cartesian.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/dmp_cartesian.f  -o $(DPO)dmp_cartesian.$(OBJ_EXT) -module $(DPO)
$(DPO)eval_usr_fct.$(OBJ_EXT) : ./cartesian_grid/eval_usr_fct.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/eval_usr_fct.f  -o $(DPO)eval_usr_fct.$(OBJ_EXT) -module $(DPO)
$(DPO)get_alpha.$(OBJ_EXT) : ./cartesian_grid/get_alpha.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)BC.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_alpha.f  -o $(DPO)get_alpha.$(OBJ_EXT) -module $(DPO)
$(DPO)get_connectivity.$(OBJ_EXT) : ./cartesian_grid/get_connectivity.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)POLYGON.mod \
            $(DPO)STL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VTK.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_connectivity.f  -o $(DPO)get_connectivity.$(OBJ_EXT) -module $(DPO)
$(DPO)get_cut_cell_flags.$(OBJ_EXT) : ./cartesian_grid/get_cut_cell_flags.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)VTK.mod \
            $(DPO)POLYGON.mod \
            $(DPO)STL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)SCALARS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_cut_cell_flags.f  -o $(DPO)get_cut_cell_flags.$(OBJ_EXT) -module $(DPO)
$(DPO)get_cut_cell_volume_area.$(OBJ_EXT) : ./cartesian_grid/get_cut_cell_volume_area.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)POLYGON.mod \
            $(DPO)STL.mod \
            $(DPO)BC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_cut_cell_volume_area.f  -o $(DPO)get_cut_cell_volume_area.$(OBJ_EXT) -module $(DPO)
$(DPO)get_delh.$(OBJ_EXT) : ./cartesian_grid/get_delh.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_delh.f  -o $(DPO)get_delh.$(OBJ_EXT) -module $(DPO)
$(DPO)get_master.$(OBJ_EXT) : ./cartesian_grid/get_master.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)BC.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_master.f  -o $(DPO)get_master.$(OBJ_EXT) -module $(DPO)
$(DPO)get_poly_data.$(OBJ_EXT) : ./cartesian_grid/get_poly_data.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)SCALARS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)RXNS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PROGRESS_BAR.mod \
            $(DPO)POLYGON.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_poly_data.f  -o $(DPO)get_poly_data.$(OBJ_EXT) -module $(DPO)
$(DPO)get_stl_data.$(OBJ_EXT) : ./cartesian_grid/get_stl_data.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)SCALARS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)RXNS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PROGRESS_BAR.mod \
            $(DPO)STL.mod \
            $(DPO)VTK.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)BC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)STL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_stl_data.f  -o $(DPO)get_stl_data.$(OBJ_EXT) -module $(DPO)
$(DPO)set_Odxyz.$(OBJ_EXT) : ./cartesian_grid/set_Odxyz.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)VTK.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/set_Odxyz.f  -o $(DPO)set_Odxyz.$(OBJ_EXT) -module $(DPO)
$(DPO)update_dashboard.$(OBJ_EXT) : ./cartesian_grid/update_dashboard.f \
            $(DPO)COMPAR.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)RUN.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)TIME_CPU.mod \
            $(DPO)RESIDUAL.mod \
            $(DPO)DASHBOARD.mod \
            $(DPO)VTK.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/update_dashboard.f  -o $(DPO)update_dashboard.$(OBJ_EXT) -module $(DPO)
$(DPO)vtk_out.$(OBJ_EXT) : ./cartesian_grid/vtk_out.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QUADRIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PGCOR.mod \
            $(DPO)VTK.mod \
            $(DPO)RXNS.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)SCALARS.mod \
            $(DPO)STL.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PARALLEL_MPI.mod \
            $(DPO)PSCOR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)CDIST.mod \
            $(DPO)POLYGON.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/vtk_out.f  -o $(DPO)vtk_out.$(OBJ_EXT) -module $(DPO)
$(DPO)write_progress_bar.$(OBJ_EXT) : ./cartesian_grid/write_progress_bar.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)SCALARS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)RXNS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)PROGRESS_BAR.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)SENDRECV.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/write_progress_bar.f  -o $(DPO)write_progress_bar.$(OBJ_EXT) -module $(DPO)
$(DPO)check_axis.$(OBJ_EXT) : ./check_data/check_axis.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_axis.f  -o $(DPO)check_axis.$(OBJ_EXT) -module $(DPO)
$(DPO)check_bc_geometry.$(OBJ_EXT) : ./check_data/check_bc_geometry.f \
            $(DPO)BC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_bc_geometry.f  -o $(DPO)check_bc_geometry.$(OBJ_EXT) -module $(DPO)
$(DPO)check_bc_inflow.$(OBJ_EXT) : ./check_data/check_bc_inflow.f \
            $(DPO)RUN.mod \
            $(DPO)SCALARS.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)BC.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_bc_inflow.f  -o $(DPO)check_bc_inflow.$(OBJ_EXT) -module $(DPO)
$(DPO)check_bc_outflow.$(OBJ_EXT) : ./check_data/check_bc_outflow.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)BC.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)SCALARS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_bc_outflow.f  -o $(DPO)check_bc_outflow.$(OBJ_EXT) -module $(DPO)
$(DPO)check_bc_walls.$(OBJ_EXT) : ./check_data/check_bc_walls.f \
            $(DPO)RUN.mod \
            $(DPO)BC.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)SCALARS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_bc_walls.f  -o $(DPO)check_bc_walls.$(OBJ_EXT) -module $(DPO)
$(DPO)check_boundary_conditions.$(OBJ_EXT) : ./check_data/check_boundary_conditions.f \
            $(DPO)PHYSPROP.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            $(DPO)BC.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARAM.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_boundary_conditions.f  -o $(DPO)check_boundary_conditions.$(OBJ_EXT) -module $(DPO)
$(DPO)check_chemical_rxns.$(OBJ_EXT) : ./check_data/check_chemical_rxns.f \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARSE.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)RXNS.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_chemical_rxns.f  -o $(DPO)check_chemical_rxns.$(OBJ_EXT) -module $(DPO)
$(DPO)check_dmp_prereqs.$(OBJ_EXT) : ./check_data/check_dmp_prereqs.f \
            $(DPO)COMPAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PARAM1.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_dmp_prereqs.f  -o $(DPO)check_dmp_prereqs.$(OBJ_EXT) -module $(DPO)
$(DPO)check_gas_phase.$(OBJ_EXT) : ./check_data/check_gas_phase.f \
            $(DPO)RUN.mod \
            $(DPO)RXNS.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PARAM1.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_gas_phase.f  -o $(DPO)check_gas_phase.$(OBJ_EXT) -module $(DPO)
$(DPO)check_geometry.$(OBJ_EXT) : ./check_data/check_geometry.f \
            $(DPO)GEOMETRY.mod \
            $(DPO)BC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARAM.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_geometry.f  -o $(DPO)check_geometry.$(OBJ_EXT) -module $(DPO)
$(DPO)check_geometry_prereqs.$(OBJ_EXT) : ./check_data/check_geometry_prereqs.f \
            $(DPO)GEOMETRY.mod \
            $(DPO)BC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARAM.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_geometry_prereqs.f  -o $(DPO)check_geometry_prereqs.$(OBJ_EXT) -module $(DPO)
$(DPO)check_ic_common_discrete.$(OBJ_EXT) : ./check_data/check_ic_common_discrete.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)IC.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PARAM.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_ic_common_discrete.f  -o $(DPO)check_ic_common_discrete.$(OBJ_EXT) -module $(DPO)
$(DPO)check_initial_conditions_dem.$(OBJ_EXT) : ./check_data/check_initial_conditions_dem.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_initial_conditions_dem.f  -o $(DPO)check_initial_conditions_dem.$(OBJ_EXT) -module $(DPO)
$(DPO)check_initial_conditions.$(OBJ_EXT) : ./check_data/check_initial_conditions.f \
            $(DPO)IC.mod \
            $(DPO)RUN.mod \
            $(DPO)PARAM.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)SCALARS.mod \
            $(DPO)DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_initial_conditions.f  -o $(DPO)check_initial_conditions.$(OBJ_EXT) -module $(DPO)
$(DPO)check_initial_conditions_mppic.$(OBJ_EXT) : ./check_data/check_initial_conditions_mppic.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)IC.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARAM.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)ERROR_MANAGER.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_initial_conditions_mppic.f  -o $(DPO)check_initial_conditions_mppic.$(OBJ_EXT) -module $(DPO)
$(DPO)check_internal_surfaces.$(OBJ_EXT) : ./check_data/check_internal_surfaces.f \
            $(DPO)IS.mod \
            $(DPO)PARAM.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)INDICES.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_internal_surfaces.f  -o $(DPO)check_internal_surfaces.$(OBJ_EXT) -module $(DPO)
$(DPO)check_numerics.$(OBJ_EXT) : ./check_data/check_numerics.f \
            $(DPO)RUN.mod \
            $(DPO)LEQSOL.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_numerics.f  -o $(DPO)check_numerics.$(OBJ_EXT) -module $(DPO)
$(DPO)check_output_control.$(OBJ_EXT) : ./check_data/check_output_control.f \
            $(DPO)OUTPUT.mod \
            $(DPO)RUN.mod \
            $(DPO)RXNS.mod \
            $(DPO)PARAM1.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_output_control.f  -o $(DPO)check_output_control.$(OBJ_EXT) -module $(DPO)
$(DPO)check_point_sources.$(OBJ_EXT) : ./check_data/check_point_sources.f \
            $(DPO)PS.mod \
            $(DPO)PARAM.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_point_sources.f  -o $(DPO)check_point_sources.$(OBJ_EXT) -module $(DPO)
$(DPO)check_run_control.$(OBJ_EXT) : ./check_data/check_run_control.f \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PARAM1.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_run_control.f  -o $(DPO)check_run_control.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_common_all.$(OBJ_EXT) : ./check_data/check_solids_common_all.f \
            $(DPO)RXNS.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)RUN.mod \
            $(DPO)DRAG.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)MFIX_PIC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_common_all.f  -o $(DPO)check_solids_common_all.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_common_discrete.$(OBJ_EXT) : ./check_data/check_solids_common_discrete.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DESGRID.mod \
            $(DPO)RUN.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)PARAM1.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)GEOMETRY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_common_discrete.f  -o $(DPO)check_solids_common_discrete.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_continuum.$(OBJ_EXT) : ./check_data/check_solids_continuum.f \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)SCALARS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)RXNS.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_continuum.f  -o $(DPO)check_solids_continuum.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_des.$(OBJ_EXT) : ./check_data/check_solids_des.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)PARAM1.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_des.f  -o $(DPO)check_solids_des.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_model_prereqs.$(OBJ_EXT) : ./check_data/check_solids_model_prereqs.f \
            $(DPO)RUN.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)SCALARS.mod \
            $(DPO)PARAM.mod \
            $(DPO)ERROR_MANAGER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_model_prereqs.f  -o $(DPO)check_solids_model_prereqs.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_mppic.$(OBJ_EXT) : ./check_data/check_solids_mppic.f \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FUNITS.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)ERROR_MANAGER.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_mppic.f  -o $(DPO)check_solids_mppic.$(OBJ_EXT) -module $(DPO)
$(DPO)check_solids_phases.$(OBJ_EXT) : ./check_data/check_solids_phases.f \
            $(DPO)RUN.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)PHYSPROP.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./check_data/check_solids_phases.f  -o $(DPO)check_solids_phases.$(OBJ_EXT) -module $(DPO)
$(DPO)check_data_odepack.$(OBJ_EXT) : ./chem/check_data_odepack.f \
            $(DPO)FUNITS.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)RXNS.mod \
            $(DPO)STIFF_CHEM.mod \
            $(DPO)COMPAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)STIFF_CHEM_DBG.mod \
            $(DPO)STIFF_CHEM_STATS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/check_data_odepack.f  -o $(DPO)check_data_odepack.$(OBJ_EXT) -module $(DPO)
$(DPO)stiff_chem_rrates.$(OBJ_EXT) : ./chem/stiff_chem_rrates.f \
            $(DPO)STIFF_CHEM_MAPS.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)STIFF_CHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/stiff_chem_rrates.f  -o $(DPO)stiff_chem_rrates.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_force_dem.$(OBJ_EXT) : ./des/calc_force_dem.f \
            $(DPO)RUN.mod \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_force_dem.f  -o $(DPO)calc_force_dem.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_force_dem_stl.$(OBJ_EXT) : ./des/calc_force_dem_stl.f \
            $(DPO)RUN.mod \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)FUNITS.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)SOFTSPRING_FUNCS_CUTCELL.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_force_dem_stl.f  -o $(DPO)calc_force_dem_stl.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_rrate_des.$(OBJ_EXT) : ./des/calc_rrate_des.f \
            $(DPO)COMPAR.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)ENERGY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INTERPOLATION.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)RXNS.mod \
            $(DPO)USR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_rrate_des.f  -o $(DPO)calc_rrate_des.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_thermo_des.$(OBJ_EXT) : ./des/calc_thermo_des.f \
            $(DPO)COMPAR.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)INTERPOLATION.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_thermo_des.f  -o $(DPO)calc_thermo_des.$(OBJ_EXT) -module $(DPO)
$(DPO)cfassign.$(OBJ_EXT) : ./des/cfassign.f \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)RUN.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)CUTCELL.mod \
            b_force1.inc                                                 \
            b_force2.inc                                                 \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfassign.f  -o $(DPO)cfassign.$(OBJ_EXT) -module $(DPO)
$(DPO)cffctowall.$(OBJ_EXT) : ./des/cffctowall.f \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffctowall.f  -o $(DPO)cffctowall.$(OBJ_EXT) -module $(DPO)
$(DPO)cffctow.$(OBJ_EXT) : ./des/cffctow.f \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffctow.f  -o $(DPO)cffctow.$(OBJ_EXT) -module $(DPO)
$(DPO)cfnewvalues.$(OBJ_EXT) : ./des/cfnewvalues.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DES_BC.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)MPPIC_WALLBC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)RANDOMNO.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfnewvalues.f  -o $(DPO)cfnewvalues.$(OBJ_EXT) -module $(DPO)
$(DPO)cfrelvel.$(OBJ_EXT) : ./des/cfrelvel.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfrelvel.f  -o $(DPO)cfrelvel.$(OBJ_EXT) -module $(DPO)
$(DPO)cfslide.$(OBJ_EXT) : ./des/cfslide.f \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslide.f  -o $(DPO)cfslide.$(OBJ_EXT) -module $(DPO)
$(DPO)cfslidewall.$(OBJ_EXT) : ./des/cfslidewall.f \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslidewall.f  -o $(DPO)cfslidewall.$(OBJ_EXT) -module $(DPO)
$(DPO)cfupdateold.$(OBJ_EXT) : ./des/cfupdateold.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfupdateold.f  -o $(DPO)cfupdateold.$(OBJ_EXT) -module $(DPO)
$(DPO)cfwallcontact.$(OBJ_EXT) : ./des/cfwallcontact.f \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DES_BC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallcontact.f  -o $(DPO)cfwallcontact.$(OBJ_EXT) -module $(DPO)
$(DPO)cfwallposvel.$(OBJ_EXT) : ./des/cfwallposvel.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DES_BC.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)MATRIX.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DRAG.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallposvel.f  -o $(DPO)cfwallposvel.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_bc.$(OBJ_EXT) : ./des/check_des_bc.f \
            $(DPO)CONSTANT.mod \
            $(DPO)DES_BC.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_bc.f  -o $(DPO)check_des_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_data.$(OBJ_EXT) : ./des/check_des_data.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)DESGRID.mod \
            $(DPO)FUNITS.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_data.f  -o $(DPO)check_des_data.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_hybrid.$(OBJ_EXT) : ./des/check_des_hybrid.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_hybrid.f  -o $(DPO)check_des_hybrid.$(OBJ_EXT) -module $(DPO)
$(DPO)check_des_rxns.$(OBJ_EXT) : ./des/check_des_rxns.f \
            $(DPO)COMPAR.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)RUN.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARSE.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RXNS.mod \
            $(DPO)STIFF_CHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_rxns.f  -o $(DPO)check_des_rxns.$(OBJ_EXT) -module $(DPO)
$(DPO)des_allocate_arrays.$(OBJ_EXT) : ./des/des_allocate_arrays.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)INDICES.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DES_BC.mod \
            $(DPO)FUNITS.mod \
            $(DPO)DESGRID.mod \
            $(DPO)DESMPI.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)ERROR_MANAGER.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_allocate_arrays.f  -o $(DPO)des_allocate_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)des_check_particle.$(OBJ_EXT) : ./des/des_check_particle.f \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DES_BC.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_check_particle.f  -o $(DPO)des_check_particle.$(OBJ_EXT) -module $(DPO)
$(DPO)des_cluster_identification.$(OBJ_EXT) : ./des/des_cluster_identification.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            $(DPO)DES_CLUSTER.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_cluster_identification.f  -o $(DPO)des_cluster_identification.$(OBJ_EXT) -module $(DPO)
$(DPO)des_functions.$(OBJ_EXT) : ./des/des_functions.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_functions.f  -o $(DPO)des_functions.$(OBJ_EXT) -module $(DPO)
$(DPO)des_granular_temperature.$(OBJ_EXT) : ./des/des_granular_temperature.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DES_BC.mod \
            $(DPO)FLDVAR.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_granular_temperature.f  -o $(DPO)des_granular_temperature.$(OBJ_EXT) -module $(DPO)
$(DPO)des_init_arrays.$(OBJ_EXT) : ./des/des_init_arrays.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)INDICES.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DES_BC.mod \
            $(DPO)RUN.mod \
            $(DPO)DESGRID.mod \
            $(DPO)DESMPI.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_arrays.f  -o $(DPO)des_init_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)des_init_bc.$(OBJ_EXT) : ./des/des_init_bc.f \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DES_BC.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_bc.f  -o $(DPO)des_init_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)des_init_namelist.$(OBJ_EXT) : ./des/des_init_namelist.f \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)DES_BC.mod \
            $(DPO)DES_IC.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod \
            des/desnamelist.inc                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_namelist.f  -o $(DPO)des_init_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)des_mass_inlet.$(OBJ_EXT) : ./des/des_mass_inlet.f \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DES_BC.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DESGRID.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod \
            function.inc                                                 \
            des/desgrid_functions.inc                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_mass_inlet.f  -o $(DPO)des_mass_inlet.$(OBJ_EXT) -module $(DPO)
$(DPO)des_physical_prop.$(OBJ_EXT) : ./des/des_physical_prop.f \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_physical_prop.f  -o $(DPO)des_physical_prop.$(OBJ_EXT) -module $(DPO)
$(DPO)des_reaction_model.$(OBJ_EXT) : ./des/des_reaction_model.f \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PARAM1.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_reaction_model.f  -o $(DPO)des_reaction_model.$(OBJ_EXT) -module $(DPO)
$(DPO)des_rrates0.$(OBJ_EXT) : ./des/des_rrates0.f \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)ENERGY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARSE.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)RXNS.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)USR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_rrates0.f  -o $(DPO)des_rrates0.$(OBJ_EXT) -module $(DPO)
$(DPO)des_set_ic.$(OBJ_EXT) : ./des/des_set_ic.f \
            $(DPO)COMPAR.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DES_IC.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)FUNITS.mod \
            $(DPO)PHYSPROP.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_set_ic.f  -o $(DPO)des_set_ic.$(OBJ_EXT) -module $(DPO)
$(DPO)des_thermo_cond.$(OBJ_EXT) : ./des/des_thermo_cond.f \
            $(DPO)CONSTANT.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)PHYSPROP.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_cond.f  -o $(DPO)des_thermo_cond.$(OBJ_EXT) -module $(DPO)
$(DPO)des_thermo_conv.$(OBJ_EXT) : ./des/des_thermo_conv.f \
            $(DPO)CONSTANT.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INTERPOLATION.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)USR.mod \
            $(DPO)COMPAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_conv.f  -o $(DPO)des_thermo_conv.$(OBJ_EXT) -module $(DPO)
$(DPO)des_thermo_newvalues.$(OBJ_EXT) : ./des/des_thermo_newvalues.f \
            $(DPO)COMPAR.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_newvalues.f  -o $(DPO)des_thermo_newvalues.$(OBJ_EXT) -module $(DPO)
$(DPO)des_thermo_rad.$(OBJ_EXT) : ./des/des_thermo_rad.f \
            $(DPO)CONSTANT.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_rad.f  -o $(DPO)des_thermo_rad.$(OBJ_EXT) -module $(DPO)
$(DPO)des_time_march.$(OBJ_EXT) : ./des/des_time_march.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PGCOR.mod \
            $(DPO)PSCOR.mod \
            $(DPO)CONT.mod \
            $(DPO)TAU_G.mod \
            $(DPO)TAU_S.mod \
            $(DPO)VISC_G.mod \
            $(DPO)VISC_S.mod \
            $(DPO)FUNITS.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)SCALARS.mod \
            $(DPO)DRAG.mod \
            $(DPO)RXNS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TIME_CPU.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DES_BC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)INTERPOLATION.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_time_march.f  -o $(DPO)des_time_march.$(OBJ_EXT) -module $(DPO)
$(DPO)des_wallbc_preprocessing.$(OBJ_EXT) : ./des/des_wallbc_preprocessing.f \
            $(DPO)PARAM1.mod \
            $(DPO)FUNITS.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)BC.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DESMPI.mod \
            $(DPO)CDIST.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)DISCRETELEMENT.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_wallbc_preprocessing.f  -o $(DPO)des_wallbc_preprocessing.$(OBJ_EXT) -module $(DPO)
$(DPO)drag_fgs.$(OBJ_EXT) : ./des/drag_fgs.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DRAG.mod \
            $(DPO)DESMPI.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)INTERPOLATION.mod \
            $(DPO)UR_FACS.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/drag_fgs.f  -o $(DPO)drag_fgs.$(OBJ_EXT) -module $(DPO)
$(DPO)gas_drag.$(OBJ_EXT) : ./des/gas_drag.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_G.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DRAG.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/gas_drag.f  -o $(DPO)gas_drag.$(OBJ_EXT) -module $(DPO)
$(DPO)generate_particle_config.$(OBJ_EXT) : ./des/generate_particle_config.f \
            $(DPO)MFIX_PIC.mod \
            $(DPO)DES_LINKED_LIST_DATA.mod \
            $(DPO)DES_LINKED_LIST_FUNCS.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)ERROR_MANAGER.mod \
            $(DPO)DES_LINKED_LIST_DATA.mod \
            $(DPO)RUN.mod \
            $(DPO)SOFTSPRING_FUNCS_CUTCELL.mod \
            $(DPO)INDICES.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)MPPIC_WALLBC.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)IC.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARAM.mod \
            $(DPO)RANDOMNO.mod \
            $(DPO)FUNITS.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)PARALLEL.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/generate_particle_config.f  -o $(DPO)generate_particle_config.$(OBJ_EXT) -module $(DPO)
$(DPO)make_arrays_des.$(OBJ_EXT) : ./des/make_arrays_des.f \
            $(DPO)PARAM1.mod \
            $(DPO)FUNITS.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)DESMPI.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)DES_IC.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_STL_FUNCTIONS.mod \
            $(DPO)CDIST.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/make_arrays_des.f  -o $(DPO)make_arrays_des.$(OBJ_EXT) -module $(DPO)
$(DPO)mppic_routines.$(OBJ_EXT) : ./des/mppic_routines.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DES_BC.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)MPPIC_WALLBC.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)DESMPI.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INTERPOLATION.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mppic_routines.f  -o $(DPO)mppic_routines.$(OBJ_EXT) -module $(DPO)
$(DPO)neighbour.$(OBJ_EXT) : ./des/neighbour.f \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DESGRID.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/neighbour.f  -o $(DPO)neighbour.$(OBJ_EXT) -module $(DPO)
$(DPO)nsquare.$(OBJ_EXT) : ./des/nsquare.f \
            $(DPO)PARAM1.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)DES_BC.mod \
            $(DPO)DES_THERMO.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/nsquare.f  -o $(DPO)nsquare.$(OBJ_EXT) -module $(DPO)
$(DPO)particles_in_cell.$(OBJ_EXT) : ./des/particles_in_cell.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)DESGRID.mod \
            $(DPO)DESMPI.mod \
            $(DPO)CUTCELL.mod \
            $(DPO)MFIX_PIC.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)BC.mod \
            $(DPO)DRAG.mod \
            $(DPO)INTERPOLATION.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)FUNITS.mod \
            $(DPO)DESMPI_WRAPPER.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/particles_in_cell.f  -o $(DPO)particles_in_cell.$(OBJ_EXT) -module $(DPO)
$(DPO)read_des_restart.$(OBJ_EXT) : ./des/read_des_restart.f \
            $(DPO)PARAM1.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            $(DPO)DES_BC.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DESMPI.mod \
            $(DPO)MACHINE.mod \
            $(DPO)CDIST.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/read_des_restart.f  -o $(DPO)read_des_restart.$(OBJ_EXT) -module $(DPO)
$(DPO)solid_drag.$(OBJ_EXT) : ./des/solid_drag.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)MATRIX.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/solid_drag.f  -o $(DPO)solid_drag.$(OBJ_EXT) -module $(DPO)
$(DPO)usr0_des.$(OBJ_EXT) : ./des/usr0_des.f \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            $(DPO)USR.mod \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr0_des.f  -o $(DPO)usr0_des.$(OBJ_EXT) -module $(DPO)
$(DPO)usr1_des.$(OBJ_EXT) : ./des/usr1_des.f \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            $(DPO)USR.mod \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr1_des.f  -o $(DPO)usr1_des.$(OBJ_EXT) -module $(DPO)
$(DPO)usr2_des.$(OBJ_EXT) : ./des/usr2_des.f \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            $(DPO)USR.mod \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr2_des.f  -o $(DPO)usr2_des.$(OBJ_EXT) -module $(DPO)
$(DPO)usr3_des.$(OBJ_EXT) : ./des/usr3_des.f \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            $(DPO)USR.mod \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr3_des.f  -o $(DPO)usr3_des.$(OBJ_EXT) -module $(DPO)
$(DPO)usr4_des.$(OBJ_EXT) : ./des/usr4_des.f \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)USR.mod \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr4_des.f  -o $(DPO)usr4_des.$(OBJ_EXT) -module $(DPO)
$(DPO)usr_rates_des.$(OBJ_EXT) : ./des/usr_rates_des.f \
            $(DPO)COMPAR.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)ENERGY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)FUNITS.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)USR.mod \
            species.inc                                                  \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/usr_rates_des.f  -o $(DPO)usr_rates_des.$(OBJ_EXT) -module $(DPO)
$(DPO)walledgecontact.$(OBJ_EXT) : ./des/walledgecontact.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)MATRIX.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DRAG.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/walledgecontact.f  -o $(DPO)walledgecontact.$(OBJ_EXT) -module $(DPO)
$(DPO)wallfacecontact.$(OBJ_EXT) : ./des/wallfacecontact.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)MATRIX.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DRAG.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/wallfacecontact.f  -o $(DPO)wallfacecontact.$(OBJ_EXT) -module $(DPO)
$(DPO)wallnodecontact.$(OBJ_EXT) : ./des/wallnodecontact.f \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)MATRIX.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)DRAG.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/wallnodecontact.f  -o $(DPO)wallnodecontact.$(OBJ_EXT) -module $(DPO)
$(DPO)write_des_data.$(OBJ_EXT) : ./des/write_des_data.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DES_BC.mod \
            $(DPO)MPI_UTILITY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DESMPI.mod \
            $(DPO)CDIST.mod \
            $(DPO)DES_THERMO.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/write_des_data.f  -o $(DPO)write_des_data.$(OBJ_EXT) -module $(DPO)
$(DPO)write_des_restart.$(OBJ_EXT) : ./des/write_des_restart.f \
            $(DPO)PARAM1.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)RUN.mod \
            $(DPO)DES_BC.mod \
            $(DPO)DES_RXNS.mod \
            $(DPO)DES_THERMO.mod \
            $(DPO)DESMPI.mod \
            $(DPO)MACHINE.mod \
            $(DPO)CDIST.mod \
            $(DPO)MPI_UTILITY.mod 
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
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)SCALARS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/source_population_eq.f  -o $(DPO)source_population_eq.$(OBJ_EXT) -module $(DPO)
$(DPO)usr_dqmom.$(OBJ_EXT) : ./dqmom/usr_dqmom.f \
            $(DPO)RUN.mod \
            $(DPO)SCALARS.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)USR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/usr_dqmom.f  -o $(DPO)usr_dqmom.$(OBJ_EXT) -module $(DPO)
$(DPO)adjust_eps_ghd.$(OBJ_EXT) : ./GhdTheory/adjust_eps_ghd.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)GHDTHEORY.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/adjust_eps_ghd.f  -o $(DPO)adjust_eps_ghd.$(OBJ_EXT) -module $(DPO)
$(DPO)bulk_viscosity.$(OBJ_EXT) : ./GhdTheory/bulk_viscosity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/bulk_viscosity.f  -o $(DPO)bulk_viscosity.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_d_ghd.$(OBJ_EXT) : ./GhdTheory/calc_d_ghd.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)SCALES.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_d_ghd.f  -o $(DPO)calc_d_ghd.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_external_forces.$(OBJ_EXT) : ./GhdTheory/calc_external_forces.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)DRAG.mod \
            $(DPO)BC.mod \
            $(DPO)SCALES.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            b_force1.inc                                                 \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_external_forces.f  -o $(DPO)calc_external_forces.$(OBJ_EXT) -module $(DPO)
$(DPO)calc_nflux.$(OBJ_EXT) : ./GhdTheory/calc_nflux.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)INDICES.mod \
            $(DPO)MFLUX.mod \
            $(DPO)COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_nflux.f  -o $(DPO)calc_nflux.$(OBJ_EXT) -module $(DPO)
$(DPO)chi_ij_GHD.$(OBJ_EXT) : ./GhdTheory/chi_ij_GHD.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/chi_ij_GHD.f  -o $(DPO)chi_ij_GHD.$(OBJ_EXT) -module $(DPO)
$(DPO)cooling_rate.$(OBJ_EXT) : ./GhdTheory/cooling_rate.f \
            $(DPO)COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/cooling_rate.f  -o $(DPO)cooling_rate.$(OBJ_EXT) -module $(DPO)
$(DPO)cooling_rate_tc.$(OBJ_EXT) : ./GhdTheory/cooling_rate_tc.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/cooling_rate_tc.f  -o $(DPO)cooling_rate_tc.$(OBJ_EXT) -module $(DPO)
$(DPO)dufour_coeff.$(OBJ_EXT) : ./GhdTheory/dufour_coeff.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/dufour_coeff.f  -o $(DPO)dufour_coeff.$(OBJ_EXT) -module $(DPO)
$(DPO)ghd.$(OBJ_EXT) : ./GhdTheory/ghd.f \
            $(DPO)DRAG.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ghd.f  -o $(DPO)ghd.$(OBJ_EXT) -module $(DPO)
$(DPO)ghdmassflux.$(OBJ_EXT) : ./GhdTheory/ghdmassflux.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)VISC_S.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)DRAG.mod \
            $(DPO)IS.mod \
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
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)MATRIX.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)DRAG.mod \
            $(DPO)FLDVAR.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/partial_elim_ghd.f  -o $(DPO)partial_elim_ghd.$(OBJ_EXT) -module $(DPO)
$(DPO)pressure.$(OBJ_EXT) : ./GhdTheory/pressure.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/pressure.f  -o $(DPO)pressure.$(OBJ_EXT) -module $(DPO)
$(DPO)shear_viscosity.$(OBJ_EXT) : ./GhdTheory/shear_viscosity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/shear_viscosity.f  -o $(DPO)shear_viscosity.$(OBJ_EXT) -module $(DPO)
$(DPO)source_ghd_granular_energy.$(OBJ_EXT) : ./GhdTheory/source_ghd_granular_energy.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)DRAG.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_S.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)TRACE.mod \
            $(DPO)INDICES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)COMPAR.mod \
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
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)VISC_S.mod \
            $(DPO)GHDTHEORY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)RUN.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)TOLERANC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/transport_coeff_ghd.f  -o $(DPO)transport_coeff_ghd.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_allocate_arrays.$(OBJ_EXT) : ./qmomk/qmomk_allocate_arrays.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)INDICES.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_allocate_arrays.f  -o $(DPO)qmomk_allocate_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_gas_drag.$(OBJ_EXT) : ./qmomk/qmomk_gas_drag.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)PARALLEL.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SCALES.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)VISC_G.mod \
            $(DPO)RXNS.mod \
            $(DPO)RUN.mod \
            $(DPO)TOLERANC.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)IS.mod \
            $(DPO)TAU_G.mod \
            $(DPO)BC.mod \
            $(DPO)COMPAR.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)DISCRETELEMENT.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)DRAG.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_gas_drag.f  -o $(DPO)qmomk_gas_drag.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_init_bc.$(OBJ_EXT) : ./qmomk/qmomk_init_bc.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)BC.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)QMOMK_QUADRATURE.mod \
            $(DPO)QMOMK_BC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_init_bc.f  -o $(DPO)qmomk_init_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_initial_conditions.$(OBJ_EXT) : ./qmomk/qmomk_initial_conditions.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)QMOMK_QUADRATURE.mod \
            $(DPO)QMOMK_PARAMETERS.mod \
            $(DPO)IC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_initial_conditions.f  -o $(DPO)qmomk_initial_conditions.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_init_namelist.$(OBJ_EXT) : ./qmomk/qmomk_init_namelist.f \
            $(DPO)PARAM1.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            qmomk/qmomknamelist.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_init_namelist.f  -o $(DPO)qmomk_init_namelist.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_make_arrays.$(OBJ_EXT) : ./qmomk/qmomk_make_arrays.f \
            $(DPO)PARAM1.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)FUNITS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_make_arrays.f  -o $(DPO)qmomk_make_arrays.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_read_restart.$(OBJ_EXT) : ./qmomk/qmomk_read_restart.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)CONT.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)INDICES.mod \
            $(DPO)RUN.mod \
            $(DPO)COMPAR.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)QMOMK_QUADRATURE.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_read_restart.f  -o $(DPO)qmomk_read_restart.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_set_bc.$(OBJ_EXT) : ./qmomk/qmomk_set_bc.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)COMPAR.mod \
            $(DPO)INDICES.mod \
            $(DPO)BC.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)QMOMK_QUADRATURE.mod \
            $(DPO)QMOMK_BC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_set_bc.f  -o $(DPO)qmomk_set_bc.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_time_march.$(OBJ_EXT) : ./qmomk/qmomk_time_march.f \
            $(DPO)PARAM.mod \
            $(DPO)PARAM1.mod \
            $(DPO)CONSTANT.mod \
            $(DPO)RUN.mod \
            $(DPO)OUTPUT.mod \
            $(DPO)PHYSPROP.mod \
            $(DPO)FLDVAR.mod \
            $(DPO)GEOMETRY.mod \
            $(DPO)CONT.mod \
            $(DPO)TAU_G.mod \
            $(DPO)TAU_S.mod \
            $(DPO)VISC_G.mod \
            $(DPO)VISC_S.mod \
            $(DPO)FUNITS.mod \
            $(DPO)VSHEAR.mod \
            $(DPO)SCALARS.mod \
            $(DPO)DRAG.mod \
            $(DPO)RXNS.mod \
            $(DPO)COMPAR.mod \
            $(DPO)TIME_CPU.mod \
            $(DPO)IS.mod \
            $(DPO)INDICES.mod \
            $(DPO)MATRIX.mod \
            $(DPO)SENDRECV.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)QMOMK_FLUXES.mod \
            $(DPO)QMOMK_QUADRATURE.mod \
            $(DPO)QMOMK_COLLISION.mod \
            $(DPO)QMOMK_PARAMETERS.mod \
            $(DPO)UR_FACS.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_time_march.f  -o $(DPO)qmomk_time_march.$(OBJ_EXT) -module $(DPO)
$(DPO)qmomk_write_restart.$(OBJ_EXT) : ./qmomk/qmomk_write_restart.f \
            $(DPO)PARAM1.mod \
            $(DPO)QMOM_KINETIC_EQUATION.mod \
            $(DPO)RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_write_restart.f  -o $(DPO)qmomk_write_restart.$(OBJ_EXT) -module $(DPO)
$(DPO)get_values.$(OBJ_EXT) : ./thermochemical/get_values.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./thermochemical/get_values.f  -o $(DPO)get_values.$(OBJ_EXT) -module $(DPO)
$(DPO)readTherm.$(OBJ_EXT) : ./thermochemical/readTherm.f \
            $(DPO)PHYSPROP.mod \
            $(DPO)DES_RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./thermochemical/readTherm.f  -o $(DPO)readTherm.$(OBJ_EXT) -module $(DPO)
