CC_FLAGS  = -c -g

.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
F_Obj1 = \
    any_more_data.$(OBJ_EXT) \
    calc_cell2.$(OBJ_EXT) \
    calc_corr_01.$(OBJ_EXT) \
    calc_corr_type_1.$(OBJ_EXT) \
    calc_distance.$(OBJ_EXT) \
    calc_dpody.$(OBJ_EXT) \
    calc_ep_g.$(OBJ_EXT) \
    calc_mu_s.$(OBJ_EXT) \
    calc_mw.$(OBJ_EXT) \
    calc_quantities.$(OBJ_EXT) \
    calc_ro_g.$(OBJ_EXT) \
    calc_vol.$(OBJ_EXT) \
    check_data_03.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    compare.$(OBJ_EXT) \
    eosg.$(OBJ_EXT) \
    error_routine.$(OBJ_EXT) \
    examine_data.$(OBJ_EXT) \
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
    make_upper_case.$(OBJ_EXT) \
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
    write_spx1.$(OBJ_EXT)

F_Obj2 = main_f.o

F_Obj3 = main_c.o               post_mfix.o

post_mfix : $(F_Obj1) $(F_Obj2)

  post_mfix : $(F_Obj1) $(F_Obj2)
	$(FORTRAN_CMD) $(LINK_FLAGS) $(F_Obj1) $(F_Obj2)                   \
           -o post_mfix

user:
	cc -c -g -I$(FORMS_DIR) main_c.c
	cc -c -g -I$(FORMS_DIR) post_mfix.c
	$(FORTRAN_CMD) $(LINK_FLAGS) $(F_Obj1) $(F_Obj3)                   \
            -L$(FORMS_DIR) -lforms -lX11 -lm \
           -o xpost_mfix

  
any_more_data.$(OBJ_EXT) : any_more_data.f \
            param.inc                                                    \
            param1.inc                                                  
calc_cell2.$(OBJ_EXT) : calc_cell2.f \
            param.inc                                                   
calc_corr_01.$(OBJ_EXT) : calc_corr_01.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            physprop.inc                                                 \
            geometry.inc                                                 \
            indices.inc                                                  \
            correl.inc                                                   \
            xforms.inc                                                   \
            function.inc                                                
calc_corr_type_1.$(OBJ_EXT) : calc_corr_type_1.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            run.inc                                                      \
            geometry.inc                                                 \
            indices.inc                                                  \
            post3d.inc                                                   \
            physprop.inc                                                 \
            correl.inc                                                   \
            xforms.inc                                                   \
            function.inc                                                
calc_distance.$(OBJ_EXT) : calc_distance.f
calc_dpody.$(OBJ_EXT) : calc_dpody.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            indices.inc                                                  \
            physprop.inc                                                 \
            geometry.inc                                                 \
            function.inc                                                
calc_ep_g.$(OBJ_EXT) : calc_ep_g.f \
            param.inc                                                    \
            param1.inc                                                   \
            physprop.inc                                                 \
            fldvar.inc                                                   \
            geometry.inc                                                 \
            indices.inc                                                  \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
calc_mu_s.$(OBJ_EXT) : calc_mu_s.f \
            param.inc                                                    \
            param1.inc                                                   \
            parallel.inc                                                 \
            physprop.inc                                                 \
            drag.inc                                                     \
            run.inc                                                      \
            geometry.inc                                                 \
            fldvar.inc                                                   \
            visc_g.inc                                                   \
            visc_s.inc                                                   \
            trace.inc                                                    \
            indices.inc                                                  \
            constant.inc                                                 \
            s_pr1.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            s_pr2.inc                                                   
calc_mw.$(OBJ_EXT) : calc_mw.f \
            param.inc                                                    \
            param1.inc                                                   \
            toleranc.inc                                                
calc_quantities.$(OBJ_EXT) : calc_quantities.f \
            param.inc                                                    \
            param1.inc                                                   \
            post3d.inc                                                  
calc_ro_g.$(OBJ_EXT) : calc_ro_g.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            indices.inc                                                  \
            physprop.inc                                                 \
            geometry.inc                                                 \
            function.inc                                                
calc_vol.$(OBJ_EXT) : calc_vol.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                 \
            indices.inc                                                  \
            fldvar.inc                                                   \
            physprop.inc                                                 \
            post3d.inc                                                   \
            function.inc                                                
check_data_03.$(OBJ_EXT) : check_data_03.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                 \
            funits.inc                                                  
check_one_axis.$(OBJ_EXT) : check_one_axis.f \
            param.inc                                                    \
            param1.inc                                                   \
            funits.inc                                                  
compare.$(OBJ_EXT) : compare.f \
            param.inc                                                    \
            param1.inc                                                   \
            toleranc.inc                                                
eosg.$(OBJ_EXT) : eosg.f \
            param.inc                                                    \
            param1.inc                                                   \
            constant.inc                                                 \
            physprop.inc                                                 \
            scales.inc                                                   \
            sc_p_g1.inc                                                  \
            sc_p_g2.inc                                                 
error_routine.$(OBJ_EXT) : error_routine.f \
            funits.inc                                                  
examine_data.$(OBJ_EXT) : examine_data.f \
            param.inc                                                    \
            param1.inc                                                   \
            constant.inc                                                 \
            physprop.inc                                                 \
            fldvar.inc                                                   \
            indices.inc                                                  \
            run.inc                                                      \
            geometry.inc                                                 \
            post3d.inc                                                   \
            xforms.inc                                                   \
            function.inc                                                
f_init_data.$(OBJ_EXT) : f_init_data.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                 \
            physprop.inc                                                 \
            constant.inc                                                 \
            xforms.inc                                                  
file_handle.$(OBJ_EXT) : file_handle.f \
            machine.inc                                                  \
            param.inc                                                    \
            param1.inc                                                   \
            post3d.inc                                                   \
            funits.inc                                                   \
            geometry.inc                                                 \
            fldvar.inc                                                   \
            physprop.inc                                                 \
            indices.inc                                                  \
            function.inc                                                 \
            run.inc                                                     
finit.$(OBJ_EXT) : finit.f \
            param.inc                                                    \
            param1.inc                                                   \
            run.inc                                                      \
            post3d.inc                                                   \
            geometry.inc                                                 \
            indices.inc                                                  \
            fldvar.inc                                                   \
            physprop.inc                                                 \
            constant.inc                                                 \
            funits.inc                                                   \
            xforms.inc                                                  
flow_gx.$(OBJ_EXT) : flow_gx.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            indices.inc                                                  \
            physprop.inc                                                 \
            geometry.inc                                                 \
            function.inc                                                
flow_gy.$(OBJ_EXT) : flow_gy.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            indices.inc                                                  \
            physprop.inc                                                 \
            geometry.inc                                                 \
            function.inc                                                
flow_gz.$(OBJ_EXT) : flow_gz.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            indices.inc                                                  \
            physprop.inc                                                 \
            geometry.inc                                                 \
            function.inc                                                
flow_sx.$(OBJ_EXT) : flow_sx.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            indices.inc                                                  \
            physprop.inc                                                 \
            geometry.inc                                                 \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
flow_sy.$(OBJ_EXT) : flow_sy.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            indices.inc                                                  \
            physprop.inc                                                 \
            geometry.inc                                                 \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
flow_sz.$(OBJ_EXT) : flow_sz.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            indices.inc                                                  \
            physprop.inc                                                 \
            geometry.inc                                                 \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
g_0.$(OBJ_EXT) : g_0.f \
            param.inc                                                    \
            param1.inc                                                   \
            physprop.inc                                                 \
            fldvar.inc                                                   \
            geometry.inc                                                 \
            indices.inc                                                  \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
gas_flux.$(OBJ_EXT) : gas_flux.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            run.inc                                                      \
            geometry.inc                                                 \
            indices.inc                                                  \
            post3d.inc                                                   \
            physprop.inc                                                 \
            xforms.inc                                                   \
            function.inc                                                
get_file_name.$(OBJ_EXT) : get_file_name.f
get_file_status.$(OBJ_EXT) : get_file_status.f
get_index.$(OBJ_EXT) : get_index.f
get_location.$(OBJ_EXT) : get_location.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                 \
            post3d.inc                                                  
get_mu_s.$(OBJ_EXT) : get_mu_s.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            geometry.inc                                                 \
            physprop.inc                                                 \
            indices.inc                                                  \
            visc_s.inc                                                   \
            run.inc                                                      \
            constant.inc                                                 \
            funits.inc                                                   \
            xforms.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
get_same_time.$(OBJ_EXT) : get_same_time.f \
            param.inc                                                    \
            param1.inc                                                  
get_selection.$(OBJ_EXT) : get_selection.f
get_substr.$(OBJ_EXT) : get_substr.f
granular_qty.$(OBJ_EXT) : granular_qty.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            geometry.inc                                                 \
            physprop.inc                                                 \
            indices.inc                                                  \
            visc_s.inc                                                   \
            run.inc                                                      \
            constant.inc                                                 \
            funits.inc                                                   \
            xforms.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
header_main.$(OBJ_EXT) : header_main.f \
            param.inc                                                    \
            param1.inc                                                   \
            post3d.inc                                                  
ik_avg.$(OBJ_EXT) : ik_avg.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            physprop.inc                                                 \
            indices.inc                                                  \
            geometry.inc                                                 \
            function.inc                                                
ik_avg_out.$(OBJ_EXT) : ik_avg_out.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            indices.inc                                                  \
            run.inc                                                      \
            physprop.inc                                                 \
            post3d.inc                                                  
init_namelist.$(OBJ_EXT) : init_namelist.f \
            param.inc                                                    \
            param1.inc                                                   \
            run.inc                                                      \
            output.inc                                                   \
            physprop.inc                                                 \
            geometry.inc                                                 \
            ic.inc                                                       \
            bc.inc                                                       \
            fldvar.inc                                                   \
            constant.inc                                                 \
            indices.inc                                                  \
            is.inc                                                       \
            toleranc.inc                                                 \
            scales.inc                                                   \
            ur_facs.inc                                                  \
            leqsol.inc                                                   \
            residual.inc                                                 \
            rxns.inc                                                     \
            namelist.inc                                                
interp_res.$(OBJ_EXT) : interp_res.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                 \
            indices.inc                                                  \
            physprop.inc                                                 \
            fldvar.inc                                                   \
            post3d.inc                                                   \
            run.inc                                                      \
            funits.inc                                                   \
            xforms.inc                                                   \
            function.inc                                                
line_too_big.$(OBJ_EXT) : line_too_big.f
main_f.$(OBJ_EXT) : main_f.f \
            xforms.inc       
machine.$(OBJ_EXT) : machine.f param.inc                                      
make_upper_case.$(OBJ_EXT) : make_upper_case.f
out_from_res.$(OBJ_EXT) : out_from_res.f \
            param.inc                                                    \
            param1.inc                                                   \
            run.inc                                                      \
            funits.inc                                                   \
            xforms.inc                                                  
out_from_spx.$(OBJ_EXT) : out_from_spx.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            geometry.inc                                                 \
            indices.inc                                                  \
            run.inc                                                      \
            funits.inc                                                   \
            post3d.inc                                                   \
            xforms.inc                                                  
out_spec_time.$(OBJ_EXT) : out_spec_time.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            run.inc                                                      \
            physprop.inc                                                 \
            indices.inc                                                  \
            geometry.inc                                                 \
            function.inc                                                
out_time.$(OBJ_EXT) : out_time.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            run.inc                                                      \
            physprop.inc                                                 \
            indices.inc                                                  \
            geometry.inc                                                 \
            function.inc                                                
parse_line.$(OBJ_EXT) : parse_line.f \
            param.inc                                                    \
            param1.inc                                                   \
            parse.inc                                                   
parse_rxn.$(OBJ_EXT) : parse_rxn.f \
            param.inc                                                    \
            param1.inc                                                   \
            parse.inc                                                    \
            rxns.inc                                                    
print_out.$(OBJ_EXT) : print_out.f \
            param.inc                                                    \
            param1.inc                                                   \
            post3d.inc                                                   \
            xforms.inc                                                  
read_namelist.$(OBJ_EXT) : read_namelist.f \
            param.inc                                                    \
            param1.inc                                                   \
            run.inc                                                      \
            output.inc                                                   \
            physprop.inc                                                 \
            geometry.inc                                                 \
            ic.inc                                                       \
            is.inc                                                       \
            bc.inc                                                       \
            fldvar.inc                                                   \
            constant.inc                                                 \
            indices.inc                                                  \
            toleranc.inc                                                 \
            funits.inc                                                   \
            scales.inc                                                   \
            ur_facs.inc                                                  \
            leqsol.inc                                                   \
            usrnlst.inc                                                  \
            residual.inc                                                 \
            rxns.inc                                                     \
            namelist.inc                                                
read_res0.$(OBJ_EXT) : read_res0.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                 \
            physprop.inc                                                 \
            run.inc                                                      \
            ic.inc                                                       \
            bc.inc                                                       \
            is.inc                                                       \
            constant.inc                                                 \
            funits.inc                                                   \
            output.inc                                                   \
            scales.inc                                                   \
            ur_facs.inc                                                  \
            toleranc.inc                                                 \
            leqsol.inc                                                  
read_res1.$(OBJ_EXT) : read_res1.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            geometry.inc                                                 \
            physprop.inc                                                 \
            run.inc                                                      \
            funits.inc                                                  
read_spx0.$(OBJ_EXT) : read_spx0.f \
            param.inc                                                    \
            param1.inc                                                   \
            run.inc                                                      \
            funits.inc                                                   \
            post3d.inc                                                  
read_spx1.$(OBJ_EXT) : read_spx1.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            geometry.inc                                                 \
            physprop.inc                                                 \
            run.inc                                                      \
            funits.inc                                                   \
            post3d.inc                                                  
remove_comment.$(OBJ_EXT) : remove_comment.f
res_from_spx.$(OBJ_EXT) : res_from_spx.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                 \
            indices.inc                                                  \
            run.inc                                                      \
            funits.inc                                                   \
            post3d.inc                                                   \
            physprop.inc                                                 \
            fldvar.inc                                                  
seek_comment.$(OBJ_EXT) : seek_comment.f
seek_end.$(OBJ_EXT) : seek_end.f
seek_time.$(OBJ_EXT) : seek_time.f \
            param.inc                                                    \
            param1.inc                                                   \
            post3d.inc                                                   \
            funits.inc                                                  
select_spx_rec.$(OBJ_EXT) : select_spx_rec.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                 \
            indices.inc                                                  \
            run.inc                                                      \
            machine.inc                                                  \
            funits.inc                                                   \
            post3d.inc                                                   \
            physprop.inc                                                 \
            fldvar.inc                                                   \
            xforms.inc                                                  
set_constants.$(OBJ_EXT) : set_constants.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            visc_s.inc                                                   \
            energy.inc                                                   \
            geometry.inc                                                 \
            indices.inc                                                  \
            physprop.inc                                                 \
            constant.inc                                                 \
            run.inc                                                      \
            funits.inc                                                   \
            drag.inc                                                    
set_dollar.$(OBJ_EXT) : set_dollar.f
set_geometry.$(OBJ_EXT) : set_geometry.f \
            param.inc                                                    \
            param1.inc                                                   \
            run.inc                                                      \
            geometry.inc                                                
set_increments.$(OBJ_EXT) : set_increments.f \
            param.inc                                                    \
            param1.inc                                                   \
            indices.inc                                                  \
            geometry.inc                                                 \
            physprop.inc                                                 \
            fldvar.inc                                                   \
            funits.inc                                                   \
            function.inc                                                
set_index1.$(OBJ_EXT) : set_index1.f \
            param.inc                                                    \
            param1.inc                                                   \
            physprop.inc                                                 \
            fldvar.inc                                                   \
            geometry.inc                                                 \
            constant.inc                                                 \
            indices.inc                                                 
set_index1a.$(OBJ_EXT) : set_index1a.f \
            param.inc                                                    \
            param1.inc                                                   \
            physprop.inc                                                 \
            geometry.inc                                                 \
            constant.inc                                                 \
            fldvar.inc                                                   \
            indices.inc                                                  \
            function.inc                                                
set_max2.$(OBJ_EXT) : set_max2.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                
set_read_spx.$(OBJ_EXT) : set_read_spx.f
shift_dxyz.$(OBJ_EXT) : shift_dxyz.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                
sol_flux.$(OBJ_EXT) : sol_flux.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            run.inc                                                      \
            geometry.inc                                                 \
            indices.inc                                                  \
            post3d.inc                                                   \
            physprop.inc                                                 \
            xforms.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
strcmp.$(OBJ_EXT) : strcmp.f
streqs.$(OBJ_EXT) : streqs.f
time_avg.$(OBJ_EXT) : time_avg.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            geometry.inc                                                
usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f \
            usrnlst.inc                                                 
usr_post.$(OBJ_EXT) : usr_post.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            run.inc                                                      \
            geometry.inc                                                 \
            indices.inc                                                  \
            post3d.inc                                                   \
            physprop.inc                                                 \
            function.inc                                                
usr_write_out1.$(OBJ_EXT) : usr_write_out1.f
write_out1.$(OBJ_EXT) : write_out1.f \
            param.inc                                                    \
            param1.inc                                                   \
            physprop.inc                                                 \
            fldvar.inc                                                   \
            run.inc                                                      \
            funits.inc                                                  
write_res0.$(OBJ_EXT) : write_res0.f \
            param.inc                                                    \
            param1.inc                                                   \
            geometry.inc                                                 \
            physprop.inc                                                 \
            run.inc                                                      \
            ic.inc                                                       \
            is.inc                                                       \
            bc.inc                                                       \
            constant.inc                                                 \
            funits.inc                                                   \
            output.inc                                                   \
            scales.inc                                                   \
            ur_facs.inc                                                  \
            leqsol.inc                                                   \
            toleranc.inc                                                
write_res1.$(OBJ_EXT) : write_res1.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            geometry.inc                                                 \
            physprop.inc                                                 \
            run.inc                                                      \
            funits.inc                                                   \
            output.inc                                                  
write_spx0.$(OBJ_EXT) : write_spx0.f \
            param.inc                                                    \
            param1.inc                                                   \
            run.inc                                                      \
            funits.inc                                                  
write_spx1.$(OBJ_EXT) : write_spx1.f \
            param.inc                                                    \
            param1.inc                                                   \
            fldvar.inc                                                   \
            geometry.inc                                                 \
            physprop.inc                                                 \
            run.inc                                                      \
            funits.inc                                                   \
            output.inc                                                  
