#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division

""" Factors for converting CGS to SI models, etc
    converting CGS to SI models"""

cgs_to_SI = {
#<start generated code>
    'aggregation_eff': 1,
    'alpha_max': 1,
    'asperities': 1e-2, # cm->m,
    'bar_resolution': 1,
    'batch_wallclock': 1,
    'bc_c_scalar': 1e2, # 1/cm -> 1/m
    'bc_c_t_g': 1e2,
    'bc_c_t_s': 1e2,
    'bc_c_theta_m': 1e2,
    'bc_c_x_g': 1e2,
    'bc_c_x_s': 1e2,
    'bc_dt_0': 1,
    'bc_dt_h': 1,
    'bc_dt_l': 1,
    'bc_e_turb_g': 1e-4, # cm^2/s^3 -> m^2/s^3
    'bc_ep_g' : 1,
    'bc_ep_s' : 1,
#   'bc_hw_g': UNDEFINED,
#   'bc_hw_s': UNDEFINED,
#   'bc_hw_scalar': UNDEFINED,
#   'bc_hw_t_g': UNDEFINED,
#   'bc_hw_t_s': UNDEFINED,
#   'bc_hw_theta_m': UNDEFINED,
#   'bc_hw_x_g': UNDEFINED,
#   'bc_hw_x_s': UNDEFINED,
#   'bc_jet_g0': UNDEFINED,
#   'bc_jet_gh': UNDEFINED,
#   'bc_jet_gl': UNDEFINED,
    'bc_k_turb_g': 1e-4, # cm^2/s^3 -> m^2/s^3
#   'bc_massflow_g': UNDEFINED,
#   'bc_massflow_s': UNDEFINED,
#   'bc_p_g': UNDEFINED,
#   'bc_pic_mi_const_statwt': UNDEFINED,
#   'bc_rop_s': UNDEFINED,
#   'bc_scalar': UNDEFINED,
#   'bc_scalarw': UNDEFINED,
#   'bc_t_g': UNDEFINED,
#   'bc_t_s': UNDEFINED,
#   'bc_theta_m': UNDEFINED,
#   'bc_thetaw_m': UNDEFINED,
#   'bc_tw_g': UNDEFINED,
#   'bc_tw_s': UNDEFINED,
#   'bc_u_g': UNDEFINED,
#   'bc_u_s': UNDEFINED,
#   'bc_uw_g': UNDEFINED,
#   'bc_uw_s': UNDEFINED,
#   'bc_v_g': UNDEFINED,
#   'bc_v_s': UNDEFINED,
#   'bc_velmag_g': UNDEFINED,
#   'bc_velmag_s': UNDEFINED,
#   'bc_volflow_g': UNDEFINED,
#   'bc_volflow_s': UNDEFINED,
#   'bc_vw_g': UNDEFINED,
#   'bc_vw_s': UNDEFINED,
#   'bc_w_g': UNDEFINED,
#   'bc_w_s': UNDEFINED,
#   'bc_ww_g': UNDEFINED,
#   'bc_ww_s': UNDEFINED,
#   'bc_x_e': UNDEFINED,
    'bc_x_g' : 1,
    'bc_x_s' : 1,
#   'bc_x_w': UNDEFINED,
    'bc_xw_g' : 1,
    'bc_xw_s' : 1,
#   'bc_y_n': UNDEFINED,
#   'bc_y_s': UNDEFINED,
#   'bc_z_b': UNDEFINED,
#   'bc_z_t': UNDEFINED,
#   'bend_r1': UNDEFINED,
#   'bend_r2': UNDEFINED,
#   'bend_theta1': UNDEFINED,
#   'bend_theta2': UNDEFINED,
#   'breakage_eff': UNDEFINED,
#   'c': UNDEFINED,
#   'c2c_r1': UNDEFINED,
#   'c2c_r2': UNDEFINED,
#   'c2c_y2': UNDEFINED,
#   'c_e': UNDEFINED,
#   'c_f': UNDEFINED,
#   'c_fac': UNDEFINED,
    'c_pg0' : 4184.0,
    'c_ps0' : 4184.0,
#   'cfl_pic': UNDEFINED,
#   'cg_ur_fac': UNDEFINED,
#   'clip_xmax': UNDEFINED,
#   'clip_xmin': UNDEFINED,
#   'clip_ymax': UNDEFINED,
#   'clip_ymin': UNDEFINED,
#   'clip_zmax': UNDEFINED,
#   'clip_zmin': UNDEFINED,
#   'cpx': UNDEFINED,
#   'cpy': UNDEFINED,
#   'cpz': UNDEFINED,
    'd_p0' : 0.01,
#   'delp_x': UNDEFINED,
#   'delp_y': UNDEFINED,
#   'delp_z': UNDEFINED,
#   'des_diffuse_width': UNDEFINED,
#   'des_em': UNDEFINED,
#   'des_en_input': UNDEFINED,
#   'des_en_wall_input': UNDEFINED,
#   'des_et_input': UNDEFINED,
#   'des_et_wall_input': UNDEFINED,
#   'des_etat_fac': UNDEFINED,
#   'des_etat_w_fac': UNDEFINED,
#   'des_interp_width': UNDEFINED,
#   'des_min_cond_dist': UNDEFINED,
    'dif_g0' : 0.0001,
    'dif_s0' : 0.0001,
    'dil_factor_vsd' : 1,
    'dil_inert_x_vsd' : 1,
#   'dquadric': UNDEFINED,
#   'drag_c1': UNDEFINED,
#   'drag_d1': UNDEFINED,
#   'dt': UNDEFINED,
#   'dt_fac': UNDEFINED,
#   'dt_max': UNDEFINED,
#   'dt_min': UNDEFINED,
#   'dx': UNDEFINED,
#   'dy': UNDEFINED,
#   'dz': UNDEFINED,
#   'e_w': UNDEFINED,
    'e_young' : 0.1,
    'e_young_actual' : 0.1,
    'ep_s_max' : 1,
    'ep_star' : 1,
    'eps_f_min' : 1,
#   'erx': UNDEFINED,
#   'ery': UNDEFINED,
#   'erz': UNDEFINED,
    'ew_young' : 0.1,
    'ew_young_actual' : 0.1,
    'fac_dim_max_cut_cell' : 1,
#   'factor_rlm': UNDEFINED,
#   'filter_size_ratio': UNDEFINED,
#   'first_dx': UNDEFINED,
#   'first_dy': UNDEFINED,
#   'first_dz': UNDEFINED,
#   'flpc': UNDEFINED,
#   'flux_g': UNDEFINED,
#   'fric_exp_pic': UNDEFINED,
    'gravity' : 0.01,
    'gravity_x' : 0.01,
    'gravity_y' : 0.01,
    'gravity_z' : 0.01,
#   'half_angle': UNDEFINED,
#   'hamaker_constant': UNDEFINED,
#   'ic_e_turb_g': UNDEFINED,
    'ic_ep_g' : 1,
    'ic_ep_s' : 1,
#   'ic_gama_rg': UNDEFINED,
#   'ic_gama_rs': UNDEFINED,
#   'ic_k_turb_g': UNDEFINED,
#   'ic_l_scale': UNDEFINED,
#   'ic_p_g': UNDEFINED,
#   'ic_p_star': UNDEFINED,
#   'ic_pic_const_statwt': UNDEFINED,
#   'ic_rop_s': UNDEFINED,
#   'ic_scalar': UNDEFINED,
#   'ic_t_g': UNDEFINED,
#   'ic_t_rg': UNDEFINED,
#   'ic_t_rs': UNDEFINED,
#   'ic_t_s': UNDEFINED,
#   'ic_theta_m': UNDEFINED,
#   'ic_u_g': UNDEFINED,
#   'ic_u_s': UNDEFINED,
#   'ic_v_g': UNDEFINED,
#   'ic_v_s': UNDEFINED,
#   'ic_w_g': UNDEFINED,
#   'ic_w_s': UNDEFINED,
#   'ic_x_e': UNDEFINED,
    'ic_x_g' : 1,
    'ic_x_s' : 1,
#   'ic_x_w': UNDEFINED,
#   'ic_y_n': UNDEFINED,
#   'ic_y_s': UNDEFINED,
#   'ic_z_b': UNDEFINED,
#   'ic_z_t': UNDEFINED,
#   'is_pc': UNDEFINED,
#   'is_vel_s': UNDEFINED,
#   'is_x_e': UNDEFINED,
#   'is_x_w': UNDEFINED,
#   'is_y_n': UNDEFINED,
#   'is_y_s': UNDEFINED,
#   'is_z_b': UNDEFINED,
#   'is_z_t': UNDEFINED,
    'k_g0' : 418.4,
    'k_s0' : 418.4,
    'kn' : 0.001,
    'kn_w' : 0.001,
#   'kt_fac': UNDEFINED,
#   'kt_w_fac': UNDEFINED,
#   'lam_hys': UNDEFINED,
#   'lambda_x': UNDEFINED,
#   'lambda_y': UNDEFINED,
#   'lambda_z': UNDEFINED,
#   'last_dx': UNDEFINED,
#   'last_dy': UNDEFINED,
#   'last_dz': UNDEFINED,
#   'leq_tol': UNDEFINED,
#   'max_inlet_vel_fac': UNDEFINED,
#   'mew': UNDEFINED,
#   'mew_w': UNDEFINED,
#   'mppic_coeff_en1': UNDEFINED,
#   'mppic_coeff_en2': UNDEFINED,
#   'mppic_coeff_en_wall': UNDEFINED,
#   'mppic_coeff_et_wall': UNDEFINED,
    'mu_g0' : 0.1,
#   'mu_gmax': UNDEFINED,
#   'mu_s0': UNDEFINED,
    'mw_avg' : 1.0,
    'mw_g' : 1.0,
    'mw_s' : 1.0,
#   'n_x': UNDEFINED,
#   'n_y': UNDEFINED,
#   'n_z': UNDEFINED,
#   'neighbor_search_rad_ratio': UNDEFINED,
#   'norm_g': UNDEFINED,
    'norm_s' : 1,
#   'out_dt': UNDEFINED,
#   'out_msh_value': UNDEFINED,
#   'out_stl_value': UNDEFINED,
#   'p_ref': UNDEFINED,
#   'p_scale': UNDEFINED,
#   'phi': UNDEFINED,
#   'phi_w': UNDEFINED,
#   'phip': UNDEFINED,
#   'phip0': UNDEFINED,
#   'piece_xmax': UNDEFINED,
#   'piece_xmin': UNDEFINED,
#   'piece_ymax': UNDEFINED,
#   'piece_ymin': UNDEFINED,
#   'piece_zmax': UNDEFINED,
#   'piece_zmin': UNDEFINED,
#   'ps_massflow_g': UNDEFINED,
#   'ps_massflow_s': UNDEFINED,
#   'ps_t_g': UNDEFINED,
#   'ps_t_s': UNDEFINED,
#   'ps_u_g': UNDEFINED,
#   'ps_u_s': UNDEFINED,
#   'ps_v_g': UNDEFINED,
#   'ps_v_s': UNDEFINED,
#   'ps_w_g': UNDEFINED,
#   'ps_w_s': UNDEFINED,
#   'ps_x_e': UNDEFINED,
    'ps_x_g' : 1,
    'ps_x_s' : 1,
#   'ps_x_w': UNDEFINED,
#   'ps_y_n': UNDEFINED,
#   'ps_y_s': UNDEFINED,
#   'ps_z_b': UNDEFINED,
#   'ps_z_t': UNDEFINED,
#   'quadric_scale': UNDEFINED,
#   'r_p': UNDEFINED,
#   'radius': UNDEFINED,
#   'reactor1_r1': UNDEFINED,
#   'reactor1_r2': UNDEFINED,
#   'reactor1_rr1': UNDEFINED,
#   'reactor1_rr2': UNDEFINED,
#   'reactor1_theta1': UNDEFINED,
#   'reactor1_theta2': UNDEFINED,
#   'reactor1_y1': UNDEFINED,
#   'reactor1_y2': UNDEFINED,
#   'reactor1_yr1': UNDEFINED,
#   'reactor1_yr2': UNDEFINED,
#   'report_mass_balance_dt': UNDEFINED,
#   'res_backup_dt': UNDEFINED,
#   'res_dt': UNDEFINED,
    'ro_g0' : 1000.0,
    'ro_s0' : 1000.0,
    'ro_xs0' : 1000.0,
#   'scale_msh': UNDEFINED,
#   'scale_stl': UNDEFINED,
#   'segregation_slope_coefficient': UNDEFINED,
    'spx_dt' : 1,
#   'stl_small_angle': UNDEFINED,
#   't_x': UNDEFINED,
#   't_y': UNDEFINED,
#   't_z': UNDEFINED,
#   'term_buffer': UNDEFINED,
#   'theta_x': UNDEFINED,
#   'theta_y': UNDEFINED,
#   'theta_z': UNDEFINED,
#   'time': UNDEFINED,
    'tol_delh' : 1,
#   'tol_diverge': UNDEFINED,
#   'tol_f': UNDEFINED,
#   'tol_msh': UNDEFINED,
#   'tol_poly': UNDEFINED,
#   'tol_resid': UNDEFINED,
#   'tol_resid_k_epsilon': UNDEFINED,
#   'tol_resid_scalar': UNDEFINED,
#   'tol_resid_t': UNDEFINED,
#   'tol_resid_th': UNDEFINED,
#   'tol_resid_x': UNDEFINED,
    'tol_small_area' : 1,
    'tol_small_cell' : 1,
    'tol_snap' : 1,
#   'tol_stl': UNDEFINED,
#   'tol_stl_dp': UNDEFINED,
#   'torus_r1': UNDEFINED,
#   'torus_r2': UNDEFINED,
#   'tstop': UNDEFINED,
#   'tx_msh': UNDEFINED,
#   'tx_stl': UNDEFINED,
#   'ty_msh': UNDEFINED,
#   'ty_stl': UNDEFINED,
#   'tz_msh': UNDEFINED,
#   'tz_stl': UNDEFINED,
#   'ucoil_r1': UNDEFINED,
#   'ucoil_r2': UNDEFINED,
#   'ucoil_y2': UNDEFINED,
#   'ur_f_gs': UNDEFINED,
#   'ur_fac': UNDEFINED,
#   'ur_kth_sml': UNDEFINED,
#   'usr_dt': UNDEFINED,
#   'usr_x_e': UNDEFINED,
#   'usr_x_w': UNDEFINED,
#   'usr_y_n': UNDEFINED,
#   'usr_y_s': UNDEFINED,
#   'usr_z_b': UNDEFINED,
#   'usr_z_t': UNDEFINED,
#   'v_ex': UNDEFINED,
#   'v_poisson': UNDEFINED,
#   'v_poisson_actual': UNDEFINED,
#   'v_sh': UNDEFINED,
#   'vdw_inner_cutoff': UNDEFINED,
#   'vdw_outer_cutoff': UNDEFINED,
#   'vtk_dt': UNDEFINED,
#   'vtk_slice_tol': UNDEFINED,
#   'vtk_x_e': UNDEFINED,
#   'vtk_y_n': UNDEFINED,
#   'vtk_z_t': UNDEFINED,
#   'vw_poisson': UNDEFINED,
#   'vw_poisson_actual': UNDEFINED,
#   'wall_hamaker_constant': UNDEFINED,
#   'wall_vdw_inner_cutoff': UNDEFINED,
#   'wall_vdw_outer_cutoff': UNDEFINED,
    'x_s0' : 1,
#   'xlength': UNDEFINED,
#   'xmin': UNDEFINED,
#   'ylength': UNDEFINED,
#   'zlength': UNDEFINED,
#<end generated code>
}

if (__name__ == '__main__'):
    import os
    import re
    from tools.general import SCRIPT_DIRECTORY
    from tools.namelistparser import buildKeywordDoc

    doc = buildKeywordDoc(os.path.join(SCRIPT_DIRECTORY, os.pardir, 'model'))
    keys = doc.keys()
    keys.sort()

    pat = re.compile(r'\[(.*) in SI.*\]')


    unit_to_SI = { 'm':    1e-2,  # from cm
                   'kg':   1e-3,  # from g
                   'J':    4.184, # from cal
                   'kmol': 1e-3,  # from mol
                   'N':    1e-5,  # from dyne
                   'Pa':   0.1,   # from barye
                   'K':    1,
                   's':    1
    }

    # Try to get units from keyword doc. Add to predefined factors above

    for key in keys[:]:
        if key in cgs_to_SI:
            continue
        if 'gravity' in key:
            factor = 1e-2 # cm -> m
            cgs_to_SI[key] = factor
            continue
        entry = doc[key]
        if not isinstance(entry, dict):
            keys.remove(key)
            del doc[key]
            continue
        dtype = entry['dtype']
        if dtype !=  'DP':
            keys.remove(key)
            continue
        desc = entry['description']
        dlower = desc.lower()

        if 'fraction' in dlower:
            cgs_to_SI[key] = 1
            continue

        match = pat.search(desc)

        if match:
            SI_unit = match.group(1)
            if SI_unit.startswith('(') and SI_unit.endswith(')'):
                SI_unit = SI_unit[1:-1]
            if '/' in SI_unit:
                num, denom = SI_unit.split('/')
            else:
                num, denom = SI_unit, ''
            if num.startswith('(') and num.endswith(')'):
                num = num[1:-1]
            if denom.startswith('(') and denom.endswith(')'):
                denom = denom[1:-1]
            factor = 1.0
            for term in num.split('.'):
                if '^' in term:
                    u, pow = term.split('^')
                    pow = int(pow)
                else:
                    u = term
                    pow = 1
                factor *= unit_to_SI[u]**pow
            if denom:
                for term in denom.split('.'):
                    if '^' in term:
                        u, pow = term.split('^')
                        pow = int(pow)
                    else:
                        u = term
                        pow = 1
                    factor /= unit_to_SI[u]**pow
            cgs_to_SI[key] = round(factor,10)
    infile = open(__file__, 'r')
    outfile = open(__file__ + '.tmp', 'w')
    skip = False
    for line in infile:
        line = line.decode('utf-8')
        if line.startswith('#<start'):
            outfile.write(line)
            for key in keys:
                factor = cgs_to_SI.get(key)
                if factor is None:
                    outfile.write("#   '%s': UNDEFINED,\n" % key)
                else:
                    outfile.write("    '%s' : %s,\n" % (key, factor))
            skip = True
        elif line.startswith('#<end'):
            outfile.write(line)
            skip = False
            continue
        elif not skip:
            outfile.write(line.encode('utf-8'))

    # We could replace the input file, but let's not go crazy
    print("Generated %s" % (__file__ + '.tmp'))



"""
Simulation Units
Simulations can be setup using the International System of Units (SI) or the centimetergram-second system (CGS). Although the majority of units are consistent with the specified systems, there are exceptions. The following table provides the SI and CGS units employed by MFIX for various quantities.


Quantity                      MFIX SI unit              MFIX CGS unit

length, position              meter(m)                  centimeter (cm)
mass                          kilogram (kg)             gram (g)
time                          second (s)                second (s)
thermal temperature           Kelvin (K)                Kelvin (K)
energy†                       Joule (J)                 calorie (cal)
amount of substance‡          kilomole (kmol)           mole (mol)
force                         Newton (1 N = 1 kg·m·s-2) dyne (1 dyn = 1 g·cm·s-2)
pressure                      Pascal (1 Pa = 1 N·m-2)   barye (1 Ba = 1 dyn·cm-2)
dynamic viscosity             Pa·s                      poise (1 P = 1 g·cm-1·s-1)
kinematic viscosity           m2·s-1                    Stokes (1 St = 1 cm2·s-1)
gas constant                  J·K-1·kmol-1              erg·K-1·mol-1
enthalpy                      J                         cal
specific heat                 J·kg-1·K-1                cal·g-1·K-1
thermal conductivity          J·m-1·K-1·s-1             cal·cm-1·K-1·s-1

† The CGS unit for energy is the ergon (1 erg = 1 dyne·cm). This is reflected in MFIX
through the gas constant. However, all thermochemical properties related to the energy
equations are based in calories for CGS units. Entries in the Burcat database are
always specified in terms of calories regardless of the simulation units. MFIX converts
the entries to joules after reading the database when SI units are used.
‡ The SI unit for the amount of a substance is the mole (mol). These units are needed
when specifying reaction rates:
• amount per time per volume for Eulerian model reactions
• amount per time for Lagrangian model reactions
"""
