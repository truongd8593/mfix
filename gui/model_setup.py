# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division

from constants import *

class ModelSetup(object):
    #    Model Setup Task Pane Window: Select MFIX solver and other conservation equations

    def init_model_setup(self):
        self.subgrid_model = 0
        ui = self.ui.model_setup

        ui.combobox_solver.activated.connect(self.set_solver)
        ui.checkbox_disable_fluid_solver.clicked.connect(self.disable_fluid_solver)

        ui.checkbox_keyword_energy_eq.clicked.connect(self.enable_energy_eq)
        ui.checkbox_enable_turbulence.clicked.connect(self.enable_turbulence)
        ui.combobox_turbulence_model.activated.connect(self.set_turbulence_model)
        ui.combobox_subgrid_model.activated.connect(self.set_subgrid_model)


    def setup_model_setup(self):
        self.disable_fluid_solver(self.fluid_solver_disabled)
        self.set_subgrid_model(0)


    def reset_model_setup(self):
        self.subgrid_model = 0


    def set_solver(self, solver):
        #    Select MFIX Solver:
        # Available selections:
        #  Single phase
        #    Selection disables 'Solids' task pane menu
        #    Selection disables 'Continuum Solids Model' task pane menu
        #    Selection disables 'Discrete Element Model' task pane menu
        #    Selection disables 'Particle-in-Cell' task pane menu
        #  MFIX-TFM
        #    Selection enables 'Solids' task pane menu
        #    Selection enables 'Continuum Solids Model' task pane menu
        #  MFIX-DEM
        #    Selection enables 'Solids' task pane menu
        #    Selection enables 'Discrete Element Model' task pane menu
        #  MFIX-PIC
        #    Selection enables 'Solids' task pane menu
        #    Selection enables 'Particle-in-Cell' task pane menu
        #  MFIX-Hybrid
        #    Selection enables 'Solids' task pane menu
        #    Selection enables 'Continuum Solids Model' task pane menu
        #    Selection enables 'Discrete Element Model' task pane menu

        self.project.solver = solver
        if solver is None: #
            return

        ui = self.ui.model_setup
        cb = ui.combobox_solver
        if cb.currentIndex != solver:
            cb.setCurrentIndex(solver)

        solver_name = {SINGLE:"MFIX Single-Phase",
                       TFM:"MFIX-TFM",
                       DEM:"MFIX-DEM",
                       PIC:"MFIX-PIC",
                       HYBRID:"MFIX-Hybrid"}.get(solver, "MFIX")

        self.print_internal("Solver: %s" % solver_name)
        self.solver_name = solver_name

        enabled = (solver != SINGLE)
        self.find_navigation_tree_item("Solids").setDisabled(not enabled)
        if enabled:
            self.solids_update_tabs()

        # Options which require TFM, DEM, or PIC
        enabled = solver in (TFM, DEM, PIC)

        # TFM only
        enabled = (solver == TFM)
        ui.combobox_subgrid_model.setEnabled(enabled)
        ui.label_subgrid_model.setEnabled(enabled)
        ui.groupbox_subgrid_params.setEnabled(enabled and
                                             self.subgrid_model > 0)

        # Solids Model selection tied to Solver
        # FIXME XXX What to do about solids that are already defined?

        valid_models = (("DEM",) if solver==DEM
                        else ("TFM",) if solver==TFM
                        else ("PIC",) if solver==PIC
                        else ("TFM", "DEM"))

        for (i,(k,v)) in enumerate(self.solids.items(), 1):
            model = v.get('model')
            if model not in valid_models:
                model = valid_models[0]
                self.update_keyword('solids_model', model, args=[i])
                v['model'] = model
        self.update_solids_table() # some of these settings are dependent on solver
        self.setup_combobox_solids_model()
        #self.update_solids_detail_pane()
        self.update_window_title()
        self.update_nav_tree()


    def disable_fluid_solver(self, disabled):
        # Option to disable the fluid phase
        # Disables the Fluid task pane menu
        # Sets keyword RO_G0 to 0.0
        self.fluid_solver_disabled = disabled
        m = self.ui.model_setup
        enabled = not disabled
        item = self.find_navigation_tree_item("Fluid")
        item.setDisabled(disabled)
        if disabled:
            self.enable_turbulence(False)
            self.saved_ro_g0 = self.project.get_value('ro_g0')
            self.update_keyword('ro_g0', 0) # issues/124
            if self.saved_ro_g0 is not None:
                self.ui.fluid.lineedit_keyword_ro_g0.setText(str(self.saved_ro_g0))
        else:
            # would be nice to restore turbulence vals
            val = self.saved_ro_g0 if self.fluid_density_model == CONSTANT else None
            self.update_keyword('ro_g0', val)

        m.checkbox_enable_turbulence.setEnabled(enabled)
        m.combobox_turbulence_model.setEnabled(enabled and
                        m.checkbox_enable_turbulence.isChecked())

        self.update_nav_tree()
        # TODO update nscalar


    def enable_energy_eq(self, enabled):
        # Additional callback on top of automatic keyword update,
        # since this has to change availability of several other GUI items
        ui = self.ui.model_setup
        ui.checkbox_keyword_energy_eq.setChecked(enabled)

        # It might not be necessary to do all this - will the fluid or
        # solid panes get updated before we display them?
        f = self.ui.fluid
        for item in (f.label_fluid_specific_heat_model,
                     f.combobox_fluid_specific_heat_model,
                     f.label_fluid_conductivity_model,
                     f.combobox_fluid_conductivity_model,
                     # more ?
                     ):
            item.setEnabled(enabled)

        # c_pg0 == specific heat for fluid phase
        lineedit = f.lineedit_keyword_c_pg0
        label = f.label_c_pg0_units
        for item in (lineedit, label):
            item.setEnabled(enabled and (self.fluid_specific_heat_model == CONSTANT))

        # k_g0 == thermal conductivity fluid phase
        lineedit = f.lineedit_keyword_k_g0
        label = f.label_k_g0_units
        for item in (lineedit, label):
            item.setEnabled(enabled and (self.fluid_conductivity_model == CONSTANT))


    def set_subgrid_model(self, index):
        ui = self.ui.model_setup
        self.subgrid_model = index
        ui.groupbox_subgrid_params.setEnabled(index > 0)


    def enable_turbulence(self, enabled):
        ui = self.ui.model_setup
        cb = ui.checkbox_enable_turbulence
        if enabled != cb.isChecked():
            cb.setChecked(enabled)
        for item in (ui.label_turbulence_model,
                     ui.combobox_turbulence_model,
                     ui.label_mu_gmax,
                     ui.lineedit_keyword_mu_gmax,
                     ui.label_mu_gmax_units):
            item.setEnabled(enabled)
        if not enabled:
            self.unset_keyword('turbulence_model')
            self.unset_keyword('mu_gmax')
        else:
            self.set_turbulence_model(ui.combobox_turbulence_model.currentIndex())


    def set_turbulence_model(self, val):
        ui = self.ui.model_setup
        cb = ui.combobox_turbulence_model
        if cb.currentIndex() != val:
            cb.setCurrentIndex(val)
        self.update_keyword('turbulence_model',
                            ['MIXING_LENGTH', 'K_EPSILON'][val])
        val = ui.lineedit_keyword_mu_gmax.value
        if val=='' or val is None:
            val = '1.e+03'
        self.update_keyword('mu_gmax', val)



    #    Option to enable thermal energy equations
    # This keyword should always be specified in the input deck
    # Sets keyword ENERGY_EQ
    # DEFAULT value of .FALSE.

    #    Option to enable turbulence
    # Selection available if fluid phase is enabled
    # Available selections:
    #  None; [DEFAULT]
    #  Mixing Length:
    #    Selection always available
    #    Sets keyword TURBULENCE_MODEL to MIXING_LENGTH
    #    Requires IC_L_SCALE for all IC regions
    #  K-Epsilon
    #    Selection always available
    #    Sets keyword TURBULENCE_MODEL to K_EPSILON
    #    Requires IC_K_TURB_G for all IC regions
    #    Requires IC_E_TURB_G for all IC regions
    #    Requires BC_K_TURB_G for inflow (MI and PI) BC regions
    #    Requires BC_E_TURB_G for inflow (MI and PI) BC regions

    #Specify maximum fluid viscosity (not shown in mockup)
    # Selection available if TURBULENCE_MODEL =/ 'NONE'
    # Sets keyword MU_GMAX
    # DEFAULT value of 1.0e3 (Pa.s)

    #Specify drag model
    # Selection requires TFM, DEM, or PIC solver
    # Sets keyword DRAG_TYPE
    # Available selections:
    #  SYAM_OBRIEN (DEFAULT)
    #    Specify model parameter: DRAG_C1
    #      DEFAULT value of 0.8
    #    Specify model parameter: DRAG_D1
    #      DEFAULT value of 2.65
    #  BVK
    #  GIDASPOW
    #  GIDASPOW_BLEND
    #  GIDASPOW_PCF
    #  GIDASPOW_BLEND_PCF
    #  HYS
    #    Specify model parameter LAM_HYS
    #      DEFAULT value of 1.0e-6 (meters)
    #  KOCH_HILL
    #  KOCH_HILL_PCF
    #  WEN_YU
    #  WEN_YU_PCF
    #  USER_DRAG

    #Specify heat transfer correlation (requires TFM, DEM, or PIC solver)
    #*This option may be premature as MFIX is limited in heat HTCs.

    #Specify momentum equation formulation; Select Model A, Model B; Jackson, Ishi

    #Select sub-grid model:
    # Selection requirements:
    #  Only available with MFIX-TFM solver
    #  DRAG_TYPE="WEN_YU"
    #  KT_TYPE="ALGEBRAIC"
    #  TURBULENCE_MODEL /= K_EPSILON
    #  BLENDING_STRESS = NONE
    #  FRICTION_MODEL /= SRIVASTAVA
    #  (There are more restrictions…)
    # Sets keyword SUBGRID_TYPE
    # Available selections
    #  NONE (DEFAULT)
    #  IGCI
    #  MILIO

    #Specify sub-grid model filter size ratio:
    # Specification requires SUBGRID_TYPE =/ NONE
    # Sets keyword FILTER_SIZE_RATIO
    # DEFAULT value of 2.0

    #Enable sub-grid wall correction model:
    # Specification requires SUBGRID_TYPE =/ NONE
    # Sets keyword SUBGRID_WALL
    # DEFAULT value of FALSE
