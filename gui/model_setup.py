# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division

class ModelSetup:
    #    Model Setup Task Pane Window: Select MFIX solver and other conservation equations

    def init_model_setup(self):
        ui = self.ui.model_setup
        ui.combobox_solver.activated.connect(self.set_solver)

        ui = self.ui
        model = ui.model_setup
        combobox = model.combobox_solver
        # activated: Only on user action, avoid recursive calls in set_solver
        combobox.activated.connect(self.set_solver)

        checkbox = model.checkbox_disable_fluid_solver
        checkbox.clicked.connect(self.disable_fluid_solver)
        self.disable_fluid_solver(self.fluid_solver_disabled)

        checkbox = model.checkbox_keyword_energy_eq
        checkbox.clicked.connect(self.enable_energy_eq)

        checkbox = model.checkbox_enable_turbulence
        checkbox.clicked.connect(self.enable_turbulence)

        combobox = model.combobox_turbulence_model
        combobox.currentIndexChanged.connect(self.set_turbulence_model)

        combobox = model.combobox_subgrid_model
        combobox.currentIndexChanged.connect(self.set_subgrid_model)
        self.set_subgrid_model(0)




    def setup_model_setup(self):
        pass


    # Top-level "Model"
    def set_solver(self, solver):
        """handler for "Solver" combobox in Model pane"""
        self.project.solver = solver
        if solver is None: #
            return

        ui = self.ui
        m = ui.model_setup
        cb = m.combobox_solver
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
        #interphase = m.groupbox_interphase
        #interphase.setEnabled(enabled)

        # TFM only
        enabled = (solver == TFM)
        m.combobox_subgrid_model.setEnabled(enabled)
        m.label_subgrid_model.setEnabled(enabled)
        m.groupbox_subgrid_params.setEnabled(enabled and
                                                       self.subgrid_model > 0)

        enabled = (self.fluid_nscalar_eq > 0)
        ui.fluid.checkbox_enable_scalar_eq.setChecked(self.fluid_nscalar_eq>0)
        ui.fluid.spinbox_nscalar_eq.setEnabled(enabled)

        # Equiv for solids is done in update_solids_detail_pane

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

        self.ui.model_setup.checkbox_keyword_energy_eq.setChecked(enabled)

        # It might not be necessary to do all this - will the fluid or
        # solid panes get updated before we display them?
        ui = self.ui
        f = ui.fluid
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
        self.subgrid_model = index
        groupbox_subgrid_params = self.ui.model_setup.groupbox_subgrid_params
        groupbox_subgrid_params.setEnabled(index > 0)


    def enable_turbulence(self, enabled):
        m = self.ui.model_setup
        if enabled != m.checkbox_enable_turbulence.isChecked():
            m.checkbox_enable_turbulence.setChecked(enabled)
        for item in (m.label_turbulence_model,
                     m.combobox_turbulence_model,
                     m.label_mu_gmax,
                     m.lineedit_keyword_mu_gmax,
                     m.label_mu_gmax_units):
            item.setEnabled(enabled)
        if not enabled:
            self.unset_keyword('turbulence_model')
            self.unset_keyword('mu_gmax')
        else:
            self.set_turbulence_model(m.combobox_turbulence_model.currentIndex())


    def set_turbulence_model(self, val):
        m = self.ui.model_setup
        cb = m.combobox_turbulence_model
        if cb.currentIndex() != val:
            cb.setCurrentIndex(val)
            return
        self.update_keyword('turbulence_model',
                            ['MIXING_LENGTH', 'K_EPSILON'][val])
        val = m.lineedit_keyword_mu_gmax.value
        if val=='' or val is None:
            val = '1.e+03'
        self.update_keyword('mu_gmax', val)





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
    #    Option to disable the fluid phase
    # Disables the 'Fluid' task pane menu
    # Sets keyword RO_G0 to 0.0
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
    # DEFAULT value of 0.8
    #    Specify model parameter: DRAG_D1
    # DEFAULT value of 2.65
    #  BVK
    #  GIDASPOW
    #  GIDASPOW_BLEND
    #  GIDASPOW_PCF
    #  GIDASPOW_BLEND_PCF
    #  HYS
    #    Specify model parameter LAM_HYS
    # DEFAULT value of 1.0e-6 (meters)
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
