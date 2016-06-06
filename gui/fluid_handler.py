# methods to deal with fluids, split off from gui.py

from __future__ import print_function, absolute_import, unicode_literals, division

from copy import deepcopy

#import Qt
from qtpy import QtWidgets, PYQT5

#local imports
from constants import *
from tools.general import set_item_noedit, get_selected_row

class FluidHandler(object):
    # Defaults
    def init_fluid_default_models(self):
        self.fluid_density_model = CONSTANT
        self.fluid_viscosity_model = CONSTANT
        self.fluid_mol_weight_model = CONSTANT
        self.fluid_specific_heat_model = CONSTANT
        self.fluid_conductivity_model = AIR
        self.fluid_diffusion_model = AIR

    ## Fluid phase methods
    def enable_fluid_species_eq(self, state):
        ui = self.ui
        fluid = ui.fluid
        for item in (fluid.combobox_fluid_diffusion_model,
                     # more ?
                     ):
            item.setEnabled(state)
        # dif_g0 == diffusion coeff model
        lineedit = fluid.lineedit_keyword_dif_g0
        if state:
            lineedit.setEnabled(self.fluid_diffusion_model == CONSTANT)
        else:
            lineedit.setEnabled(False)

    def enable_fluid_scalar_eq(self, state):
        spinbox = self.ui.fluid.spinbox_fluid_nscalar_eq
        spinbox.setEnabled(state)
        if state:
            val = spinbox.value()
            self.set_fluid_nscalar_eq(val)
        else:
            # Don't call set_fluid_nscalar_eq(0) b/c that will clobber spinbox
            prev_nscalar = self.fluid_nscalar_eq + self.solid_nscalar_eq
            self.fluid_nscalar_eq = 0
            self.update_scalar_equations(prev_nscalar)

    def set_fluid_nscalar_eq(self, value):
        # This *sums into* nscalar - not a simple keyword
        prev_nscalar = self.fluid_nscalar_eq + self.solid_nscalar_eq
        self.fluid_nscalar_eq = value
        spinbox = self.ui.fluid.spinbox_fluid_nscalar_eq
        if value != spinbox.value():
            spinbox.setValue(value)
        self.update_scalar_equations(prev_nscalar)

    def init_fluid_handler(self):
        ui = self.ui
        ui.fluid.lineedit_fluid_phase_name.default_value = "Fluid"
        self.saved_fluid_species = None
        self.init_fluid_default_models()
        # Handle a number of cases which are essentially the same
        # see 'set_fluid_mol_weight_model' below to help understand this
        def make_fluid_model_setter(self, name, key):
            def setter(model):
                setattr(self, name, model) # self.fluid_<name>_model = model
                combobox = getattr(self.ui.fluid, 'combobox_' + name)
                prev_model = combobox.currentIndex()
                if model != prev_model:
                    combobox.setCurrentIndex(model)

                # Enable lineedit for constant model
                key_g0 = 'c_pg0' if key=='c_p' else key + '_g0'
                key_usr = 'usr_cpg' if key=='c_p' else 'usr_' + key + 'g'
                lineedit = getattr(self.ui.fluid,
                                   'lineedit_keyword_%s' % key_g0)
                lineedit.setEnabled(model==CONSTANT)

                if model == CONSTANT:
                    value = lineedit.value # Possibly re-enabled gui item
                    if self.project.get_value(key_g0) != value:
                        self.set_keyword(key_g0, value) # Restore keyword value
                elif model == UDF:
                    self.unset_keyword(key_g0)
                    self.set_keyword(key_usr, True)
                else: # Ideal gas law, Sutherland, etc
                    self.unset_keyword(key_g0)
                    self.unset_keyword(key_usr)
                    # anything else to do in this case? validation?
            return setter


        # Create setters for the cases which are similar (mol. wt. handled separately)
        for (name, key) in (
                ('density', 'ro'),
                ('viscosity', 'mu'),
                ('specific_heat', 'c_p'),
                ('conductivity', 'k'),
                ('diffusion', 'dif')):
            model_name = 'fluid_%s_model' % name
            setattr(self, 'set_'+model_name, make_fluid_model_setter(self, model_name, key))

            # Set the combobox default value
            combobox = getattr(self.ui.fluid, 'combobox_'+model_name)
            combobox.default_value = getattr(self, model_name)
            #print(model_name, combobox.default_value)

        combobox = self.ui.fluid.combobox_fluid_mol_weight_model
        combobox.default_value = self.fluid_mol_weight_model

        # more stuff moved from gui.__init__
        checkbox = ui.fluid.checkbox_keyword_species_eq_args_0
        checkbox.clicked.connect(self.enable_fluid_species_eq)

        ui.fluid.lineedit_fluid_phase_name.editingFinished.connect(
            self.handle_fluid_phase_name)
        ui.fluid.checkbox_enable_fluid_scalar_eq.clicked.connect(
            self.enable_fluid_scalar_eq)
        ui.fluid.spinbox_fluid_nscalar_eq.valueChanged.connect(
            self.set_fluid_nscalar_eq)

        # Fluid phase models
        # Density
        for name in ('density', 'viscosity', 'specific_heat', 'mol_weight',
                     'conductivity', 'diffusion'):
            combobox = getattr(ui.fluid, 'combobox_fluid_%s_model' % name)
            setter = getattr(self,'set_fluid_%s_model' % name)
            combobox.currentIndexChanged.connect(setter)

        # Fluid species
        f = ui.fluid
        tb = f.toolbutton_fluid_species_add
        tb.clicked.connect(self.fluid_species_add)
        tb = f.toolbutton_fluid_species_copy # misnomer
        tb.clicked.connect(self.fluid_species_edit)
        tb.setEnabled(False)
        tb = f.toolbutton_fluid_species_delete
        tb.setEnabled(False)
        tb.clicked.connect(self.fluid_species_delete)
        tw = f.tablewidget_fluid_species
        tw.itemSelectionChanged.connect(self.handle_fluid_species_selection)


    # molecular wt model only has 2 choices, and the key names don't
    # follow the same pattern, so create its setter specially
    def set_fluid_mol_weight_model(self, model):
        self.fluid_mol_weight_model = model
        combobox = self.ui.fluid.combobox_fluid_mol_weight_model
        prev_model = combobox.currentIndex()
        if model != prev_model:
            combobox.setCurrentIndex(model)
        # Enable lineedit for constant mol_weight model
        lineedit = self.ui.fluid.lineedit_keyword_mw_avg
        lineedit.setEnabled(model==CONSTANT)
        if model == CONSTANT:
            value = lineedit.value # Possibly re-enabled gui item
            if self.project.get_value("mw_avg") != value:
                self.set_keyword("mw_avg", value) # Restore keyword value
        else: # Mixture
            # TODO: validate, require mw for all component species
            self.unset_keyword("mw_avg")


    def handle_fluid_phase_name(self): # editingFinished signal does not include value
        value = self.ui.fluid.linedit_fluid_phase_name.text()
        self.set_fluid_phase_name(value)

    def set_fluid_phase_name(self, value):
        if value != self.ui.fluid.lineedit_fluid_phase_name.text():
            self.ui.fluid.lineedit_fluid_phase_name.setText(value)
        self.project.mfix_gui_comments['fluid_phase_name'] = value
        self.set_unsaved_flag()


    # --- fluid species methods ---
    def fluid_species_revert(self):
        if self.saved_fluid_species is None:
            return
        self.fluid_species = self.saved_fluid_species
        self.species_popup.defined_species = deepcopy(self.fluid_species)
        self.update_fluid_species_table()

    def fluid_species_save(self):
        self.fluid_species = deepcopy(self.species_popup.defined_species)
        self.update_fluid_species_table()

    def update_fluid_species_table(self):
        """Update table in fluids pane.  Also set nmax_g, species_g and species_alias_g keywords,
        which are not tied to a single widget"""

        hv = QtWidgets.QHeaderView
        table = self.ui.fluid.tablewidget_fluid_species
        if PYQT5:
            resize = table.horizontalHeader().setSectionResizeMode
        else:
            resize = table.horizontalHeader().setResizeMode
        for n in range(5):
            resize(n, hv.ResizeToContents if n>0
                   else hv.Stretch)

        table.clearContents()
        if self.fluid_species is None:
            return
        nrows = len(self.fluid_species)
        table.setRowCount(nrows)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        old_nmax_g = self.project.get_value('nmax_g')
        nmax_g = len(self.fluid_species)
        if nmax_g > 0:
            self.update_keyword('nmax_g', nmax_g)
        else:
            self.unset_keyword('nmax_g')
        for (row, (species,data)) in enumerate(self.fluid_species.items()):
            for (col, key) in enumerate(('alias', 'phase', 'mol_weight',
                                        'h_f', 'source')):
                alias = data.get('alias', species) # default to species if no alias
                data['alias'] = alias # for make_item
                table.setItem(row, col, make_item(data.get(key)))
                self.update_keyword('species_g', species, args=row+1)
                self.update_keyword('species_alias_g', alias, args=row+1)
                # We're avoiding mw_g in favor of the settings in THERMO DATA
                #self.update_keyword('mw_g', data['mol_weight'], args=row+1)#

        # Clear any keywords with indices above nmax_g
        if old_nmax_g is None:
            old_nmax_g = 0
        for i in range(nmax_g+1, old_nmax_g+1):
            self.unset_keyword('species_g', i)
            self.unset_keyword('species_alias_g', i)
            #self.unset_keyword('mw_g', i) # TBD

        self.project.update_thermo_data(self.fluid_species)

    def handle_fluid_species_selection(self):
        row = get_selected_row(self.ui.fluid.tablewidget_fluid_species)
        enabled = (row is not None)
        self.ui.fluid.toolbutton_fluid_species_delete.setEnabled(enabled)
        self.ui.fluid.toolbutton_fluid_species_copy.setEnabled(enabled)

    def fluid_species_add(self):
        sp = self.species_popup
        sp.phases='GL' # ? is this correct
        # how to avoid this if dialog open already?
        self.saved_fluid_species = deepcopy(self.fluid_species) # So we can revert
        sp.cancel.connect(self.fluid_species_revert)
        sp.save.connect(self.fluid_species_save)
        sp.defined_species = self.fluid_species
        sp.update_defined_species()
        sp.setWindowTitle("Fluid Species")
        sp.show()
        sp.raise_()
        sp.activateWindow()

    def fluid_species_delete(self):
        # XXX FIXME this is potentially a big problem since
        # it results in species being renumbered, or a hole in
        # the sequence - either way is trouble.  Have to warn
        # user, if species is referenced elsewhere.
        table = self.ui.fluid.tablewidget_fluid_species
        row = get_selected_row(table)
        if row is None: # No selection
            return
        table.clearSelection()
        key = self.fluid_species.keys()[row]
        del self.fluid_species[key]
        self.update_fluid_species_table()
        # Sigh, we have to update the row in the popup too.
        # Should the popup just be modal, to avoid this?
        sp = self.species_popup
        sp.defined_species = self.fluid_species
        sp.update_defined_species()

    def fluid_species_edit(self):
        table = self.ui.fluid.tablewidget_fluid_species
        row = get_selected_row(table)
        sp = self.species_popup
        self.saved_fluid_species = deepcopy(self.fluid_species) # So we can revert
        sp.cancel.connect(self.fluid_species_revert)
        sp.save.connect(self.fluid_species_save)
        sp.defined_species = self.fluid_species
        sp.update_defined_species()
        if row is None:
            sp.tablewidget_defined_species.clearSelection()
        else:
            sp.tablewidget_defined_species.setCurrentCell(row, 0)
        sp.setWindowTitle("Fluid Species")
        sp.show()
        sp.raise_()
        sp.activateWindow()

    def reset_fluids(self):
        # Set all fluid-related state back to default
        self.saved_fluid_species = None
        self.fluid_species.clear()
        self.init_fluid_default_models()
