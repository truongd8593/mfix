# methods to deal with fluid phase

from __future__ import print_function, absolute_import, unicode_literals, division

from copy import deepcopy
from collections import OrderedDict

import logging
log = logging.getLogger(__name__)

#import Qt
from qtpy import QtWidgets, PYQT5

#local imports
from mfix.constants import *
from mfix.tools.general import set_item_noedit, get_selected_row
from mfix.species_handler import SpeciesHandler

class FluidHandler(SpeciesHandler):
    # Defaults
    def init_fluid_default_models(self):
        self.fluid_density_model = CONSTANT
        self.fluid_viscosity_model = CONSTANT
        self.fluid_mol_weight_model = CONSTANT
        self.fluid_specific_heat_model = CONSTANT
        self.fluid_conductivity_model = AIR
        self.fluid_diffusion_model = AIR


    ## Fluid phase methods
    def enable_fluid_species_eq(self, enabled):
        ui = self.ui.fluid
        for item in (ui.combobox_fluid_diffusion_model,
                     ui.label_fluid_diffusion_model,
                     # more ?
                     ):
            item.setEnabled(enabled)
        # dif_g0 == diffusion coeff model
        items = (ui.lineedit_keyword_dif_g0,
                 ui.label_dif_g0_units)
        for item in items:
            item.setEnabled(enabled and (self.fluid_diffusion_model == CONSTANT))


    def enable_fluid_scalar_eq(self, state):
        ui = self.ui.fluid
        ui.checkbox_enable_scalar_eq.setChecked(state)
        spinbox = ui.spinbox_nscalar_eq
        spinbox.setEnabled(state)
        if state:
            val = spinbox.value()
            self.set_fluid_nscalar_eq(val)
        else:
            # Don't call set_fluid_nscalar_eq(0) b/c that will clobber spinbox
            prev_nscalar = self.fluid_nscalar_eq + self.solids_nscalar_eq
            self.fluid_nscalar_eq = 0
            self.update_scalar_equations(prev_nscalar)


    def set_fluid_nscalar_eq(self, value):
        # This *sums into* nscalar - not a simple keyword
        prev_nscalar = self.fluid_nscalar_eq + self.solids_nscalar_eq
        self.fluid_nscalar_eq = value
        spinbox = self.ui.fluid.spinbox_nscalar_eq
        if value != spinbox.value():
            spinbox.setValue(value)
        self.update_scalar_equations(prev_nscalar)


    def init_fluid_handler(self):
        self.fluid_species = OrderedDict() # keyed by ALIAS

        ui = self.ui.fluid

        ui.lineedit_fluid_phase_name.default_value = self.fluid_phase_name = "Fluid"
        self.saved_fluid_species = None
        self.init_fluid_default_models()
        # Handle a number of cases which are essentially the same
        # see 'set_fluid_mol_weight_model' below to help understand this
        def make_fluid_model_setter(self, name, key):
            def setter(model):
                ui = self.ui.fluid
                setattr(self, name, model) # self.fluid_<name>_model = model
                combobox = getattr(ui, 'combobox_' + name)
                prev_model = combobox.currentIndex()
                if model != prev_model:
                    combobox.setCurrentIndex(model)
                # Make tooltip match setting (for longer names which are truncated)
                combobox.setToolTip(combobox.currentText())

                # Enable lineedit for constant model
                key_g0 = 'c_pg0' if key=='c_p' else key + '_g0'
                key_usr = 'usr_cpg' if key=='c_p' else 'usr_' + key + 'g'
                lineedit = getattr(ui, 'lineedit_keyword_%s' % key_g0)
                label = getattr(ui, 'label_%s_units' % key_g0)

                for item in (lineedit, label):
                    item.setEnabled(model==CONSTANT)

                if model == CONSTANT:
                    value = lineedit.value # Possibly re-enabled gui item
                    if value != '' and self.project.get_value(key_g0) != value:
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

            # Set the combobox default value (?)
            combobox = getattr(ui, 'combobox_'+model_name)
            combobox.default_value = getattr(self, model_name)
            #print(model_name, combobox.default_value)

        combobox = ui.combobox_fluid_mol_weight_model
        combobox.default_value = self.fluid_mol_weight_model

        # more stuff moved from mfix.__init__
        checkbox = ui.checkbox_keyword_species_eq_args_0
        checkbox.clicked.connect(self.enable_fluid_species_eq)

        ui.lineedit_fluid_phase_name.editingFinished.connect(
            self.handle_fluid_phase_name)
        ui.checkbox_enable_scalar_eq.clicked.connect(
            self.enable_fluid_scalar_eq)
        ui.spinbox_nscalar_eq.valueChanged.connect(
            self.set_fluid_nscalar_eq)

        # Fluid phase models
        for name in ('density', 'viscosity', 'specific_heat', 'mol_weight',
                     'conductivity', 'diffusion'):
            model_name = 'fluid_%s_model' % name
            combobox = getattr(ui, 'combobox_%s' % model_name)
            setter = getattr(self,'set_%s' % model_name)
            combobox.currentIndexChanged.connect(setter)

        # Fluid species
        tb = ui.toolbutton_fluid_species_add
        tb.clicked.connect(self.fluid_species_add)
        tb = ui.toolbutton_fluid_species_copy # misnomer
        tb.clicked.connect(self.fluid_species_edit)
        tb.setEnabled(False)
        tb = ui.toolbutton_fluid_species_delete
        tb.setEnabled(False)
        tb.clicked.connect(self.fluid_species_delete)
        tw = ui.tablewidget_fluid_species
        tw.itemSelectionChanged.connect(self.handle_fluid_species_selection)
        self.fixup_fluid_table()


    def fixup_fluid_table(self):
        hv = QtWidgets.QHeaderView
        ui = self.ui.fluid
        table = ui.tablewidget_fluid_species
        if PYQT5:
            resize = table.horizontalHeader().setSectionResizeMode
        else:
            resize = table.horizontalHeader().setResizeMode
        for n in range(table.columnCount()):
            resize(n, hv.ResizeToContents if n>0
                   else hv.Stretch)


    # molecular wt model only has 2 choices, and the key names don't
    # follow the same pattern, so create its setter specially
    def set_fluid_mol_weight_model(self, model):
        ui = self.ui.fluid
        self.fluid_mol_weight_model = model
        combobox = ui.combobox_fluid_mol_weight_model
        # Make tooltip match setting (for longer names which are truncated)
        combobox.setToolTip(combobox.currentText())
        prev_model = combobox.currentIndex()
        if model != prev_model:
            combobox.setCurrentIndex(model)
        # Enable lineedit for constant mol_weight model
        lineedit = ui.lineedit_keyword_mw_avg
        label = ui.label_mw_avg_units
        for item in (lineedit, label):
            item.setEnabled(model==CONSTANT)
        if model == CONSTANT:
            value = lineedit.value # Possibly re-enabled gui item
            if value != '' and self.project.get_value("mw_avg") != value:
                self.set_keyword("mw_avg", value) # Restore keyword value
        else: # Mixture
            # TODO: validate, require mw for all component species
            self.unset_keyword("mw_avg")


    def handle_fluid_phase_name(self): # editingFinished signal does not include value
        value = self.ui.fluid.lineedit_fluid_phase_name.text()
        self.set_fluid_phase_name(value)


    def set_fluid_phase_name(self, value):
        self.fluid_phase_name = value
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
        """Update table in fluid pane.  Also set nmax_g, species_g and species_alias_g keywords,
        which are not tied to a single widget"""
        table = self.ui.fluid.tablewidget_fluid_species
        table.clearContents()
        if self.fluid_species is None:
            self.fixup_fluid_table()
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
            self.unset_keyword('species_g', args=i)
            self.unset_keyword('species_alias_g', args=i)
            #self.unset_keyword('mw_g', i) # TBD

        self.project.update_thermo_data(self.fluid_species)
        self.fixup_fluid_table()


    def handle_fluid_species_selection(self):
        ui = self.ui.fluid
        table = ui.tablewidget_fluid_species
        row = get_selected_row(table)
        enabled = (row is not None)
        ui.toolbutton_fluid_species_delete.setEnabled(enabled)
        ui.toolbutton_fluid_species_copy.setEnabled(enabled)
        if enabled:
            table.doubleClicked.connect(self.fluid_species_edit)
        else:
            try:
                table.doubleClicked.disconnect() #self.fluid_species_edit)
            except:
                pass


    def fluid_species_add(self):
        sp = self.species_popup
        sp.set_phases('GL')
        sp.do_search('') # Init to full db
        # how to avoid this if dialog open already?
        self.saved_fluid_species = deepcopy(self.fluid_species) # So we can revert
        sp.reset_signals()
        sp.cancel.connect(self.fluid_species_revert)
        sp.save.connect(self.fluid_species_save)
        sp.defined_species = self.fluid_species
        sp.extra_aliases = self.fluid_make_extra_aliases()
        sp.update_defined_species()
        sp.setWindowTitle("Fluid Species")
        sp.enable_density(False)
        sp.popup()


    def fluid_make_extra_aliases(self):
        # Construct the 'extra_aliases' set to pass to the species popup
        # Exclude the fluid phase
        aliases = set()
        for ss in self.solids_species.values():
            aliases.update(set(s['alias'] for s in ss.values()))
        return aliases


    def fluid_species_delete(self):
        # XXX FIXME this is potentially a big problem since
        # it results in species being renumbered, or a hole in
        # the sequence - either way is trouble.  Have to warn
        # user, if species is referenced elsewhere.
        ui = self.ui.fluid
        tw = ui.tablewidget_fluid_species
        row = get_selected_row(tw)
        if row is None: # No selection
            return

        table.clearSelection() #?
        alias = tw.item(row,0).text()
        del self.fluid_species[alias]

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
        sp.set_phases('GL')
        self.saved_fluid_species = deepcopy(self.fluid_species) # So we can revert
        sp.reset_signals()
        sp.cancel.connect(self.fluid_species_revert)
        sp.save.connect(self.fluid_species_save)
        sp.defined_species = self.fluid_species
        sp.extra_aliases = self.fluid_make_extra_aliases()
        sp.update_defined_species()
        if row is None:
            sp.tablewidget_defined_species.clearSelection()
        else:
            sp.tablewidget_defined_species.setCurrentCell(row, 0)
        sp.setWindowTitle("Fluid Species")
        sp.enable_density(False)
        sp.popup()


    def setup_fluid(self):
        # Called whenever we switch to fluid tab
        self.P = 0
        ui = self.ui.fluid
        tw = ui.tablewidget_fluid_species
        # Autoselect if unique row
        if get_selected_row(tw) is None and tw.rowCount() == 1:
            tw.setCurrentCell(0,0)

        enabled = (self.fluid_nscalar_eq > 0)
        ui.checkbox_enable_scalar_eq.setChecked(self.fluid_nscalar_eq>0)
        ui.spinbox_nscalar_eq.setEnabled(enabled)


    def reset_fluid(self):
        # Set all fluid-related state back to default
        self.fluid_phase_name = 'Fluid'
        self.saved_fluid_species = None
        self.fluid_species.clear()
        self.init_fluid_default_models()



#Fluid phase Task Pane Window: (unavailable if fluid phase was disable)
#    Option to rename the phase (e.g, air, gas)

#    Option to disable Momentum Equations (enabled by default)
# Sets keyword: MOMENTUM_X/Y/Z_EQ(0)

#    Option to enable Species Equations
# Sets keyword: SPECIES_EQ(0)

#    Option to enable scalar equations
# Define the number of scalar equations
# Value sums into keyword NSCALAR
# Sets keyword PHASE4SCALAR(*)=0 for total listed scalars

#    Select Density Model:
# Selection always available
# Available selections:
#  Constant: [DEFAULT]
#    Selection always available
#    Specify constant gas density, RO_G0
#  Ideal gas law:
#    Selection always available
#    Keyword RO_G0 must be undefined

#    Requires a fluid phase molecular weight
#    Requires temperature field for full domain
#  UDF
#    Selection is always available
#    Sets keyword USR_ROg
#    MFIX runtime check verifies UDF was provided

#Select Viscosity Model:
# Selection always available
# Available selections:
#  Constant: [DEFAULT]
#    Selection always available
#    Specify constant gas viscosity, MU_G0
#  Sutherland's law
#    Selection always available
#    Keyword MU_G0 must be undefined
#    Requires temperature field for full domain
#  UDF
#    Selection always available
#    Sets keyword USR_MUg
#    MFIX runtime check verifies UDF was provided

#Select Molecular Weight Model:
# Selection always available
# Available selections:
#  Constant; [DEFAULT]
#    Specification always available
#    Specify constant molecular weight, MW_AVG
#  Mixture:
#    Selection always available
#    Requires molecular weights for all species components

#Select Specific Heat Model:
# Selection available only when solving thermal energy equations
# Available selections:
#  Constant; [DEFAULT]
#    Selection always available
#    Specify constant fluid phase specific heat, C_PG0
#  Mixture:
#    Selection always available
#    Keyword C_PG0 must be undefined
#    Requires specific heats for all species components
#  UDF
#    Selection always available
#    Sets keyword USR_CPg
#    MFIX runtime check verifies UDF was provided

#Select Thermal Conductivity Model:
# Selection only available when solving thermal energy equations
# Available selections:
#  Constant
#    Selection always available
#    Specify constant thermal conductivity, K_G0
#  Temperature dependent (air); [DEFAULT]
#    Selection always available
#    Keyword K_G0 must be undefined
#  UDF
#    Selection always available
#    Set keyword USR_KG
#    MFIX runtime check verifies UDF was provided

#Select Diffusion Coefficient Model:
# Selection only available when solving species equations
# Available selections:
#  Constant
#    Selection always available
#    Specify a constant diffusion coefficient, DIF_G0
#  Dilute Mixture Approximation (air); [DEFAULT]
#    Selection always available
#    Keyword DIF_G0 must be undefined
#    Requires temperature field for full domain
#  UDF
#    Selection always available
#    Sets keyword USR_DIFG
#    MFIX runtime check verifies UDF was provided

#Fluid phase species selection:
# Species data required under any of the following conditions:
#  Solving species equations
#  Density model is the ideal gas law with mixture molecular weight model
#  Energy equations are solved with mixture specific heat model
# Specification panel operates as a popup window triggered by an Add/Edit button
# Summary window provides a list of the species and an overview of some properties

#Fluid phase Material Database window (popup):
#    Select database (BURCAT); later could link in other databases.
#    Capability to search selected database for chemical name
#    Import from database copies the usable information from the database into a new entry in the 'run database'
#    New creates a new 'blank' species in the 'run database' where the user must supply all the thermochemical data.
#    Delete removes an entry from the 'run database'

#NOTE: The gas phase species molecular weights, MW_G(#) cannot be directly specified. This
#keyword is not needed because users can edit the molecular weight in the material database popup
#window.
