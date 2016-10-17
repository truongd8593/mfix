# methods to deal with solids, split off from mfixgui.py

from __future__ import print_function, absolute_import, unicode_literals, division

from copy import deepcopy
from collections import OrderedDict

import logging
log = logging.getLogger(__name__)

#import Qt
from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtCore import Qt
QT_MAX = 16777215
UserRole = QtCore.Qt.UserRole

# ./model/param_mod.f:67:      INTEGER, PARAMETER :: DIM_M = 10 # max # of solids phases
DIM_M = 10

#local imports
from mfixgui.constants import *
from mfixgui.tools.general import (set_item_noedit, get_selected_row,
                           widget_iter,
                           format_key_with_args,
                           get_combobox_item, set_item_enabled,
                           drop_row_column_triangular, append_row_column_triangular )


from mfixgui.tools import keyword_args

from mfixgui.widgets.base import LineEdit

from mfixgui.solids_tfm import SolidsTFM
from mfixgui.solids_dem import SolidsDEM
from mfixgui.solids_pic import SolidsPIC

from mfixgui.species_handler import SpeciesHandler

class SolidsHandler(SolidsTFM, SolidsDEM, SolidsPIC, SpeciesHandler):

    def init_solids_default_models(self):
        self.solids_density_model = CONSTANT
        self.solids_viscosity_model = OTHER
        self.solids_mol_weight_model = MIXTURE
        self.solids_specific_heat_model = CONSTANT
        self.solids_conductivity_model = OTHER


    def init_solids_handler(self):
        self.solids = OrderedDict() # keyed by name of solids phase
        self.solids_current_phase = None
        self.solids_species = {} #dict of OrderedDict, keyed by phase num
        #  The value dict  is keyed by species, but probably should
        #  be keyed by alias FIXME.

        self.solids_current_tab = 0 # Materials

        ui = self.ui.solids
        tb = ui.toolbutton_solids_add
        tb.clicked.connect(self.solids_add)
        tb = ui.toolbutton_solids_delete
        tb.clicked.connect(self.solids_delete)
        tb.setEnabled(False)
        tw_solids, tw_solids_species = ui.tablewidget_solids, ui.tablewidget_solids_species
        # Force summary table to update on kw updates
        class TableWidgetProxy:
            def objectName(self):
                return "proxy"
            def updateValue(*args):
                self.update_solids_table()
        # TODO the solids_table should just be directly editable
        self.project.register_widget(TableWidgetProxy(),
                             ['solids_model', 'd_p0', 'ro_s0'], args='*')
        tw_solids.itemSelectionChanged.connect(self.handle_solids_table_selection)
        cb = ui.combobox_solids_model
        cb.activated.connect(self.handle_combobox_solids_model)
        cb.setToolTip(cb.currentText())
        ui.lineedit_solids_phase_name.editingFinished.connect(self.handle_solids_phase_name)
        ui.checkbox_enable_scalar_eq.stateChanged.connect(self.enable_solids_scalar_eq)
        ui.spinbox_nscalar_eq.valueChanged.connect(self.set_solids_nscalar_eq)

        self.saved_solids_species = None
        self.init_solids_default_models()
        ui.checkbox_keyword_species_eq_args_P.clicked.connect(
            self.handle_solids_species_eq)

        # Handle a number of cases which are essentially the same
        # avoid repetition in set_solids_*_model methods
        def make_solids_model_setter(self, name, key):
            def setter(model=None): # If not specified, determine correct model from keywords
                ui = self.ui.solids
                phase = self.solids_current_phase
                key_s0 = 'c_ps0' if key=='c_p' else key + '_s0'
                key_usr = 'usr_cps' if key=='c_p' else 'usr_' + key + 's'
                lineedit = getattr(ui, 'lineedit_keyword_%s_args_P' % key_s0)
                units_label = getattr(ui, 'label_%s_units' % key_s0)

                if model is None and phase is not None:
                    val_s0 = self.project.get_value(key_s0, args=[phase])
                    val_usr = self.project.get_value(key_usr, args=[phase])
                    model = (CONSTANT if val_s0 is not None
                             else UDF if val_usr is not None
                             else VARIABLE)
                    if val_s0 is not None:
                        lineedit.setText(str(val_s0))

                setattr(self, name, model) # self.solids_<name>_model = model
                combobox = getattr(ui, 'combobox_' + name)
                if model != combobox.currentIndex():
                    combobox.setCurrentIndex(model)
                # Make tooltip match setting (for longer names which are truncated)
                combobox.setToolTip(combobox.currentText())

                # Enable lineedit for constant model
                enabled = (model==CONSTANT)
                for item in (lineedit, units_label):
                    item.setEnabled(enabled)

                if phase is None:
                    return

                if model == CONSTANT:
                    # FIXME: we could have a saved value from a different phase
                    #  (Clear out lineedits on phase change?)
                    value = lineedit.value # Possibly re-enabled gui item
                    if value=='None':
                        value = ''
                    if value != '' and self.project.get_value(key_s0, args=phase) != value:
                        self.update_keyword(key_s0, value, args=[phase]) # Restore keyword value
                    if value == '':
                        val_s0 = self.project.get_value(key_s0, args=[phase])
                        if val_s0 == 'None': # Should not happen, but filter out if found
                            self.unset_keyword(key_s0, args=[phase])
                            val_s0 = None
                        if val_s0 is not None:
                            lineedit.setText(str(val_s0))
                elif model == UDF:
                    self.unset_keyword(key_s0, args=phase)
                    self.update_keyword(key_usr, True, args=phase)
                else: # Continuum, mixture, etc
                    self.unset_keyword(key_s0, args=phase)
                    self.unset_keyword(key_usr, args=phase)

                if name == 'solids_density_model': # Extra handling: Density column changes
                    self.update_solids_baseline_groupbox(model) # availability
                    self.update_solids_table() # 'density' column changes

            return setter


        for (name, key) in (
                ('density', 'ro'),
                ('viscosity', 'mu'),
                ('specific_heat', 'c_p'),
                ('conductivity', 'k')):

            model_name = 'solids_%s_model' % name
            setter = make_solids_model_setter(self, model_name, key)
            setattr(self, 'set_%s' % model_name, setter)
            combobox = getattr(ui, 'combobox_'+model_name)
            combobox.activated.connect(setter)
            combobox.default_value = getattr(self, model_name)
            #print(model_name, combobox.default_value)


        # Solids species
        tb = ui.toolbutton_solids_species_add
        tb.clicked.connect(self.solids_species_add)
        tb = ui.toolbutton_solids_species_copy # misnomer - edit
        tb.clicked.connect(self.solids_species_edit)
        tb.setEnabled(False)
        tb = ui.toolbutton_solids_species_delete
        tb.setEnabled(False)
        tb.clicked.connect(self.solids_species_delete)
        tw_solids_species.itemSelectionChanged.connect(self.handle_solids_species_selection)

        # Advanced
        cb = ui.checkbox_disable_close_pack
        cb.clicked.connect(self.disable_close_pack)
        self.add_tooltip(cb, 'close_packed')

        cb = ui.checkbox_enable_added_mass_force
        cb.clicked.connect(self.enable_added_mass_force)
        self.add_tooltip(cb, 'added_mass')

        # connect solid tab buttons
        for i, btn in enumerate((ui.pushbutton_solids_materials,
                                 ui.pushbutton_solids_tfm,
                                 ui.pushbutton_solids_dem,
                                 ui.pushbutton_solids_pic)):
            btn.pressed.connect(lambda i=i, btn=btn: self.solids_change_tab(i, btn))

        for tw in (ui.tablewidget_solids, ui.tablewidget_solids_species,
                   ui.tablewidget_solids_baseline):
            self.fixup_solids_table(tw)

        self.init_solids_tfm()
        self.init_solids_dem()
        self.init_solids_pic()

    def solids_update_tabs(self):
        ui = self.ui.solids
        solver = self.project.solver
        if solver == SINGLE: # We shouldn't be here
            return

        items = (ui.pushbutton_solids_tfm,
                 ui.pushbutton_solids_dem,
                 ui.pushbutton_solids_pic)

        enable_dict = {TFM: (True, False, False),
                       DEM: (False, True, False),
                       PIC: (False, False, True),
                       HYBRID: (True, True, False)}

        enabled = enable_dict[solver] if self.solids else (False, False, False)

        for (item, e) in zip(items, enabled):
            item.setDisabled(not e)

        # Don't stay on a disabled tab!
        # (Do we ever disable "Materials"? - don't think so)
        active_tab = ui.stackedwidget_solids.currentIndex()
        if active_tab > 0 and not enabled[active_tab-1]:
            if solver==SINGLE:
                i, p = 0, ui.pushbutton_solids_materials
            elif solver in (TFM, HYBRID):
                i, p = 1, ui.pushbutton_solids_tfm
            elif solver==DEM:
                i, p = 2, ui.pushbutton_solids_dem
            elif solver==PIC:
                i, p = 3, ui.pushbutton_solids_pic
            self.solids_change_tab(i, p)

    # Solids sub-pane navigation
    def solids_change_tab(self, tabnum, to_btn):
        ui = self.ui.solids
        self.solids_current_tab = tabnum
        self.animate_stacked_widget(
            ui.stackedwidget_solids,
            ui.stackedwidget_solids.currentIndex(),
            tabnum,
            direction='horizontal',
            line=ui.line_solids,
            to_btn=to_btn,
            btn_layout=ui.gridlayout_solid_tab_btns)
        self.setup_solids_tab(tabnum)
        for btn in (ui.pushbutton_solids_pic, ui.pushbutton_solids_dem, ui.pushbutton_solids_tfm,
                    ui.pushbutton_solids_materials):
            btn.setChecked(btn==to_btn)
            font = btn.font()
            font.setBold(btn==to_btn)
            btn.setFont(font)


    def setup_solids(self):
        self.setup_solids_tab(self.solids_current_tab)


    def setup_solids_tab(self, tabnum):
        if tabnum == 1:
            self.setup_tfm_tab()
        elif tabnum == 2:
            self.setup_dem_tab()
        elif tabnum == 3:
            self.setup_pic_tab()
        else: # Materials
            self.update_solids_table()
            self.update_solids_detail_pane()

    # Advanced
    def disable_close_pack(self, val):
        ui = self.ui.solids
        cb = ui.checkbox_disable_close_pack
        if val != cb.isChecked(): # not from a click action
            cb.setChecked(val)
            return
        phase = self.solids_current_phase
        if phase is None:
            return
        close_packed = self.project.get_value('close_packed', default=None, args=phase)
        if (close_packed is not False) and val: # Disabling - popup as per SRS p15
            resp=self.message(text="Disabling close-packing for %s\nAre you sure?" % self.solids_current_phase_name,
                              buttons=['yes','no'],
                              default = 'no')
            if resp != 'yes':
                cb.setChecked(False)
                return

        if close_packed is None and not val:
            return # We will leave the keyword at its default unset value

        self.update_keyword('close_packed', not val, args=phase)

    def enable_added_mass_force(self, val):
        ui = self.ui.solids
        cb = ui.checkbox_enable_added_mass_force
        if val != cb.isChecked(): # not from a click action
            cb.setChecked(val)
            return
        phase = self.solids_current_phase
        if phase is None:
            return
        # Warn when unsetting added mass for other phases
        prev_phase = self.project.get_value('m_am')
        if prev_phase is not None and prev_phase != phase and val:
            prev_phase_name = list(self.solids.keys())[prev_phase-1]
            resp=self.message(text="Disabling added mass force for %s\nAre you sure?" % prev_phase_name,
                              buttons=['yes','no'],
                              default = 'no')
            if resp == 'no':
                cb.setChecked(False)
                return

        if val:
            self.update_keyword('m_am', phase)
            self.update_keyword('added_mass', True)

        else:
            self.unset_keyword('m_am')
            self.unset_keyword('added_mass')


    def fixup_solids_table(self, tw, stretch_column=0):
        # Should we just hide the entire table (including header) if no rows?
        ui = self.ui.solids
        hv = QtWidgets.QHeaderView
        if PYQT5:
            resize = tw.horizontalHeader().setSectionResizeMode
        else:
            resize = tw.horizontalHeader().setResizeMode
        ncols = tw.columnCount()
        for n in range(0, ncols):
            resize(n, hv.Stretch if n==stretch_column else hv.ResizeToContents)

        # trim excess vertical space - can't figure out how to do this in designer
        header_height = tw.horizontalHeader().height()

        # TODO FIXME scrollbar handling is not right - scrollbar status can change
        # outside of this function.  We need to call this everytime window geometry changes
        scrollbar_height = tw.horizontalScrollBar().isVisible() * (4+tw.horizontalScrollBar().height())
        nrows = tw.rowCount()

        if nrows==0:
            height = header_height+scrollbar_height
        else:
            height =  (header_height+scrollbar_height
                       + nrows*tw.rowHeight(0) + 4) # extra to avoid unneeded scrollbar
        if tw == ui.tablewidget_solids: # In a splitter
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(header_height)
            ui.top_frame.setMaximumHeight(height+40)
            ui.top_frame.setMinimumHeight(header_height+40)
            ui.top_frame.updateGeometry()
        else:
            tw.setMaximumHeight(height) # Works for tablewidget inside groupbox
            tw.setMinimumHeight(height) #? needed for tablewidget_des_en_input. should we allow scrollbar?
        tw.updateGeometry() #? needed?

    def handle_solids_species_eq(self, enabled):
        ui = self.ui.solids
        if not enabled:
            self.set_solids_density_model(CONSTANT)
            #ui.combobox_solids_density_model.setCurrentIndex(CONSTANT)
        set_item_enabled(get_combobox_item(ui.combobox_solids_density_model,
                                                VARIABLE), enabled)
        self.update_solids_species_groupbox() # availability

        # TODO integrate with BCs
        #When solving solids species equations:
        #    Set keyword BC_HW_X_S(#,#,#) to 0.0
        #    Set keyword BC_C_X_S(#,#,#) to 0.0
        #    Set keyword BC_XW_S(#,#,#) to UNDEFINED
        # getting this right is tricky - how about when we add phases & species?

    def setup_combobox_solids_model(self):
        """solids model combobox is tied to solver setting"""
        solver = self.project.solver
        if solver == SINGLE:
            # Note, if Single-Phase solver is enabled, this pane is disabled
            return
        ui = self.ui.solids
        cb = ui.combobox_solids_model
        model = cb.model()
        #          TFM,  DEM,  PIC
        enabled = [False, False, False]
        enabled[0] = (solver==TFM or solver==HYBRID)
        enabled[1] = (solver==DEM or solver==HYBRID)
        enabled[2] = (solver==PIC)
        for (i, e) in enumerate(enabled):
            set_item_enabled(get_combobox_item(cb,i), e)

        # Set for current phase
        if self.solids_current_phase is not None:
            mod = self.project.get_value("solids_model", args=self.solids_current_phase)
            i = 0 if mod=='TFM' else 1 if mod=='DEM' else 2 if mod=='PIC' else None
            if i is None:
                i = 0 if solver in (TFM, HYBRID) else 1 if solver == DEM else 2

            if not enabled[i]:
                # Current selection no longer valid, so pick first valid choice
                # Don't leave a non-enabled item selected!
                i = enabled.index(True)
            cb.setCurrentIndex(i)
            cb.setToolTip(cb.currentText())
        else: # Set based on overall solver
            i = 0 if solver in (TFM, HYBRID) else 1 if solver == DEM else 2
            cb.setCurrentIndex(i)
            cb.setToolTip(cb.currentText())


    def handle_combobox_solids_model(self, index):
        ## NB:  Solids model is not the same as solver!
        # Solver values are enumerated in constants.py.  Solids models are strings, 'TFM', 'DEM', 'PIC'

        if self.solids_current_phase is None:
            return # shouldn't get here

        phase = self.solids_current_phase
        name, data = list(self.solids.items())[phase-1]
        model = ('TFM', 'DEM', 'PIC')[index]
        self.update_keyword('solids_model', model, args=self.solids_current_phase)
        data['model'] = model
        self.update_solids_table()
        self.update_solids_detail_pane()
        self.update_nav_tree() # PSs depends on solids_model


    def make_solids_name(self, n):
        while True:
            name = 'Solid %d' % n
            if name not in self.solids:
                break
            n += 1
        return name

    def solids_add(self):
        ui = self.ui.solids
        tw = ui.tablewidget_solids
        nrows = tw.rowCount()
        n = nrows + 1 # index of new solid
        tw.setRowCount(n)
        name = self.make_solids_name(n)
        if self.project.solver == SINGLE: # Should not get here! this pane is disabled.
            return
        else:
            model = [None, 'TFM', 'DEM', 'PIC', 'TEM'][self.project.solver]
        diameter = 0.0
        density = 0.0
        self.update_keyword('solids_model', model, args=n)
        self.solids[name] = {'model': model,
                             'diameter': diameter, # TODO: diameter is REQUIRED
                             'density': density} # more?
        self.solids_species[n] = OrderedDict()
        self.update_keyword('mmax', len(self.solids))
        #self.update_keyword('nmax_s', 0, args=[n])
        self.update_keyword('species_eq', False, args=[n]) # Disable species eq by default, per SRS
        self.update_solids_table()
        tw.setCurrentCell(nrows, 0) # Select new item

        # ICs enabled/disabled depends on number of solids
        self.update_nav_tree()

        # Tabs enable/disable depending on number of solids
        self.solids_update_tabs()
        if len(self.solids) >= DIM_M:
            ui.toolbutton_solids_add.setEnabled(False)


    def handle_solids_table_selection(self):
        ui = self.ui.solids
        self.P = self.solids_current_phase # should already be set
        tw = ui.tablewidget_solids
        row = get_selected_row(tw)
        enabled = (row is not None)
        ui.toolbutton_solids_delete.setEnabled(enabled)
        #ui.toolbutton_solids_copy.setEnabled(enabled) - removed from ui
        self.solids_current_phase_name = None if row is None else tw.item(row,0).text()
        self.solids_current_phase = (row+1) if row is not None else None
        self.update_solids_detail_pane()


    def update_solids_detail_pane(self):
        """update the solids detail pane for currently selected solids phase.
        if no solid phase # selected, pane is cleared and disabled"""
        ui = self.ui.solids
        phase = self.P = self.solids_current_phase
        if phase is None: #name is None or phase is None: # current solid phase name.
            # Disable all inputs
            self.update_solids_species_table()
            self.fixup_solids_table(ui.tablewidget_solids)
            ui.detail_pane.setEnabled(False)
            for item in widget_iter(ui.detail_pane):
                if isinstance(item, QtWidgets.QCheckBox):
                    item.setChecked(False)
                if isinstance(item, QtWidgets.QLineEdit):
                    # Surprise, a spinbox includes a lineedit!
                    if 'spinbox' in item.objectName(): # 'qt_spinbox_lineedit':
                        item.setText('1')
                    else:
                        item.setText('')
            return
        name = list(self.solids.keys())[phase-1]
        solid = self.solids[name]
        model = solid['model']

        # Enable the input areas, initialize to values for current solid
        ui.detail_pane.setEnabled(True)
        ui.lineedit_solids_phase_name.setText(self.solids_current_phase_name)
        self.setup_combobox_solids_model()

        # Inialize all the line edit widgets
        def as_str(x):
            return '' if x is None else str(x)
        # Why not get these from keywords?
        def get_widget(key):
            return getattr(ui, 'lineedit_keyword_%s_args_P' % key)
        for key in ('d_p0', 'ro_s0', 'des_em'): # use a widget iterator?
            get_widget(key).setText(as_str(self.project.get_value(key, args=phase)))

        # And the checkboxes
        for key in ('momentum_x_eq', 'momentum_y_eq', 'momentum_z_eq'):
            cb = getattr(ui, 'checkbox_keyword_%s_args_P'%key)
            if model == 'TFM': # momentum eq only avail w/ TFM solid model
                val = self.project.get_value(key, default=True, args=phase)
                cb.setEnabled(True)
                cb.setChecked(val)
            else:
                cb.setEnabled(False)
                cb.setChecked(False)

        # Set species eq checkbox to correct value
        key = 'species_eq'
        species_eq = self.project.get_value(key, default=True, args=[phase])
        cb = getattr(ui, 'checkbox_keyword_%s_args_P'%key)
        cb.setChecked(species_eq)

        # Deal with scalar eq
        sb = ui.spinbox_nscalar_eq
        nscalar = self.project.get_value('nscalar', default=0)
        nscalar_phase = sum(1 for i in range(1, nscalar+1)
                            if self.project.get_value('phase4scalar', args=i) == phase)
        saved_nscalar_eq = solid.get('saved_nscalar_eq', 1)

        if sb.value() != saved_nscalar_eq:
            # set value without triggering callback
            # Maybe this is too much trouble to go through,
            # just to save/restore per-phase nscalar
            sb.valueChanged.disconnect()
            sb.setValue(saved_nscalar_eq)
            sb.valueChanged.connect(self.set_solids_nscalar_eq)
        enabled = (nscalar_phase > 0)
        ui.checkbox_enable_scalar_eq.setChecked(enabled)
        sb.setEnabled(enabled)
        if enabled:
            sb.setValue(nscalar_phase)
        #self.enable_solids_scalar_eq(enabled)

        ### Restrictions (see SRS p13)

        # Variable density model requires species equations
        cb_density =  ui.combobox_solids_density_model
        enabled = bool(species_eq)
        set_item_enabled(get_combobox_item(cb_density, VARIABLE), enabled)
        ro_s0 = self.project.get_value('ro_s0', args=[phase], default=None)
        self.set_solids_density_model(VARIABLE if (enabled and ro_s0 is None) else CONSTANT)

        # Viscosity only available for TFM solids
        self.set_solids_viscosity_model()
        enabled = (model=='TFM')
        for item in (ui.label_solids_viscosity_model,
                     ui.combobox_solids_viscosity_model):
            item.setEnabled(enabled)
        if not enabled:
            for item in (ui.lineedit_keyword_mu_s0_args_P,
                         ui.label_mu_s0_units):
                item.setEnabled(False)

        # Mol. wt is locked to MIXTURE

        # Specific heat only available when solving energy eq
        self.set_solids_specific_heat_model()
        energy_eq = self.project.get_value('energy_eq', default=True)
        enabled = (energy_eq==True)
        for item in (ui.combobox_solids_specific_heat_model,
                     ui.label_solids_specific_heat_model):
            item.setEnabled(enabled)
        if not enabled:
            for item in (ui.lineedit_keyword_c_ps0_args_P,
                         ui.label_c_ps0_units):
                item.setEnabled(False)

        # Thermal Conductivity Model:
        # Selection only available for MFIX-TFM solids model
        # Selection only available when solving thermal energy equations
        cb = ui.combobox_solids_conductivity_model
        if model != 'TFM':
            self.set_solids_conductivity_model(CONSTANT)
            for (i,e) in enumerate((True,False,False)):
                set_item_enabled(get_combobox_item(cb,i), e)
        else:
            self.set_solids_conductivity_model()
            for (i,e) in enumerate((True,True,True)):
                set_item_enabled(get_combobox_item(cb,i), e)
        enabled = bool(energy_eq)
        for item in (ui.combobox_solids_conductivity_model,
                     ui.label_solids_conductivity_model):
            item.setEnabled(enabled)
        if not enabled:
            for item in (ui.lineedit_keyword_k_s0_args_P,
                         ui.label_k_s0_units):
                item.setEnabled(False)

        # Specify solids phase emissivity
        # Selection only available for MFIX-DEM solids model
        enabled = (model=='DEM')
        for item in (ui.label_des_em,
                     ui.lineedit_keyword_des_em_args_P):
            item.setEnabled(enabled)

        # Species input is in its own function
        self.update_solids_species_groupbox()
        # as is the baseline table
        self.update_solids_baseline_groupbox(self.solids_density_model)

        self.update_solids_species_table()
        self.fixup_solids_table(ui.tablewidget_solids)

        # Advanced
        enabled = (model=='TFM')
        ui.groupbox_advanced.setEnabled(enabled)
        if enabled:
            close_packed = self.project.get_value('close_packed', default=True, args=phase)
            self.disable_close_pack(not(close_packed))
        else:
            ui.checkbox_disable_close_pack.setChecked(False)
            ui.checkbox_enable_added_mass_force.setChecked(False)

        added_mass = self.project.get_value('added_mass', default=False)
        m_am = self.project.get_value('m_am', default=None)
        self.enable_added_mass_force(added_mass and (m_am == phase))


    def update_solids_table(self):
        ui = self.ui.solids
        table = ui.tablewidget_solids
        if self.solids is None:
            table.clearContents()
            self.unset_keyword('mmax')
            return

        nrows = len(self.solids)
        table.setRowCount(nrows)

        # helper fn
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item

        #Update the internal table from keywords
        for (i, (k,v)) in enumerate(self.solids.items(), 1):
            for (myname, realname) in (('model', 'solids_model'),
                                       ('diameter', 'd_p0'),
                                       ('density', 'ro_s0')):
                self.solids[k][myname] = self.project.get_value(realname, args=i)
                #if self.solids[k]['density'] is None and self.solids_density_model==VARIABLE:
                #    self.solids[k]['density'] = "Variable"


        for (row,(k,v)) in enumerate(self.solids.items()):
            table.setItem(row, 0, make_item(k))
            for (col, key) in enumerate(('model', 'diameter', 'density'), 1):
                table.setItem(row, col, make_item(v[key]))

        if nrows == 1: # If there's only 1 let's auto-select it for the user's convenience
            table.setCurrentCell(0,0)

        # trim excess horizontal space - can't figure out how to do this in designer
        header_height = table.horizontalHeader().height()
        if nrows==0:
            # 34 px is empirical, fixme, should calc. row height
            table.setMaximumHeight(header_height + 34)
        else:
            table.setMaximumHeight(header_height+nrows*table.rowHeight(0) + 4)
            # a little extra to avoid horiz scrollbar when not needed


    def handle_solids_phase_name(self):
        ui = self.ui.solids
        new_name = ui.lineedit_solids_phase_name.text()
        phase = self.solids_current_phase
        if phase is None:
            return
        old_name = list(self.solids.keys())[phase-1]
        if new_name in self.solids: # Reject the input
            self.ui.solids.lineedit_solids_phase_name.setText(old_name)
            return
        self.solids_current_phase_name = new_name
        # rewriting dict to change key while preserving order - hack
        d = OrderedDict()
        for (k,v) in self.solids.items():
            if k==old_name:
                k = new_name
            d[k] = v
        self.solids = d
        self.update_solids_table()
        self.project.mfix_gui_comments['solids_phase_name(%s)'%phase] = new_name
        self.set_unsaved_flag()
        return True

    def solids_delete(self):
        ui = self.ui.solids
        tw = ui.tablewidget_solids
        row = get_selected_row(tw)

        if row is None: # No selection
            return

        phase = row+1
        phase_name = tw.item(row,0).text()

        key = self.solids_phase_in_use(phase)
        if key:
            self.message(text="%s is referenced by %s" %
                         (phase_name, key))
            return

        self.solids_current_phase = self.P = None
        self.solids_current_phase_name = None

        tw.itemSelectionChanged.disconnect() # Avoid selection callbacks during delete.  Reenabled below
        try:
            name = tw.item(row, 0).text()
            tw.removeRow(row)
            del self.solids[name] # This is an ordered dict, keyed by name - no 'hole' problem

            if len(self.solids) < DIM_M:
                ui.toolbutton_solids_add.setEnabled(True)

            self.update_keyword('mmax', len(self.solids))

            # Clear out all keywords related to deleted phase
            self.solids_delete_phase_keys(phase)

            # Fix holes
            # Fixup phase names in mfix_gui_comments
            # TODO make some util functions to manage mfix_gui_comments in a less ad-hoc manner
            for (k,v) in list(self.project.mfix_gui_comments.items()):
                if k.startswith('solids_phase_name('):
                    del self.project.mfix_gui_comments[k]
            for (i,name) in enumerate(self.solids.keys(), 1):
                self.project.mfix_gui_comments['solids_phase_name(%s)'%i] = name

            # species dictionary (also need to remap keys)
            for n in range(phase, len(self.solids)+1):
                self.solids_species[n] = self.solids_species[n+1]
            del self.solids_species[len(self.solids_species)]

            # fix nscalar
            nscalar = self.project.get_value('nscalar', default=0)
            key = 'phase4scalar'
            vals = [self.project.get_value(key, default=0, args=i) for i in range(1, nscalar+1)]
            new_vals = [v if v<phase else v-1 for v in vals if v != phase]
            for (i, val) in enumerate(new_vals, 1):
                self.update_keyword(key, val, args=i)
            for i in range(len(new_vals)+1, nscalar+1):
                self.unset_keyword(key, args=i)
            self.update_keyword('nscalar', nscalar)
            #  TODO fix initial conditions for scalars

            # Fix hole in restitution coeffs
            n = len(self.solids)
            for key in  ('des_et_input', 'des_en_input'):
                prev_size = ((n+1)*(n+2))//2 # Size before row deleted
                vals = [self.project.get_value(key, args=i)
                        for i in range(1, 1+prev_size)]
                if any(v is not None for v in vals):
                    new_vals = drop_row_column_triangular(vals, n+1, phase)
                    for (i, val) in enumerate(new_vals, 1):
                        self.update_keyword(key, val, args=i)
                    for i in range(len(new_vals)+1,  len(vals)+1):
                        self.unset_keyword(key, args=i)
            for key in ('des_et_wall_input', 'des_en_wall_input'):
                prev_size = n+1
                vals = [self.project.get_value(key, args=i)
                        for i in range(1, 1+prev_size)]
                if any(v is not None for v in vals):
                    new_vals = vals[:]
                    del new_vals[phase-1] # 1-based
                    for (i, val) in enumerate(new_vals, 1):
                        self.update_keyword(key, val, args=i)
                    for i in range(len(new_vals)+1,  len(vals)+1):
                        self.unset_keyword(key, args=i)

            self.update_solids_table()

            # ICs enabled/disabled depends on nscalar & number of solids
            self.update_nav_tree()

            # Tabs enable/disable depending on number of solids
            self.solids_update_tabs()
        finally:
            # selection callbacks were disabled
            tw.itemSelectionChanged.connect(self.handle_solids_table_selection)
            self.handle_solids_table_selection()
            self.set_unsaved_flag()


    def enable_solids_scalar_eq(self, state):
        ui = self.ui.solids
        spinbox = ui.spinbox_nscalar_eq
        spinbox.setEnabled(state)
        phase = self.solids_current_phase
        if phase is None:
            return

        name = list(self.solids.keys())[phase-1]
        solid = self.solids[name]
        if state:
            value = spinbox.value()
            self.set_solids_nscalar_eq(value)
        else:
            # Don't call set_solids_nscalar_eq(0) b/c that will clobber spinbox
            prev_nscalar = self.fluid_nscalar_eq + self.solids_nscalar_eq
            solid['saved_nscalar_eq'] = solid.get('nscalar_eq', 0)
            solid['nscalar_eq'] = 0
            self.solids_nscalar_eq = sum(s.get('nscalar_eq', 0) for s in self.solids.values())
            self.update_scalar_equations(prev_nscalar)


    def update_solids_species_groupbox(self):
        """enable/disable species tables based on state"""
        # Species data required under any of the following conditions:
        #  Solving species equations
        #  Energy equations are solved with mixture specific heat model
        phase = self.solids_current_phase
        ui = self.ui.solids

        if phase is None:
            enabled = False
        else:
            species_eq = self.project.get_value('species_eq', args=phase, default=True)
            energy_eq = self.project.get_value('energy_eq', default=True)
            enabled = (species_eq) or (energy_eq and self.solids_specific_heat_model == MIXTURE)

        # Is it a good idea to have hidden items? No.
        #ui.groupbox_species.setVisible(enabled)
        ui.groupbox_species.setEnabled(enabled)
        # Buttons seem to take up a lot of space when table is shrunk
        ui.frame_add_delete_copy_species.setVisible(enabled)

        #tw = ui.tablewidget_solids_species
        #if not enabled: # Hide species?  shrink input area?
        #    tw.clearContents()
        #    tw.setMaximumHeight(32) #?
        #else:
        #    tw.setMaximumHeight(QT_MAX)
        #    self.update_solids_species_table()

    def update_solids_baseline_groupbox(self, density_model):
        #Baseline (unreacted) composition selection:
        # Available only for variable solids density model
        ui = self.ui.solids
        enabled = density_model==VARIABLE
        ui.groupbox_baseline.setEnabled(enabled)
        if not enabled:
            ui.tablewidget_solids_baseline.clearContents()
            ui.tablewidget_solids_baseline.setRowCount(0)
        else:
            self.update_solids_baseline_table()

    def handle_solids_mass_fraction(self, widget, value_dict, args):
        phase = self.solids_current_phase
        if phase is None:
            return
        ui = self.ui.solids
        key = 'x_s0'
        val = value_dict[key]

        table = ui.tablewidget_solids_baseline
        widget.updateValue(key, val)
        if val == '':
            self.unset_keyword(key, args=[phase]+args)
        else:
            self.update_keyword(key, val, args=[phase]+args)
        self.update_solids_mass_fraction_total()

    def update_solids_mass_fraction_total(self):
        ui = self.ui.solids
        key = 'x_s0'
        phase = self.solids_current_phase
        if phase is None:
            return
        table = ui.tablewidget_solids_baseline
        if table.rowCount() == 0:
            return
        total = sum(float(self.project.get_value(key, default=0.0, args=[phase,i]))
                    for i in range(1,len(self.solids_species[phase])+1))
        item = table.item(table.rowCount()-1, 1)
        font = item.font()
        font.setBold(True)
        item.setFont(font)
        item.setText(str(total))

    def handle_solids_inert_species(self, species_index, val):
        phase = self.solids_current_phase
        if phase is None:
            return
        if val:
            self.update_keyword('inert_species', species_index, args=phase)
        else:
            self.unset_keyword('inert_species', args=phase)

        for (i, cb) in enumerate(self.solids_inert_species_checkboxes, 1):
            if i != species_index:
                cb.setChecked(False)


    def update_solids_baseline_table(self):
        ui = self.ui.solids
        table = ui.tablewidget_solids_baseline
        phase = self.solids_current_phase
        if phase is None:
            return
        table.clearContents()
        if self.solids_species[phase]:
            nrows = len(self.solids_species[phase])+1
        else:
            nrows = 0
        table.setRowCount(nrows) # "Total" row at end
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item

        inert_species = self.project.get_value('inert_species', default=None, args=phase)
        self.solids_inert_species_checkboxes = []
        for (row, (species,data)) in enumerate(self.solids_species[phase].items()):
            alias = data.get('alias', species) # default to species if no alias
            table.setItem(row, 0, make_item(alias))

            # mass fraction
            le = LineEdit()
            le.setdtype('dp')
            le.setValInfo(min=0.0, max=1.0)
            key = 'x_s0'
            le.key = key
            le.args = [row+1]
            val = self.project.get_value(key, args=[phase, row+1], default=None)
            if val is not None:
                le.updateValue(key, val)

            le.value_updated.connect(self.handle_solids_mass_fraction)
            table.setCellWidget(row, 1, le)

            # "Inert" checkbox
            cb = QtWidgets.QCheckBox()
            self.solids_inert_species_checkboxes.append(cb)
            cb.setChecked(row+1 == inert_species)
            cb.clicked.connect(lambda val, species_index=row+1:
                               self.handle_solids_inert_species(species_index, val))
            # center checkbox in cell
            layout =QtWidgets.QHBoxLayout()
            layout.addWidget(cb)
            layout.setAlignment(Qt.AlignCenter);
            layout.setContentsMargins(0,0,0,0);
            widget = QtWidgets.QWidget()
            widget.setLayout(layout);
            table.setCellWidget(row, 2, widget)

        if nrows > 0:
            table.setItem(nrows-1, 0, make_item("Total"))
            table.setItem(nrows-1, 1, make_item(''))
            item = table.item(nrows-1, 0)
            font = item.font()
            font.setBold(True)
            item.setFont(font)
            self.update_solids_mass_fraction_total()

        self.fixup_solids_table(table)

    def set_solids_nscalar_eq(self, value):
        # This *sums into* nscalar - not a simple keyword
        phase = self.solids_current_phase
        if phase is None:
            return
        ui = self.ui.solids
        name = list(self.solids.keys())[phase-1]
        solid = self.solids[name]

        nscalar = self.project.get_value('nscalar', default=0)
        prev_nscalar = self.fluid_nscalar_eq + self.solids_nscalar_eq

        solid['nscalar_eq'] = value
        solid['saved_nscalar_eq'] = value
        self.solids_nscalar_eq = sum(s.get('nscalar_eq', 0) for s in self.solids.values())

        spinbox = ui.spinbox_nscalar_eq
        if value != spinbox.value():
            spinbox.valueChanged.disconnect()
            spinbox.setValue(value)
            spinbox.valueChanged.connect(self.set_solids_nscalar_eq)
            return
        self.update_scalar_equations(prev_nscalar)



    # --- solids species methods ---
    def solids_species_revert(self):
        if self.saved_solids_species is None:
            return
        phase = self.solids_current_phase
        if phase is None:
            return
        self.solids_species[phase] = self.saved_solids_species
        self.species_popup.defined_species = deepcopy(self.solids_species[phase])
        self.update_solids_species_table()
        self.update_solids_baseline_groupbox(self.solids_density_model)

    def solids_species_save(self):
        phase = self.solids_current_phase
        if phase is None:
            return
        ui = self.ui.solids
        self.solids_species[phase] = deepcopy(self.species_popup.defined_species)
        self.update_solids_species_table()
        self.update_solids_baseline_groupbox(self.solids_density_model)
        self.fixup_solids_table(ui.tablewidget_solids_species)
        self.fixup_solids_table(ui.tablewidget_solids_baseline)

    def update_solids_species_table(self):
        """Update table in solids pane.  Also sets nmax_s, species_s and species_alias_s keywords,
        which are not tied to a single widget"""
        ui = self.ui.solids
        table = ui.tablewidget_solids_species
        table.clearContents()
        phase = self.solids_current_phase
        if phase is None:
            return
        if self.solids_species.get(phase) is None:
            return

        nrows = len(self.solids_species[phase])
        table.setRowCount(nrows)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        old_nmax_s = self.project.get_value('nmax_s', default=0, args=phase)
        nmax_s = len(self.solids_species[phase])
        self.update_keyword('nmax_s', nmax_s, args=phase)

        for (row, (species,data)) in enumerate(self.solids_species[phase].items()):
            # order must match table headers
            args = [phase,row+1]
            for (col, key) in enumerate(('alias', 'phase', 'density', 'mol_weight',
                                        'h_f', 'source')):
                alias = data.get('alias', species) # default to species if no alias
                data['alias'] = alias # for make_item
                table.setItem(row, col, make_item(data.get(key)))
                self.update_keyword('species_s', species, args=args)
                self.update_keyword('species_alias_s', alias, args=args)
                density = data.get('density')
                if density is not None:
                    self.update_keyword('ro_xs0', density, args=args)
                else:
                    self.unset_keyword('ro_xs0', args=args)

                # We're avoiding mw_s in favor of the settings in THERMO DATA
                #self.update_keyword('mw_s', data['mol_weight'], args=row+1)#

        # Clear any keywords with indices above nmax_s NEEDED?
        # TODO use keyword_args here
        for i in range(nmax_s+1, old_nmax_s+1):
            for kw in ('species_s', 'species_alias_s', 'x_s0', 'ro_xs0'):
                self.unset_keyword(kw, args=[phase,i])
            #self.unset_keyword('mw_s', args=[phase,i]) # TBD FIXME

        # FIXME, what's the right place for this?
        #self.project.update_thermo_data(self.solids_species[phase])
        self.fixup_solids_table(ui.tablewidget_solids_species)
        self.fixup_solids_table(ui.tablewidget_solids_baseline)

        # Autoselect unique row
        for tw in (ui.tablewidget_solids_species, ui.tablewidget_solids_baseline):
            if tw.rowCount()==1 and get_selected_row(tw) is None:
                tw.setCurrentCell(0, 0)

    def handle_solids_species_selection(self):
        ui = self.ui.solids
        table = self.ui.solids.tablewidget_solids_species
        row = get_selected_row(table)
        enabled = (row is not None)
        ui.toolbutton_solids_species_delete.setEnabled(enabled)
        ui.toolbutton_solids_species_copy.setEnabled(enabled) # edit
        if enabled:
            table.doubleClicked.connect(self.solids_species_edit)
        else:
            try:
                table.doubleClicked.disconnect() #self.solids_species_edit)
            except:
                pass

    def solids_species_add(self):
        phase = self.solids_current_phase
        if phase is None:
            return
        # Workaround for deleted phase
        if phase not in self.solids_species:
            self.solids_species[phase] = OrderedDict()
        sp = self.species_popup
        sp.set_phases('SC')
        sp.do_search('')
        # how to avoid this if dialog open already?
        self.saved_solids_species = deepcopy(self.solids_species[phase]) # So we can revert
        sp.reset_signals()
        sp.cancel.connect(self.solids_species_revert)
        sp.save.connect(self.solids_species_save)
        sp.defined_species = self.solids_species[phase]
        sp.extra_aliases = self.solids_make_extra_aliases(phase)
        sp.update_defined_species()
        sp.setWindowTitle("Species for %s" %self.solids_current_phase_name)
        sp.enable_density(True)
        sp.popup()

    def solids_make_extra_aliases(self, phase):
        # Construct the 'extra_aliases' set to pass to the species popup
        # Exclude the specified phase
        aliases = set(f['alias'] for f in self.fluid_species.values())
        for (p, ss) in self.solids_species.items():
            if p == phase:
               continue
            aliases.update(s['alias'] for s in ss.values())
        return aliases

    def solids_species_delete(self):
        phase = self.solids_current_phase
        if phase is None:
            return
        ui = self.ui.solids

        tw = ui.tablewidget_solids_species
        row = get_selected_row(tw)
        if row is None: # No selection
            return
        tw.clearSelection() #?
        alias = tw.item(row,0).text()
        if alias in self.solids_species[phase]:
            del self.solids_species[phase][alias] # XXXX Fixme

        self.update_solids_species_table()
        self.update_solids_baseline_groupbox(self.solids_density_model)
        self.fixup_solids_table(ui.tablewidget_solids_species)
        self.fixup_solids_table(ui.tablewidget_solids_baseline)

        # Sigh, we have to update the row in the popup too.
        # Should the popup just be modal, to avoid this?
        sp = self.species_popup
        sp.defined_species = self.solids_species[phase]
        sp.update_defined_species()

    def solids_species_edit(self):
        phase = self.solids_current_phase
        if phase is None:
            return
        ui = self.ui.solids
        table = ui.tablewidget_solids_species
        row = get_selected_row(table)
        sp = self.species_popup
        sp.set_phases('SC')

        self.saved_solids_species = deepcopy(self.solids_species[phase]) # So we can revert
        sp.reset_signals()
        sp.cancel.connect(self.solids_species_revert)
        sp.save.connect(self.solids_species_save)
        sp.defined_species = self.solids_species[phase]
        sp.extra_aliases = self.solids_make_extra_aliases(phase)
        sp.update_defined_species()
        if row is None:
            sp.do_search('')
            sp.tablewidget_defined_species.clearSelection()
        else:
            # Initialize search box to current species (?)
            sp.do_search(list(self.solids_species[phase].keys())[row])
            sp.tablewidget_defined_species.setCurrentCell(row, 0)
        sp.setWindowTitle("Species for %s" % self.solids_current_phase_name)
        sp.enable_density(True)
        sp.popup()

    def solids_phase_in_use(self, phase):
        for k in keyword_args.keys_by_type['phase']:
            indices = self.project.get_key_indices(k)
            if not indices:
                continue
            arg_types = keyword_args.keyword_args[k]

            if arg_types == ['phase']:
                # Single reference is OK.  We will
                # delete these keys along with the phase in solids_delete_phase_keys
                continue
            if arg_types == ['ic', 'phase']:
                for args in indices:
                    if args[1] != phase:
                        continue
                    # TODO check if value is default, if so allow delete
                    if self.project.get_value(k, args=args) is not None:
                        return format_key_with_args(k,args)

        return None # Ok to delete phase, no refs

    def solids_delete_phase_keys(self, phase):
        prev_size = len(self.solids) + 1 # Size before row deleted
        for key in keyword_args.keys_by_type['phase']:
            arg_types = keyword_args.keyword_args[key]
            indices = self.project.get_key_indices(key)
            if not indices:
                continue

            if arg_types == ['phase']:
                # Copy/slide/trim
                vals = [self.project.get_value(key, args=i)
                        for i in range(1, 1+prev_size)]
                del vals[phase-1] # Slide (1-based)
                for (i, val) in enumerate(vals, 1):
                    self.update_keyword(key, val, args=i)
                self.unset_keyword(key, args=prev_size) #Trim off the end

            # the only phase-phase keyword is R_P,  this code doesn't handle that yet
            elif arg_types == ['phase', 'phase']: # FIXME
                raise TypeError(key, arg_types)

            else:
                # Multidimensional copy-and-slide, using dict instead of list
                phase_pos = arg_types.index('phase')
                new_vals = {}
                for args in indices:
                    args_phase = args[phase_pos]
                    if args_phase == phase:
                        continue # skip the phase we're deleting
                    new_args = list(args)
                    if args_phase > phase:
                        new_args[phase_pos] -= 1 #Slide along 'phase_pos' axis
                    new_vals[tuple(new_args)] = self.project.get_value(key, args=args)
                for (args, val) in new_vals.items():
                    self.update_keyword(key, val, args=args)
                for args in indices: # Trim
                    key_phase = args[phase_pos]
                    if key_phase == prev_size:
                        self.unset_keyword(key, args)


    def reset_solids(self):
        ui = self.ui.solids
        # Set all solid-related state back to default
        self.solids_current_phase = self.P = None
        self.solids_current_phase_name = None
        self.saved_solids_species = None
        self.solids.clear()
        self.solids_species.clear()
        self.update_solids_table()
        self.update_solids_detail_pane()
        self.solids_current_tab = 0
        self.solids_change_tab(0, ui.pushbutton_solids_materials)
        # TODO (?)  reset THERMO_DATA ?