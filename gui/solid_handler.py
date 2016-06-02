# methods to deal with solids, split off from gui.py

from __future__ import print_function, absolute_import, unicode_literals, division

from collections import OrderedDict

#import Qt
from qtpy import QtCore, QtWidgets, PYQT5

#local imports
from constants import *
from tools.general import set_item_noedit, get_selected_row, widget_iter

class SolidHandler(object):

    def init_solid_handler(self):
        pass

    def setup_combobox_solids_model(self):
        """solids model combobox is tied to solver setting"""
        solver = self.project.solver
        if solver == SINGLE:
            # Note, if Single-Phase solver is enabled, this pane is disabled
            return

        cb = self.ui.solids.combobox_solids_model
        model = cb.model()
        #          TFM,  DEM,  PIC
        enabled = [False, False, False]
        enabled[0] = (solver==TFM or solver==HYBRID)
        enabled[1] = (solver==DEM or solver==HYBRID)
        enabled[2] = (solver==PIC)
        for (i, e) in enumerate(enabled):
            item = model.item(i, 0)
            flags = item.flags()
            if e:
                flags |= QtCore.Qt.ItemIsEnabled
            else:
                flags &= ~QtCore.Qt.ItemIsEnabled
            item.setFlags(flags)
        # Set for current phase
        if self.solids_current_phase is not None:
            model = self.project.get_value("solids_model", args = self.solids_current_phase)
            i = 0 if model=='TFM' else 1 if model=='DEM' else 2 if model=='PIC' else None
            if i is None:
                i = 0 if solver in (TFM, HYBRID) else 1 if solver == DEM else 2

            if not enabled[i]:
                # Current selection no longer valid, so pick first valid choice
                # Don't leave a non-enabled item selected!
                i = enabled.index(True)
            cb.setCurrentIndex(i)
        else: # Set based on overall solver
            i = 0 if solver in (TFM, HYBRID) else 1 if solver == DEM else 2
            cb.setCurrentIndex(i)

    def handle_combobox_solids_model(self, index):
        if self.solids_current_phase is None:
            return # shouldn't get here

        phase = self.solids_current_phase
        name, data = self.solids.items()[phase-1] # FIXME, use SpeciesCollection not OrderedDict here
        model = ('TFM', 'DEM', 'PIC')[index]
        self.update_keyword('solids_model', model, args=self.solids_current_phase)
        data['model'] = model
        self.update_solids_table()
        self.update_solids_detail_pane()

    def make_solids_name(self, n):
        while True:
            name = 'Solid %d' % n
            if name not in self.solids:
                break
            n += 1
        return name

    def solids_add(self):
        tw = self.ui.solids.tablewidget_solids
        nrows = tw.rowCount()
        n = nrows + 1
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
        self.update_solids_table()
        tw.setCurrentCell(nrows, 0) # Select new item

    def handle_solids_table_selection(self):
        s = self.ui.solids
        tw = s.tablewidget_solids
        row = get_selected_row(tw)
        enabled = (row is not None)
        s.toolbutton_solids_delete.setEnabled(enabled)
        s.toolbutton_solids_copy.setEnabled(enabled)
        name = None if row is None else tw.item(row,0).text()
        self.solids_current_phase = (row+1) if row is not None else None
        self.update_solids_detail_pane()

    def update_solids_detail_pane(self):
        """update the solids detail pane for currently selected solids phase.
        if no solid phase # selected, pane is cleared and disabled"""
        s = self.ui.solids
        sa = s.scrollarea_solids_detail
        phase = self.solids_current_phase
        if phase is None: #name is None or phase is None: # current solid phase name.
            # Disable all inputs
            sa.setEnabled(False)
            for item in widget_iter(sa):
                if isinstance(item, QtWidgets.QCheckBox):
                    item.setChecked(False)
                # Clear out all values?
            return

        name = list(self.solids.keys())[phase-1] # ugh
        solid = self.solids[name]
        model = solid['model']

        # Enable the input areas, initialize to values for current solid
        sa.setEnabled(True)
        s.lineedit_solids_phase_name.setText(name)
        self.setup_combobox_solids_model()

        # Inialize all the line edit widgets
        def as_str(x):
            return '' if x is None else str(x)
        s.lineedit_keyword_d_p0_args_S.setText(as_str(solid['diameter']))
        s.lineedit_keyword_ro_s0_args_S.setText(as_str(solid['density']))
        # And the checkboxes
        for key in ('momentum_x_eq', 'momentum_y_eq', 'momentum_z_eq'):
            cb = getattr(s, 'checkbox_keyword_%s_args_S'%key)
            if model == 'TFM': # momentum eq only avail w/ TFM solid model
                val = self.project.get_value(key, default=True, args=phase)
                cb.setEnabled(True)
                cb.setChecked(val)
            else:
                cb.setEnabled(False)
                cb.setChecked(False)
        key = 'species_eq'
        val = self.project.get_value(key, default=False, args=phase)
        cb = getattr(s, 'checkbox_keyword_%s_args_S'%key)
        cb.setChecked(val)

        nscalar = self.project.get_value('nscalar', 0)
        nscalar_phase = sum(1 for i in range(1, nscalar+1)
                            if self.project.get_value('phase4scalar', args=i) == phase)
        saved_nscalar_eq = solid.get('saved_nscalar_eq', 0)
        s.spinbox_nscalar_eq.setValue(saved_nscalar_eq)
        enabled = solid.get('enable_scalar_eq', False) # (nscalar_phase > 0)
        s.checkbox_enable_scalar_eq.setChecked(enabled)
        self.enable_solid_scalar_eq(enabled)


    def update_solids_table(self):
        hv = QtWidgets.QHeaderView
        table = self.ui.solids.tablewidget_solids
        if PYQT5:
            resize = table.horizontalHeader().setSectionResizeMode
        else:
            resize = table.horizontalHeader().setResizeMode
        for n in range(4):
            resize(n, hv.ResizeToContents if n>0
                   else hv.Stretch)

        if self.solids is None:
            table.clearContents()
            return
        nrows = len(self.solids)
        table.setRowCount(nrows)

        # helper fn
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item

        #Update the internal table from keywords
        # TODO: remove this, use SolidsCollection
        for (i, (k,v)) in enumerate(self.solids.items(), 1):
            for (myname, realname) in (('model', 'solids_model'),
                                       ('diameter', 'd_p0'),
                                       ('density', 'ro_s0')):
                self.solids[k][myname] = self.project.get_value(realname, args=i)

        for (row,(k,v)) in enumerate(self.solids.items()):
            table.setItem(row, 0, make_item(k))
            for (col, key) in enumerate(('model', 'diameter', 'density'), 1):
                table.setItem(row, col, make_item(v[key]))

        if nrows == 1: # If there's only 1 let's select it for the user's convenience
            table.setCurrentCell(0,0)
    def handle_solids_phase_name(self):
        new_name = self.ui.solids.lineedit_solids_phase_name.text()
        phase = self.solids_current_phase
        if phase is None:
            return
        old_name = self.solids.keys()[phase-1] # Ugh
        if new_name in self.solids: # Reject the input
            self.ui.solids.lineedit_solids_phase_name.setText(old_name)
            return

        # Rewriting dict to change key while preserving order
        d = OrderedDict() # Ugh.  Use SpeciesCollection!
        for (k,v) in self.solids.iteritems():
            if k==old_name:
                k = new_name
            d[k] = v
        self.solids = d
        self.update_solids_table()
        self.project.mfix_gui_comments['solids_phase_name(%s)'%phase] = new_name
        self.set_unsaved_flag()

    def solids_delete(self):
        tw = self.ui.solids.tablewidget_solids
        row = get_selected_row(tw)
        if row is None: # No selection
            return
        name = tw.item(row, 0).text()
        del self.solids[name]
        tw.removeRow(row)
        tw.clearSelection()
        self.set_unsaved_flag()

    def enable_solid_scalar_eq(self, state):
        spinbox = self.ui.solids.spinbox_nscalar_eq
        spinbox.setEnabled(state)
        phase = self.solids_current_phase
        if phase is None:
            return

        name = list(self.solids.keys())[phase-1] # ugh
        solid = self.solids[name]

        current_state = solid.get('enable_scalar_eq')

        if state == current_state: # avoid clobbering saved values
            return

        solid['enable_scalar_eq'] = state
        if state:
            value = solid.get('saved_nscalar_eq', 0)
            self.set_solid_nscalar_eq(value)
        else:
            # Don't call set_solid_nscalar_eq(0) b/c that will clobber spinbox
            prev_nscalar = self.fluid_nscalar_eq + self.solid_nscalar_eq
            solid['saved_nscalar_eq'] = solid.get('nscalar_eq', 0)
            solid['nscalar_eq'] = 0
            self.solid_nscalar_eq = sum(s.get('nscalar_eq', 0) for s in self.solids.values())
            self.update_scalar_equations(prev_nscalar)

    def set_solid_nscalar_eq(self, value):
        # This *sums into* nscalar - not a simple keyword
        phase = self.solids_current_phase
        if phase is None:
            return
        name = list(self.solids.keys())[phase-1] # ugh
        solid = self.solids[name]

        nscalar = self.project.get_value('nscalar', 0)
        prev_nscalar = self.fluid_nscalar_eq + self.solid_nscalar_eq

        solid['nscalar_eq'] = value
        solid['saved_nscalar_eq'] = value
        # would it be better to just use nscalar - fluid_nscalar_eq?
        self.solid_nscalar_eq = sum(s.get('nscalar_eq', 0) for s in self.solids.values())

        spinbox = self.ui.solids.spinbox_nscalar_eq
        if value != spinbox.value():
            spinbox.setValue(value)
            return
        self.update_scalar_equations(prev_nscalar)


    def solids_change_tab(self, tabnum, btn):
        """ switch solids stacked widget based on selected """
        self.animate_stacked_widget(
            self.ui.solids.stackedwidget_solids,
            self.ui.solids.stackedwidget_solids.currentIndex(),
            tabnum,
            direction='horizontal',
            line=self.ui.solids.line_solids,
            to_btn=btn,
            btn_layout=self.ui.solids.gridlayout_solid_tab_btns)

    def reset_solids(self):
        # Set all solid-related state back to default
        self.solids_current_phase = None
        self.solids.clear()
        self.update_solids_table()
        self.update_solids_detail_pane()
