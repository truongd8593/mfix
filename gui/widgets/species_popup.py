#!/usr/bin/env python
from __future__ import print_function, absolute_import, unicode_literals, division

"""Species selector dialog for MFIX GUI, includes stand-alone test"""

import os
import sys
import signal
from collections import OrderedDict
import pickle

from qtpy import QtCore, QtWidgets, PYQT5, uic
from qtpy.QtWidgets import QTableWidgetItem, QTableWidget, QLineEdit
from qtpy.QtGui import QValidator, QDoubleValidator
from qtpy.QtCore import QObject, QEvent
UserRole = QtCore.Qt.UserRole

from tools.general import (set_item_noedit, set_item_enabled,
                           get_selected_row, get_selected_rows,
                           widget_iter)

if PYQT5:
    def resize_column(table, col, flags):
        table.horizontalHeader().setSectionResizeMode(col, flags)
else:
    def resize_column(table, col, flags):
        table.horizontalHeader().setResizeMode(col, flags)


class SpeciesPopup(QtWidgets.QDialog):

    save = QtCore.Signal()
    cancel = QtCore.Signal()

    def load_burcat(self, path):
        if not os.path.exists(path):
            print("%s not found, create it by running read_burcat.py" % path)
            sys.exit(-1)
        with open(path, 'rb') as f:
            db = pickle.load(f)
        by_phase = {}

        for k,v in db.items():
            phase = k[1]
            if phase not in by_phase:
                by_phase[phase] = {}
            key = k[0], k[2], k[3]
            by_phase[phase][key] = v
        self.db = by_phase

        # build search list, lowercased
        self.haystack = []
        self.comments = {}
        for phase in 'GLSC':
            htmp = [((k[0].lower(), v[2].lower()), k, phase) for (k,v) in self.db[phase].items()]
            htmp.sort()
            self.haystack.extend(htmp)
            # comment fields
            self.comments[phase] = dict((k, v[2])
                                        for (k,v) in self.db[phase].items())

    def do_search(self, string):
        lineedit = self.ui.lineedit_search
        if string != lineedit.text():
            lineedit.setText(string)
            return

        results = {}
        self.ui.tablewidget_search.clearContents()
        results = []
        match_empty = True
        if match_empty or string:
            needle = string.lower()
            for ((k, key, phase)) in self.haystack:
                if (phase in self.phases and
                    (needle in k[0] or
                     (self.include_comments and needle in k[1]))):
                    results.append((key, phase))
        table = self.ui.tablewidget_search
        nrows = len(results)
        self.ui.tablewidget_search.clearContents()
        self.ui.tablewidget_search.setRowCount(nrows)
        self.search_results = [None]*nrows

        # http://stackoverflow.com/questions/10192579/
        table.model().blockSignals(True);
        for (i,r) in enumerate(results):
            key, phase = r
            comment = self.comments[phase][key]
            item = QTableWidgetItem(key[0])
            item.setToolTip(comment)
            set_item_noedit(item)
            table.setItem(i, 0, item)
            item = QTableWidgetItem(phase)
            set_item_noedit(item)
            table.setItem(i, 1, item)
            self.search_results[i] = (key, phase)
        table.model().blockSignals(False);

    def get_species_data(self, key, phase):
        """exposes species database to external clients"""
        db = self.db.get(phase)
        if not db:
            return None
        # FIXME, this is inefficient.  remove tmin/tmax from key tuple.
        #  also, it's possible that there are multiple definitions for the
        #  same species, with different temp. ranges.  This just returns
        #  the first one
        for (keytuple, data) in db.items():
            (species, tmin, tmax) = keytuple
            if species == key:
                (coeffs, mol_weight, comment) = data
                a_low = coeffs[:7]
                a_high = coeffs[7:14]
                h_f = coeffs[14]
                return {'source': 'BURCAT',
                        'phase': phase,
                        'mol_weight': mol_weight,
                        'h_f': h_f,
                        'tmin':  tmin,
                        'tmax': tmax,
                        'a_low': a_low,
                        'a_high': a_high}

    def handle_search_selection(self):
        row = get_selected_row(self.tablewidget_search)
        self.ui.pushbutton_import.setEnabled(row is not None)

    def handle_include_comments(self, val):
        self.include_comments = val
        self.do_search(self.ui.lineedit_search.text())

    def clear_species_panel(self):
        for item in self.species_panel_items:
            item.setEnabled(False)
            if hasattr(item, 'setText'):
                item.setText('')
        table = self.ui.tablewidget_params
        for row in range(8):
            for col in range(2):
                w = table.cellWidget(row, col)
                if w:
                    w.setText('')
                    #table.cellWidget(i,j).setText('')

    def enable_species_panel(self):
        for item in self.species_panel_items:
            item.setEnabled(True)
        species = self.current_species
        data = self.defined_species.get(species)
        ui = self.ui
        def make_item(val, key=None):
            item = QLineEdit(table)
            item.setText(str(val))
             #item = QTableWidgetItem(str(val))
            item.setValidator(QDoubleValidator(item))
            item.setFrame(False)
            if key:
                item.editingFinished.connect(make_handler(item=item, key=key))
            return item
        self.ui.combobox_phase.setEnabled(False)
        def make_handler(item, key):
            def handler(item=item,key=key):
                if not self.current_species:
                    print("Error, no current species")
                val = item.text()
                try:
                    data = self.defined_species[self.current_species]
                    val = float(val)
                    if type(key) is tuple:
                        data[key[0]][key[1]] = val
                    else:
                        data[key] = val
                    if key != 'density':
                        data['source'] = 'User Defined' # Data has been modified
                        # Density is not in BURCAT so we don't consider this a mod (?)
                    ui.label_species_source.setText(data['source'])
                except ValueError:
                    # reset field to prev. value
                    pass

            return handler
        ui.label_species_source.setText(data['source'])
        ui.label_species.setText(species)
        ui.combobox_phase.setCurrentIndex('GLSC'.index(data['phase']))
        ui.lineedit_alias.setText(data['alias'])
        ui.lineedit_mol_weight.setText(str(data['mol_weight']))
        ui.lineedit_mol_weight.editingFinished.connect(make_handler(ui.lineedit_mol_weight,'mol_weight'))
        ui.lineedit_h_f.setText(str(data['h_f']))
        ui.lineedit_h_f.editingFinished.connect(make_handler(ui.lineedit_h_f,'h_f'))
        if self.density_enabled:
            density = data.get('density')
            ui.lineedit_density.setText('' if density is None else str(density))
            ui.lineedit_density.editingFinished.connect(make_handler(ui.lineedit_density,'density'))

        table = ui.tablewidget_params
        table.setCellWidget(0, 0, make_item(data['tmin'], key='tmin'))
        table.setCellWidget(0, 1, make_item(data['tmax'], key='tmax'))
        for (i,x) in enumerate(data['a_low']):
            table.setCellWidget(i, 0, make_item(x, key=('a_low', i)))
        for (i,x) in enumerate(data['a_high']):
            table.setCellWidget(i, 1, make_item(x, key=('a_high', i)))


    def enable_density(self, enabled):
        self.density_enabled = enabled
        ui = self.ui
        for w in (ui.label_density, ui.label_density_units):
            w.setEnabled(enabled)
        if not enabled:
            if ui.lineedit_density in self.species_panel_items:
                self.species_panel_items.remove(ui.lineedit_density)
            ui.lineedit_density.clear()
        else:
            self.species_panel_items.append(ui.lineedit_density)

    def handle_defined_species_selection(self):
        self.ui.tablewidget_search.clearSelection() # is this right?
        table = self.tablewidget_defined_species
        row = get_selected_row(table)

        if row is None:
            self.current_species = None
            self.clear_species_panel()
            self.ui.pushbutton_delete.setEnabled(False)
            self.ui.combobox_phase.setEnabled(False)
        else:
            self.ui.pushbutton_delete.setEnabled(True)
            self.current_species = table.item(row, 0).data(UserRole)
            self.enable_species_panel()

    def make_alias(self, species):
        return species.split(' ', 1)[0]

    def make_user_species_name(self):
        n=1
        while ("Species %d" % n) in self.user_species_names:
            n += 1
        name = "Species %d" % n
        self.user_species_names.add(name)
        return name

    def do_import(self):
        rows = get_selected_rows(self.tablewidget_search)
        for row in rows:
            self.do_import_row(row)

    def do_import_row(self, row):
        self.ui.combobox_phase.setEnabled(False)
        rowdata = self.search_results[row]
        key, phase = rowdata
        data = self.db[phase][key]
        (species, tmin, tmax) = key
        (coeffs, mol_weight, comment) = data
        # Do we want to allow same species in different phases, etc?
        if species in self.defined_species:
            return # Don't allow duplicates

        alias = self.make_alias(species)
        a_low = coeffs[:7]
        a_high = coeffs[7:14]
        h_f = coeffs[14]

        species_data = {'source': 'BURCAT',
                        'phase': phase,
                        'alias': alias,
                        'mol_weight': mol_weight,
                        'h_f': h_f,
                        'tmin':  tmin,
                        'tmax': tmax,
                        'a_low': a_low,
                        'a_high': a_high}
        if self.density_enabled:
            species_data['density'] = None # ? where do we get this?
        self.defined_species[species] = species_data
        self.add_defined_species_row(species, select=True)
        self.set_save_button(True)

    def update_defined_species(self):
        self.tablewidget_defined_species.clearSelection()
        self.tablewidget_defined_species.setRowCount(0)
        for s in self.defined_species.keys():
            self.add_defined_species_row(s, select=False)

    def add_defined_species_row(self, species, select=False):
        species_data = self.defined_species[species]
        ui = self.ui
        table = ui.tablewidget_defined_species
        nrows = table.rowCount()
        table.setRowCount(nrows+1)
        alias = species_data['alias']
        phase = species_data['phase']
        item = QTableWidgetItem(alias)
        set_item_noedit(item)
        item.setData(UserRole, species)
        table.setItem(nrows, 0, item)
        item = QTableWidgetItem(phase)
        set_item_noedit(item)
        table.setItem(nrows, 1, item)

        if select:
            table.setCurrentCell(nrows, 0) # Cause the new row to be selected

    def handle_delete(self):
        table = self.ui.tablewidget_defined_species
        row = get_selected_row(table)
        if row is None: # No selection
            return
        current_species = self.current_species # will be reset when selection cleared
        table.removeRow(row)
        table.clearSelection()
        self.ui.tablewidget_search.clearSelection()

        if current_species:
            del self.defined_species[current_species]
        if current_species in self.user_species_names:
            self.user_species_names.remove(current_species)
        self.current_species = None
        self.clear_species_panel()
        self.set_save_button(True)
        self.ui.combobox_phase.setEnabled(False)

    def handle_new(self):
        phase = self.default_phase
        species = self.make_user_species_name()
        alias = species.replace(' ', '') #?
        mol_weight = 0
        density = None
        h_f = 0
        tmin = 200.0 # ?
        tmax = 600.0 # ?
        a_low = [0.0]*7
        a_high = [0.0]*7

        species_data = {'source': 'User Defined',
                        'phase': phase,
                        'alias': alias,
                        'mol_weight': mol_weight,
                        'density': density,
                        'h_f': h_f,
                        'tmin':  tmin,
                        'tmax': tmax,
                        'a_low': a_low,
                        'a_high': a_high}

        self.defined_species[species] = species_data
        self.current_species = species
        self.enable_species_panel()
        self.add_defined_species_row(species, select=True)
        lineedit = self.ui.lineedit_alias
        lineedit.selectAll()
        lineedit.setFocus()
        self.ui.combobox_phase.setEnabled(True)

    def handle_alias(self):
        val = self.ui.lineedit_alias.text()
        table = self.ui.tablewidget_defined_species
        row = get_selected_row(table)
        if row is None: # No selection
            return
        #NB making a new item here, instead of changing item inplace
        item = QTableWidgetItem(val)
        item.setData(UserRole, self.current_species)
        set_item_noedit(item)
        table.setItem(row, 0, item)
        self.defined_species[self.current_species]['alias'] = val

    def set_save_button(self, state):
        self.ui.pushbutton_save.setEnabled(state)

    def handle_combobox_phase(self, index):
        phase = 'GLSC'[index]
        if not self.current_species:
            return
        species =  self.defined_species[self.current_species]
        species['phase'] = phase

    def reset_signals(self):
        # todo:  fix this so it's not the caller's responsibility
        #  (make a util function that calls this & pops up dialog)
        for sig in (self.cancel, self.save):
            try:
                sig.disconnect()
            except:
                pass # isSignalConnected only exists in qt5.


    def handle_phase(self):
        phases = ''
        for phase in 'GLSC':
            button = getattr(self.ui, 'pushbutton_%s' % phase)
            if button.isChecked():
                phases += phase
        if phases == self.phases:
            return
        self.phases = phases
        self.default_phase = phases[0] if phases else ''
        self.do_search(self.ui.lineedit_search.text())

    def __init__(self, app, parent=None, phases='GLCS'):
        super(SpeciesPopup, self).__init__(parent)
        self.app = app
        self.phases = phases
        self.include_comments = False
        self.default_phase = phases[0] if phases else ''
        self.density_enabled = True
        thisdir = os.path.abspath(os.path.dirname(__file__))
        datadir = thisdir
        self.load_burcat(os.path.join(datadir, 'burcat.pickle'))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        ui = self.ui = uic.loadUi(os.path.join(uidir, 'species_popup.ui'), self)

        # key=species, val=data tuple.  can add phase to key if needed
        self.defined_species = OrderedDict()
        self.extra_aliases = set() # To support enforcing
        # uniqueness of aliases across all phases

        self.search_results = []
        self.user_species_names = set()

        # Set up UI
        ui.lineedit_search.textChanged.connect(self.do_search)
        ui.pushbutton_import.clicked.connect(self.do_import)
        ui.pushbutton_import.setEnabled(False)
        ui.tablewidget_search.itemSelectionChanged.connect(
            self.handle_search_selection)
        ui.tablewidget_defined_species.itemSelectionChanged.connect(
            self.handle_defined_species_selection)

        ui.pushbutton_new.clicked.connect(self.handle_new)
        ui.pushbutton_delete.clicked.connect(self.handle_delete)
        ui.checkbox_include_comments.clicked.connect(self.handle_include_comments)

        for phase in 'GLSC':
            button = getattr(self.ui, 'pushbutton_%s' % phase)
            button.clicked.connect(self.handle_phase)

        ui.combobox_phase.currentIndexChanged.connect(self.handle_combobox_phase)
        #http://stackoverflow.com/questions/15845487/how-do-i-prevent-the-enter-key-from-closing-my-qdialog-qt-4-8-1
        # Do not use buttonbox.  https://mfix.netl.doe.gov/gitlab/develop/mfix/issues/101
        buttons = (ui.pushbutton_save, ui.pushbutton_cancel)
        buttons[0].clicked.connect(lambda: (self.save.emit(), self.close()))
        buttons[1].clicked.connect(lambda: (self.cancel.emit(), self.close()))

        class AliasValidator(QValidator):

            # Make sure aliases are unique
            def __init__(self, parent=None):
                super(AliasValidator, self).__init__()
                self.parent = parent

            def validate(self, text, pos):
                if text=="":
                    self.parent.set_save_button(False)
                    return (QValidator.Intermediate, text, pos)
                if " " in text: # ?is space allowed?
                    self.parent.set_save_button(False)
                    # More visual indication of invalid alias?
                    return (QValidator.Invalid, text, pos)
                current_species = self.parent.current_species
                current_alias = self.parent.defined_species[current_species]['alias']
                for (k,v) in self.parent.defined_species.items():
                    if v['alias'] == current_alias: # Skip selected item
                        continue
                    if v['alias'] == text:
                        self.parent.set_save_button(False)
                        # More visual indication of invalid alias?
                        return (QValidator.Intermediate, text, pos)
                if text in self.parent.extra_aliases:
                    self.parent.set_save_button(False)
                    # More visual indication of invalid alias?
                    return (QValidator.Intermediate, text, pos)
                self.parent.set_save_button(True)
                return (QValidator.Acceptable, text, pos)

        lineedit = ui.lineedit_alias
        lineedit.setValidator(AliasValidator(parent=self))
        lineedit.editingFinished.connect(self.handle_alias)

        for l in (ui.lineedit_mol_weight,
                  ui.lineedit_h_f,
                  ui.lineedit_density):
            l.setValidator(QDoubleValidator())

        self.species_panel_items=[
            ui.label_species_source,
            ui.label_species,
            ui.lineedit_alias,
            ui.lineedit_mol_weight,
            ui.lineedit_h_f,
            ui.combobox_specific_heat_model,
            ui.lineedit_density,
            ui.tablewidget_params]

        hv = QtWidgets.QHeaderView
        for tw in (self.tablewidget_search, self.tablewidget_defined_species):
            resize_column(tw, 0, hv.Stretch)
            resize_column(tw, 1, hv.ResizeToContents)
        tw = self.tablewidget_params
        for i in (0, 1):
            resize_column(tw, i, hv.Stretch)

        self.set_save_button(False) # nothing to Save
        self.clear_species_panel()

    def set_phases(self, phases):
        if phases == self.phases:
            return
        self.phases = phases
        for phase in 'GLSC':
            button = getattr(self.ui, 'pushbutton_%s' % phase)
            button.setChecked(phase in phases)
        self.default_phase = phases[0] if phases else ''
        self.do_search(self.ui.lineedit_search.text())


    def popup(self):
        self.show()
        self.raise_()
        self.activateWindow()


if __name__ == '__main__':
    args = sys.argv
    qapp = QtWidgets.QApplication(args)
    dialog = QtWidgets.QDialog()
    species_popup = SpeciesPopup(dialog, phases='GL')
    species_popup.show()
    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    qapp.exec_()
    qapp.deleteLater()

    sys.exit()
