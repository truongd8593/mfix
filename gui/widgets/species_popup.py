#!/usr/bin/env python

"""Species selector dialog for MFIX GUI, includes stand-alone test"""

import os
import sys
import signal
import time
import cPickle

from qtpy import QtCore, QtWidgets, QtGui
from qtpy.QtWidgets import QTableWidgetItem, QLineEdit
from qtpy.QtGui import QValidator, QDoubleValidator
UserRole = QtCore.Qt.UserRole

try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic


def set_item_noedit(item):
    item.setFlags(item.flags() ^ QtCore.Qt.ItemIsEditable)

class SpeciesPopup(QtWidgets.QDialog):

    save = QtCore.Signal()
    cancel = QtCore.Signal()

    def load_burcat(self, path):
        if not os.path.exists(path):
            print("%s not found, create it by running read_burcat.py" % path)
            sys.exit(-1)
        with open(path) as f:
            db = cPickle.load(f)
        by_phase = {}

        for k,v in db.items():
            phase = k[1]
            if phase not in by_phase:
                by_phase[phase] = {}
            key = k[0], k[2], k[3]
            by_phase[phase][key] = v
        self.db = by_phase

        def format_comment(comment):
            lines = []
            while comment:
                line, comment = comment[:80], comment[80:]
                lines.append(line)
            return '\n'.join(lines)

        # build search list, lowercased
        self.haystack = []
        self.comments = {}
        for phase in 'GLCS':
            # TODO - also match comments?
            htmp = [(k[0].lower(), k, phase) for k in self.db[phase].keys()]
            htmp.sort()
            self.haystack.extend(htmp)
            # comment fields
            self.comments[phase] = dict((k, format_comment(v[2]))
                                        for (k,v) in self.db[phase].items())

    def do_search(self, string):
        results = {}
        self.ui.tablewidget_search.clearContents()
        needle = string.lower()
        results = []
        for ((key_low, key, phase)) in self.haystack:
            if needle in key_low and phase in self.phases:
                results.append((key, phase))
        table = self.ui.tablewidget_search
        nrows = len(results)
        self.ui.tablewidget_search.setRowCount(nrows)
        self.search_results = [None]*nrows
        for (i,r) in enumerate(results):
            key, phase = r
            comment = self.comments[phase][key]
            item = QTableWidgetItem(key[0])
            item.setToolTip(comment)
            table.setItem(i, 0, item)
            item = QTableWidgetItem(phase)
            table.setItem(i, 1, item)
            self.search_results[i] = (key, phase)

    def handle_search_selection(self):
        table = self.tablewidget_search
        rows = set(i.row() for i in (table.selectedItems()))
        if len(rows) == 1: # Multiple selection is disabled
            self.selected_search_row = rows.pop()
        else:
            self.selected_search_row = None

        if self.selected_search_row is None:
            self.ui.pushbutton_import.setEnabled(False)
        else:
            self.ui.pushbutton_import.setEnabled(True)

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
                item.textEdited.connect(make_handler(key=key))
            return item

        def make_handler(key):
            def handler(val, key=key):
                if not self.current_species:
                    print("Error, no current species")
                try:
                    data = self.defined_species[self.current_species]
                    val = float(val)
                    if type(key) is tuple:
                        data[key[0]][key[1]] = val
                    else:
                        data[key] = val
                    data['source'] = 'User Defined' # Data has been modified
                    ui.label_species_source.setText(data['source'])
                except ValueError:
                    # reset field to prev. value
                    pass

            return handler
        ui.label_species_source.setText(data['source'])
        ui.label_species.setText(species)
        ui.lineedit_alias.setText(data['alias'])
        ui.lineedit_molecular_weight.setText(str(data['molecular_weight']))
        ui.lineedit_molecular_weight.textEdited.connect(make_handler('molecular_weight'))
        ui.lineedit_heat_of_formation.setText(str(data['heat_of_formation']))
        ui.lineedit_heat_of_formation.textEdited.connect(make_handler('heat_of_formation'))
        # Density - disabled

        table = ui.tablewidget_params
        table.setCellWidget(0, 0, make_item(data['tmin'], key='tmin'))
        table.setCellWidget(0, 1, make_item(data['tmax'], key='tmax'))
        for (i,x) in enumerate(data['a_low'], 1):
            table.setCellWidget(i, 0, make_item(x, key=('a_low', i)))
        for (i,x) in enumerate(data['a_high'], 1):
            table.setCellWidget(i, 1, make_item(x, key=('a_high', i)))

    def handle_defined_species_selection(self):
        self.ui.tablewidget_search.clearSelection() # is this right?
        table = self.tablewidget_defined_species
        rows = set(d.row() for d in table.selectedItems())
        species = set(d.data(UserRole) for d in table.selectedItems())
        if None in species:
            species.remove(None)

        if len(rows) != 1:
            self.current_species = None
            self.selected_species_row = None
            self.clear_species_panel()
            self.pushbutton_delete.setEnabled(False)
        else:
            self.pushbutton_delete.setEnabled(True)
            self.current_species = species.pop()
            self.selected_species_row = rows.pop()
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
        rowdata = self.search_results[self.selected_search_row]
        key, phase = rowdata
        data = self.db[phase][key]
        (species, tmin, tmax) = key
        (coeffs, molecular_weight, comment) = data
        # Do we want to allow same species in different phases, etc?
        if species in self.defined_species:
            return # Don't allow duplicates

        alias = self.make_alias(species)
        a_low = coeffs[:7]
        a_high = coeffs[7:14]
        heat_of_formation = coeffs[14]

        species_data = {'source': 'BURCAT',
                        'phase': phase,
                        'alias': alias,
                        'molecular_weight': molecular_weight,
                        'heat_of_formation': heat_of_formation,
                        'tmin':  tmin,
                        'tmax': tmax,
                        'a_low': a_low,
                        'a_high': a_high}

        self.defined_species[species] = species_data
        self.add_defined_species_row(species)
        self.set_save_button(True)

    def add_defined_species_row(self, species):
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

        self.selected_species_row = nrows
        table.setCurrentCell(nrows, 0) # Cause the new row to be selected

    def handle_delete(self):
        table = self.ui.tablewidget_defined_species
        rows = set(item.row() for item in table.selectedItems())
        if len(rows)!=1: # multiple selection disabled
            return
        row = rows.pop()
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

    def handle_new(self):
        phase = self.default_phase # Note - no way for user to specify phase - FIXME
        species = self.make_user_species_name()
        alias = species.replace(' ', '') #?
        molecular_weight = 0
        heat_of_formation = 0
        tmin = 200.0 # ?
        tmax = 600.0 # ?
        a_low = [0.0]*7
        a_high = [0.0]*7

        species_data = {'source': 'User Defined',
                        'phase': phase,
                        'alias': alias,
                        'molecular_weight': molecular_weight,
                        'heat_of_formation': heat_of_formation,
                        'tmin':  tmin,
                        'tmax': tmax,
                        'a_low': a_low,
                        'a_high': a_high}

        self.defined_species[species] = species_data
        self.current_species = species
        self.enable_species_panel()
        self.add_defined_species_row(species)
        lineedit = self.ui.lineedit_alias
        lineedit.selectAll()
        lineedit.setFocus()

    def handle_alias(self, val):
        ui = self.ui
        row = self.selected_species_row
        if row is None:
            return
        item = QTableWidgetItem(val)
        item.setData(UserRole, self.current_species)
        ui.tablewidget_defined_species.setItem(row, 0, item)
        self.defined_species[self.current_species]['alias'] = val

    def set_save_button(self, state):
        self.ui.buttonbox.buttons()[0].setEnabled(state)

    def __init__(self, app, parent=None, phases='GLCS'):
        super(SpeciesPopup, self).__init__(parent)
        self.app = app
        self.phases = phases
        self.default_phase = phases[0]
        thisdir = os.path.abspath(os.path.dirname(__file__))
        datadir = thisdir
        self.load_burcat(os.path.join(datadir, 'burcat.pickle'))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        ui = self.ui = uic.loadUi(os.path.join(uidir, 'species_popup.ui'), self)

        self.defined_species = {} # key=species, val=data tuple.  can add phase to key if needed
        self.search_results = []
        self.user_species_names = set()
        self.selected_search_row = None
        self.selected_species_row = None

        # Set up UI
        ui.lineedit_search.textChanged.connect(self.do_search)
        self.do_search('')
        ui.pushbutton_import.clicked.connect(self.do_import)
        ui.pushbutton_import.setEnabled(False)
        ui.tablewidget_search.itemSelectionChanged.connect(
            self.handle_search_selection)
        ui.tablewidget_defined_species.itemSelectionChanged.connect(
            self.handle_defined_species_selection)

        ui.pushbutton_new.clicked.connect(self.handle_new)
        ui.pushbutton_delete.clicked.connect(self.handle_delete)

        buttons = ui.buttonbox.buttons()
        buttons[0].clicked.connect(lambda: self.save.emit())
        buttons[1].clicked.connect(lambda: self.cancel.emit())

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
                self.parent.set_save_button(True)
                return (QValidator.Acceptable, text, pos)

        lineedit = ui.lineedit_alias
        lineedit.setValidator(AliasValidator(parent=self))
        lineedit.textEdited.connect(self.handle_alias)

        for l in (ui.lineedit_molecular_weight,
                  ui.lineedit_heat_of_formation,
                  ui.lineedit_density):
            l.setValidator(QDoubleValidator())

        self.species_panel_items=[
            ui.label_species_source,
            ui.label_species,
            ui.lineedit_alias,
            ui.lineedit_molecular_weight,
            ui.lineedit_heat_of_formation,
            ui.combobox_specific_heat_model,
            ui.tablewidget_params]

        self.set_save_button(False) # nothing to Save
        self.clear_species_panel()


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
