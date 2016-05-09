#!/usr/bin/env python

"""Stand-alone test/demo for species selector dialog in MFiX GUI"""

import os
import sys
import signal
import time
import cPickle

from qtpy import QtCore, QtWidgets, QtGui
from qtpy.QtWidgets import QTableWidgetItem

try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic


class SpeciesPopup(QtWidgets.QDialog):

    def load_burcat(self, path):
        if not os.path.exists(path):
            print("%s not found, create it by running read_burcat.py")
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
            if needle in key_low:
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

    def handle_tablewidget_search_selection(self):
        table = self.tablewidget_search
        selected_rows = set(i.row() for i in (table.selectedIndexes()))
        if len(selected_rows) == 1: # Multiple selection is disabled
            self.selected_row = selected_rows.pop()
        else:
            self.selected_row = None

        if self.selected_row is None:
            self.ui.pushbutton_import.setEnabled(False)
        else:
            self.ui.pushbutton_import.setEnabled(True)


    def mkalias(self, species):
        return species.split(' ', 1)[0]

    def do_import(self):
        rowdata = self.search_results[self.selected_row]
        key, phase = rowdata
        data = self.db[phase][key]
        (species, tmin, tmax) = key
        (coeffs, mol_weight, comment) = data
        # Do we want to allow same species in different phases, etc?
        if speciesmk in self.defined_species:
            return # Don't allow duplicates
        self.defined_species.add(species)

        table = self.ui.tablewidget_defined_species
        nrows = table.rowCount()
        table.setRowCount(nrows+1)
        alias = self.mkalias(species)
        item = QTableWidgetItem(alias)
        table.setItem(nrows, 0, item)
        item = QTableWidgetItem(phase)
        table.setItem(nrows, 1, item)


    def __init__(self, app, parent=None):
        super(SpeciesPopup, self).__init__(parent)
        self.app = app

        self.load_burcat('tools/burcat.pickle')
        self.ui = uic.loadUi('uifiles/species_popup.ui', self)

        self.defined_species = set()
        self.search_results = []

        # Set up UI
        ui = self.ui
        ui.lineedit_search.textChanged.connect(self.do_search)
        ui.pushbutton_import.clicked.connect(self.do_import)
        ui.pushbutton_import.setEnabled(False)
        ui.tablewidget_search.itemSelectionChanged.connect(
            self.handle_tablewidget_search_selection)

if __name__ == '__main__':
    args = sys.argv
    qapp = QtWidgets.QApplication(args)
    dialog = QtWidgets.QDialog()
    species_popup = SpeciesPopup(dialog)
    species_popup.show()
    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    qapp.exec_()
    qapp.deleteLater()

    sys.exit()
