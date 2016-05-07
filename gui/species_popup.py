#!/usr/bin/env python

"""Stand-alone test/demo for species selector dialog in MFiX GUI"""

import os
import sys
import signal
import time
import cPickle

from qtpy import QtCore, QtWidgets, QtGui

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
        # todo: move this into preprocessor
        for k,v in db.items():
            phase = k[1]
            if phase not in by_phase:
                by_phase[phase] = {}
            key = k[0], k[2], k[3]
            by_phase[phase][key] = v
        self.db = by_phase

    def do_search(self, string):
        results = {}
        self.ui.tablewidget_search.clearContents()
        string = string.lower()
        for phase in 'GSLC':
            results[phase] = []
            db = self.db[phase]
            for (k,v) in db.items():
                species = k[0]
                if string in species.lower():
                    # TODO - also match comments, options for search params
                    results[phase].append(species)

        table = self.ui.tablewidget_search
        nrows = sum(len(r) for r in results.values())
        self.ui.tablewidget_search.setRowCount(nrows)
        i = 0
        for phase in 'GSLC':
            results[phase].sort()
            for r in results[phase]:
                table.setItem(i, 0, QtWidgets.QTableWidgetItem(r))
                table.setItem(i, 1, QtWidgets.QTableWidgetItem(phase))
                i += 1

    def __init__(self, app, parent=None):
        super(SpeciesPopup, self).__init__(parent)
        self.app = app

        self.load_burcat('tools/burcat.pickle')
        self.ui = uic.loadUi('uifiles/species_popup.ui', self)

        ui = self.ui
        ui.lineedit_search.textChanged.connect(self.do_search)


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
