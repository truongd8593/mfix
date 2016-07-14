#!/usr/bin/env python

"""Regions popup for MFIX-GUI
used for initial & boundary conditions"""

import os
import sys
import signal
from collections import OrderedDict
import pickle

from qtpy import QtCore, QtWidgets, PYQT5, uic
from qtpy.QtWidgets import QTableWidgetItem, QLineEdit
from qtpy.QtGui import QValidator, QDoubleValidator
UserRole = QtCore.Qt.UserRole

def set_item_noedit(item):
    item.setFlags(item.flags() ^ QtCore.Qt.ItemIsEditable)

def get_selected_row(table):
    # note, currentRow can return  >0 even when there is no selection
    rows = set(i.row() for i in table.selectedItems())
    return None if not rows else rows.pop()

if PYQT5:
    def resize_column(table, col, flags):
        table.horizontalHeader().setSectionResizeMode(col, flags)
else:
    def resize_column(table, col, flags):
        table.horizontalHeader().setResizeMode(col, flags)

class RegionsPopup(QtWidgets.QDialog):

    save = QtCore.Signal()
    cancel = QtCore.Signal()

    def handle_regions_selection(self):
        table = self.tablewidget_regions
        row = get_selected_row(table)

        if row is None:
            self.current_species = None
            self.clear_species_panel()
            self.pushbutton_delete.setEnabled(False)
            self.ui.combobox_phase.setEnabled(False)
        else:
            self.pushbutton_delete.setEnabled(True)
            self.current_species = table.item(row, 0).data(UserRole)
            self.enable_species_panel()

    def reset_signals(self):
        # todo:  fix this so it's not the caller's responsibility
        #  (make a util function that calls this & pops up dialog)
        for sig in (self.cancel, self.save):
            try:
                sig.disconnect()
            except:
                pass # isSignalConnected only exists in qt5.

    def __init__(self, app, parent=None):
        super(RegionsPopup, self).__init__(parent)
        self.app = app
        ui = self.ui = uic.loadUi(os.path.join(uidir, 'regions_popup.ui'), self)

        # key=species, val=data dict
        self.defined_regions = OrderedDict()

        buttons = ui.buttonbox.buttons()
        buttons[0].clicked.connect(lambda: self.save.emit())
        buttons[1].clicked.connect(lambda: self.cancel.emit())

    def popup(self):
        self.show()
        self.raise_()
        self.activateWindow()
