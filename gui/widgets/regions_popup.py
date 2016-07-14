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
        pass

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
        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        ui = self.ui = uic.loadUi(os.path.join(uidir, 'regions_popup.ui'), self)

        # key=species, val=data dict
        self.defined_regions = OrderedDict()

        buttons = ui.buttonbox.buttons()
        buttons[0].clicked.connect(lambda: self.save.emit())
        buttons[1].clicked.connect(lambda: self.cancel.emit())


    def clear(self):
        self.ui.table.clear()

    def add_row(self, row):
        table = self.ui.table
        row_count = table.rowCount()
        table.setRowCount(row_count+1)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        table.setItem(row_count, 0, make_item(row))

    def popup(self):
        self.show()
        self.raise_()
        self.activateWindow()
