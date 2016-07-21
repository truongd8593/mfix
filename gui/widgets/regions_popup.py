from __future__ import print_function, absolute_import, unicode_literals, division

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

from tools.general import (set_item_noedit, set_item_enabled, get_selected_row)

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
        self.ui.table.clearContents()
        self.ui.table.setRowCount(0)

    def add_row(self, row):
        table = self.ui.table
        nrows = table.rowCount()
        table.setRowCount(nrows+1)
        def make_item(val, enabled):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            set_item_enabled(item, enabled)
            return item
        (name, shape, available) = row
        table.setItem(nrows, 0, make_item(name, available))
        table.setItem(nrows, 1, make_item(shape, available))
        table.setItem(nrows, 2, make_item('Yes' if available else 'No', available))

    def popup(self, text="Select regions"):
        self.ui.label.setText(text)
        self.show()
        self.raise_()
        self.activateWindow()

    def get_selection_list(self):
        """return list of selected region names"""
        rows = set([i.row() for i in self.ui.table.selectedItems()])
        rows = list(rows)
        rows.sort()
        names = [self.ui.table.item(r,0).text() for r in rows]
        return names
