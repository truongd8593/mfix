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

from tools.general import (set_item_noedit, set_item_enabled,
                           get_combobox_item, get_selected_rows)

if PYQT5:
    def resize_column(table, col, flags):
        table.horizontalHeader().setSectionResizeMode(col, flags)
else:
    def resize_column(table, col, flags):
        table.horizontalHeader().setResizeMode(col, flags)

# move to "constants"
# must match order in combobox_bc_type
bc_types = ['MI', 'PO', 'NSW', 'FSW', 'PSW', 'PI', 'MO']
default_bc_type = bc_types.index('NSW')  # No Slip Wall

class RegionsPopup(QtWidgets.QDialog):

    save = QtCore.Signal()
    cancel = QtCore.Signal()

    def handle_regions_selection(self):
        tw = self.ui.table
        cb = self.ui.combobox_bc_type

        if self.boundary:
            #  Pressure Inflow: Not available for STL regions
            pi_index = bc_types.index('PI')
            pi_item = get_combobox_item(cb, pi_index)
            disable = any(tw.item(x,1).text()=='STL'
                          for x in get_selected_rows(tw))
            set_item_enabled(pi_item, not disable)

            # Don't stay on disabled item
            if disable and cb.currentIndex() == pi_index:
                cb.setCurrentIndex(default_bc_type)
        self.update_available_regions()


    def handle_bc_type(self, val):
        self.update_available_regions()


    def update_available_regions(self):
        tw = self.ui.table
        cb = self.ui.combobox_bc_type
        bc_type = bc_types[cb.currentIndex()]
        selections = get_selected_rows(tw)
        if self.boundary:
            # Wall type boundary
            if bc_type.endswith('W'):
                self.reset_available()

            # For inflows/outflows, only allow compatible orientation
            elif bc_type.endswith('I') or bc_type.endswith('O'):
                if len(selections) == 1:
                    region_type = tw.item(selections[0],1).text()
                    for i in range(0, tw.rowCount()):
                        if i == selections[0]:
                            continue
                        enable = (tw.item(i,1).text() == region_type)
                        for j in (0,1,2):
                            set_item_enabled(tw.item(i,j), enable)
                            tw.item(i,2).setText('Yes' if enable else 'No')
                elif len(selections) == 0:
                    self.reset_available()

    def reset_available(self):
        tw = self.ui.table
        for i in range(0, tw.rowCount()):
            enable = tw.item(i,2).data(UserRole)
            for j in (0,1,2):
                set_item_enabled(tw.item(i,j), enable)
                tw.item(i,2).setText('Yes' if enable else 'No')

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
        #self.ui.table.doubleClicked.connect(self.save.emit) # Too easy to double-click accidentally
        #self.ui.table.doubleClicked.connect(self.accept)
        self.ui.table.itemSelectionChanged.connect(self.handle_regions_selection)
        self.ui.combobox_bc_type.currentIndexChanged.connect(self.handle_bc_type)
        self.rejected.connect(self.cancel.emit)


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
        # Keep track of original 'available' so we can restore
        table.item(nrows, 2).setData(UserRole, available)


    def popup(self, boundary=False):
        # Widget is shared by ICs and BCs pane
        text = "Select region(s) for %s coundition" % ('boundary' if boundary else 'initial')
        self.boundary = boundary
        for item in (self.ui.combobox_bc_type, self.ui.label_bc_type):
            if boundary:
                item.show()
            else:
                item.hide()
        self.ui.combobox_bc_type.setCurrentIndex(default_bc_type)
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
