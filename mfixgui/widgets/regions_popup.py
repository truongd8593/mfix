"""Regions popup for MFIX-GUI
used for initial & boundary conditions"""

# Note, throughout this module, the abbreviation 'is_' means "internal surface", not "is"

from __future__ import print_function, absolute_import, unicode_literals, division

import os
from collections import OrderedDict

from qtpy import uic

from qtpy.QtCore import (
    Qt,
    Signal,
)

from qtpy.QtWidgets import (
    QHeaderView,
    QTableWidgetItem,
    QAbstractItemView,
    QLineEdit,
    QTableWidgetItem,
    QDialog,
)

from mfixgui.constants import (
    IMPERMEABLE,
    IS_NAMES,
    BC_TYPES,
    BC_NAMES,
    DEFAULT_IS_TYPE,
    DEFAULT_BC_TYPE,
    IS_TYPES,
    MASS_INFLOW,
    MASS_OUTFLOW,
    PRESSURE_INFLOW,
    PRESSURE_OUTFLOW,
    SEMIPERMEABLE,
    X_IMPERMEABLE,
    X_SEMIPERMEABLE,
    Y_IMPERMEABLE,
    Y_SEMIPERMEABLE,
    Z_IMPERMEABLE,
    Z_SEMIPERMEABLE,
)

from mfixgui.tools.general import (
    get_combobox_item,
    get_selected_rows,
    item_enabled,
    set_item_enabled,
    set_item_noedit,
)

USER_ROLE = Qt.UserRole

def resize_column(table, col, flags):
    table.horizontalHeader().setSectionResizeMode(col, flags)


class RegionsPopup(QDialog):

    save = Signal()
    cancel = Signal()

    def handle_regions_selection(self):
        tw = self.ui.table
        cb = self.ui.combobox
        selections = get_selected_rows(tw)

        buttonbox = self.ui.buttonbox
        buttonbox.button(buttonbox.Ok).setEnabled(bool(selections))

        region_types = [tw.item(x, 1).text() for x in selections]
        types_match = len(set(region_types)) < 2

        if self.boundary:
            #  Pressure Inflow
            #     Not available for STL regions
            #     Not available for volume regions
            #  (also, chained regions must match type)

            disable = any(x in ('STL', 'box') for x in region_types) or not types_match
            item = get_combobox_item(cb, PRESSURE_INFLOW)
            set_item_enabled(item, not disable)

            # Volume BCs are only allowed for walls
            disable = any(x == 'box' for x in region_types) or not types_match
            for idx in (MASS_INFLOW, PRESSURE_OUTFLOW, MASS_OUTFLOW):
                item = get_combobox_item(cb, idx)
                set_item_enabled(item, not disable)

            # Don't stay on disabled item
            if not item_enabled(get_combobox_item(cb, cb.currentIndex())):
                cb.setCurrentIndex(BC_TYPES.index(DEFAULT_BC_TYPE))

            bc_type = BC_TYPES[cb.currentIndex()]
            # Wall type boundary
            if bc_type.endswith('W'):
                self.handle_type(self.ui.combobox.currentIndex())
            # For inflows/outflows, only allow compatible orientation
            else:
                if len(selections) == 1:
                    region_type = tw.item(selections[0], 1).text()
                    for i in range(0, tw.rowCount()):
                        if i == selections[0]:
                            continue
                        enable = (tw.item(i, 1).text() == region_type)
                        self.enable_row(i, enable)
                elif len(selections) == 0:
                    self.handle_type(self.ui.combobox.currentIndex())

                else:
                    pass


        elif self.surface:
            #  (would be nicer if this were in iss.py)
            # IS regions can be planes or volumes (not points or STLs)
            # *-Axis permeable/semipermeable not available for planes
            disable = selections and any('plane' in tw.item(x, 1).text()
                                         for x in selections)
            for index in (X_SEMIPERMEABLE, Y_SEMIPERMEABLE, Z_SEMIPERMEABLE,
                          X_IMPERMEABLE, Y_IMPERMEABLE, Z_IMPERMEABLE):
                item = get_combobox_item(cb, index)
                set_item_enabled(item, (not selections) or (not disable))
            for index in (SEMIPERMEABLE, IMPERMEABLE):
                item = get_combobox_item(cb, index)
                set_item_enabled(item, (not selections) or disable)

            # Available selections:
            #  Impermeable
            #    Selection only available for plane regions
            #  X-Axis Impermeable
            #    Selection only available for volume regions
            #  Y-Axis Impermeable
            #    Selection only available for volume regions
            #  z-Axis Impermeable
            #    Selection only available for volume regions

            #  Semi-permeable
            #    Selection only available for plane regions
            #  X-Axis semi-permeable
            #    Selection only available for volume regions
            #  Y-Axis semi-permeable
            #    Selection only available for volume regions
            #  Z-Axis semi-permeable
            #    Selection only available for volume regions
            # DEFAULT - Impermeable

            is_type = IS_TYPES[cb.currentIndex()]

            if len(selections) == 1:
                region_type = tw.item(selections[0], 1).text()
                plane = 'plane' in region_type
                for i in range(0, tw.rowCount()):
                    if i == selections[0]:
                        continue
                    text = tw.item(i, 1).text()
                    enable = (('plane' in text or text == 'box')
                              or (plane == ('plane' in tw.item(i, 1).text())))
                    self.enable_row(i, enable)
            elif len(selections) == 0:
                self.handle_type(self.ui.combobox.currentIndex())


    def handle_type(self, val):
        tw = self.ui.table
        cb = self.ui.combobox
        selections = get_selected_rows(tw)
        target = tw.item(selections[0], 1).text() if selections else None

        if self.boundary:
            bc_type = BC_TYPES[val]
            if bc_type.endswith('W'):
                for i in range(0, tw.rowCount()):
                    text = tw.item(i, 1).text()
                    enable = 'cyclic' not in text
                    self.enable_row(i, enable)

            elif bc_type == 'PI':
                #    Not available for volume regions
                #    Not available for STL regions
                for i in range(tw.rowCount()):
                    text = tw.item(i, 1).text()
                    enable = ('plane' in text) and (target is None or text == target)
                    self.enable_row(i, enable)

            elif bc_type in ('MI', 'PO', 'MO'):
                #    Not available for volume regions
                for i in range(tw.rowCount()):
                    text = tw.item(i, 1).text()
                    enable = (('plane' in text or text == 'STL')
                              and (target is None or text==target))
                    self.enable_row(i, enable)
            elif bc_type == 'CYCLIC':
                tw.clearSelection()
                for i in range(tw.rowCount()):
                    text = tw.item(i, 1).text()
                    enable = 'cyclic' in text
                    self.enable_row(i, enable)
            else:
                self.error("Unknown bc_type %s" % bc_type)
            tw.setSelectionMode(QAbstractItemView.SingleSelection if bc_type == 'CYCLIC'
                                else QAbstractItemView.MultiSelection)

        elif self.surface:
            tw.setSelectionMode(QAbstractItemView.MultiSelection)
            is_type = IS_TYPES[val]
            plane = is_type in ('IMPERMEABLE', 'SEMIPERMEABLE')
            #enable planes, disable boxes
            for i in range(tw.rowCount()):
                text = tw.item(i, 1).text()
                enable = ((text == 'box' or 'plane' in text)
                          and (('plane' in tw.item(i, 1).text()) == plane))
                self.enable_row(i, enable)


    def enable_row(self, i, enable):
        tw = self.ui.table
        enable = enable and tw.item(i, 2).data(USER_ROLE) # Initial enable setting
        tw.item(i, 2).setText('Yes' if enable else 'No')
        for j in (0, 1, 2):
            set_item_enabled(tw.item(i, j), enable)


    def reset_signals(self):
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
        self.ui = uic.loadUi(os.path.join(uidir, 'regions_popup.ui'), self)

        # key=species, val=data dict
        self.defined_regions = OrderedDict()

        buttons = self.ui.buttonbox.buttons()
        buttons[0].clicked.connect(lambda: self.save.emit())
        buttons[1].clicked.connect(lambda: self.cancel.emit())
        #self.ui.table.doubleClicked.connect(self.save.emit) # Too easy to double-click accidentally
        #self.ui.table.doubleClicked.connect(self.accept)
        self.ui.table.itemSelectionChanged.connect(self.handle_regions_selection)
        self.ui.combobox.currentIndexChanged.connect(self.handle_type)
        self.rejected.connect(self.cancel.emit)


    def clear(self):
        self.ui.table.clearContents()
        self.ui.table.setRowCount(0)


    def add_row(self, row):
        table = self.ui.table
        nrows = table.rowCount()
        table.setRowCount(nrows+1)
        def make_item(val, enabled):
            item = QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            set_item_enabled(item, enabled)
            return item
        (name, shape, available) = row
        table.setItem(nrows, 0, make_item(name, available))
        table.setItem(nrows, 1, make_item(shape, available))
        table.setItem(nrows, 2, make_item('Yes' if available else 'No', available))
        # Keep track of original 'available' so we can restore
        table.item(nrows, 2).setData(USER_ROLE, available)


    def popup(self, label_text):
        # Widget is shared by ICs/BCs/PSs/ISs/VTK output
        # NB, we distinguish the caller based on the label text
        ui = self.ui
        ui.label_top.setText(label_text)

        buttonbox = self.ui.buttonbox
        buttonbox.button(buttonbox.Ok).setEnabled(False)

        self.boundary = boundary = ('boundary condition' in label_text)
        self.surface = surface = ('internal surface' in label_text)
        self.vtk = vtk = ('VTK output' in label_text)

        tw = ui.table
        if vtk: # Issues/264 - don't combine vtk output regions,
                # because filenames will collide
            tw.setSelectionMode(QAbstractItemView.SingleSelection)
        else:
            tw.setSelectionMode(QAbstractItemView.MultiSelection)

        # setup the combobox appropriately
        cb = ui.combobox
        label = ui.label_type
        if boundary or surface:
            label.setText('Boundary type' if boundary else 'Surface type')
            cb.clear()
            cb.addItems(BC_NAMES if boundary else IS_NAMES)
            index = BC_TYPES.index(DEFAULT_BC_TYPE) if boundary else IS_TYPES.index(DEFAULT_IS_TYPE)
            cb.setCurrentIndex(index)
            self.handle_type(index)
            ui.frame_object_type.show()
        elif vtk:
            label.setText('Output type')
            cb.clear()
            cb.addItems(['Cell data', 'Particle data'])
            # slight hack to set tooltips
            get_combobox_item(cb, 0).setToolTip('Cell data (VTU file)')
            get_combobox_item(cb, 1).setToolTip('Particle data (VTP file)')
            ui.frame_object_type.show()
        else:
            ui.frame_object_type.hide()

        self.show()
        self.raise_()
        self.activateWindow()
        # Table fixup

        resize = tw.horizontalHeader().setSectionResizeMode
        ncols = tw.columnCount()
        stretch_column = 0
        for col_index in range(0, ncols):
            if col_index == stretch_column:
                resize(col_index, QHeaderView.Stretch)
            else:
                resize(col_index, QHeaderView.ResizeToContents)

        table_height = 4 + tw.horizontalHeader().height() + (tw.rowHeight(0)*tw.rowCount())

        min_table_height = 150
        tw.setMinimumHeight(min(min_table_height, table_height))
        tw.setMaximumHeight(table_height)

        width = (tw.horizontalHeader().width()
                 + tw.verticalScrollBar().isVisible() * (4+tw.verticalScrollBar().width()))
        self.setMaximumWidth(max(width, ui.label_top.size().width()))
        #self.setMinimumWidth(width)

        # Height of all "stuff" other than the main table
        stuff_height = 3 * self.combobox.size().height() + 20 # hack, size() of label returns weird value

        height = (table_height + stuff_height)

        self.setMinimumHeight(min(min_table_height+stuff_height, height))
        self.setMaximumHeight(height)


    def get_selection_list(self):
        """return list of selected region names"""
        rows = list(set([i.row() for i in self.ui.table.selectedItems()]))
        rows.sort()
        names = [self.ui.table.item(r, 0).text() for r in rows]
        return names
