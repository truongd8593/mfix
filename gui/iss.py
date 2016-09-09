# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtWidgets import (QLabel, QLineEdit, QPushButton, QGridLayout,
                            QHBoxLayout, QWidget, QGroupBox)

UserRole = QtCore.Qt.UserRole

from constants import *
from widgets.regions_popup import RegionsPopup
from widgets.base import LineEdit, ComboBox

from project import Equation, FloatExp

from tools.general import (set_item_noedit, set_item_enabled,
                           get_combobox_item, get_selected_row,
                           widget_iter)

from tools.keyword_args import mkargs

from json import JSONDecoder, JSONEncoder

def safe_float(val):
    try:
        return float(val)
    except ValueError:
        return 0.0

class ISS(object):
    #Internal Surfaces Task Pane Window: This section allows a user to define internal surfaces for
    #the described model. This section relies on regions named in the Regions section.

    def init_iss(self):
        ui = self.ui.internal_surfaces

        self.iss = {} # key: index.  value: data dictionary for internal surface
        self.iss_current_indices = [] # List of IS indices
        self.iss_current_regions = [] # And the names of the regions which define them
        self.iss_region_dict = None

        #The top of the task pane is where users define/select IS regions
        #Icons to add/remove/duplicate boundary conditions are given at the top
        #Clicking the 'add' and 'duplicate' buttons triggers a popup window where the user must select
        #a region to apply the internal surface.
        ui.toolbutton_add.clicked.connect(self.iss_show_regions_popup)
        ui.toolbutton_delete.clicked.connect(self.iss_delete_regions)
        ui.toolbutton_delete.setEnabled(False) # Need a selection

        ui.tablewidget_regions.itemSelectionChanged.connect(self.handle_iss_region_selection)


    def iss_show_regions_popup(self):
        # Users cannot select inapplicable regions.
        # IS regions can be planes or volumes (not points or STLs)
        # No region can define more than one internal surface.

        #Select internal surface type
        # Selection is required
        # Available selections:
        #  Impermeable
        #    Selection only available for plane regions
        #    set keyword IS_TYPE(#) to 'IMPERMEABLE'
        #  X-Axis Impermeable
        #    Selection only available for volume regions
        #    set keyword IS_TYPE(#) to 'X_IMPERMEABLE'
        #  Y-Axis Impermeable
        #    Selection only available for volume regions
        #    set keyword IS_TYPE(#) to 'Y_IMPERMEABLE'
        #  z-Axis Impermeable
        #    Selection only available for volume regions
        #    set keyword IS_TYPE(#) to 'Z_IMPERMEABLE'
        #  Semi-permeable
        #    Selection only available for plane regions
        #    set keyword IS_TYPE(#) to 'SEMIPERMEABLE'
        #  X-Axis semi-permeable
        #    Selection only available for volume regions
        #    set keyword IS_TYPE(#) to 'X_SEMIPERMEABLE'
        #  Y-Axis semi-permeable
        #    Selection only available for volume regions
        #    set keyword IS_TYPE(#) to 'Y_SEMIPERMEABLE'
        #  Z-Axis semi-permeable
        #    Selection only available for volume regions
        #    set keyword IS_TYPE(#) to 'Z_SEMIPERMEABLE'
        # DEFAULT - Impermeable
        # (selection logic implemented in regions_popup.py)

        ui = self.ui.internal_surfaces
        rp = self.regions_popup
        rp.clear()
        for (name,data) in self.iss_region_dict.items():
            shape = data.get('type', '---')
            # Assume available if unmarked
            available = data.get('available', True) and (shape=='box' or 'plane' in shape)
            row = (name, shape, available)
            rp.add_row(row)
        rp.reset_signals()
        rp.save.connect(self.iss_add_regions)
        rp.cancel.connect(self.iss_cancel_add)
        for item in (ui.tablewidget_regions,
                     ui.detail_pane,
                     ui.toolbutton_add,
                     ui.toolbutton_delete):
            item.setEnabled(False)
        rp.popup('internal surface')


    def iss_cancel_add(self):
        ui = self.ui.internal_surfaces

        for item in (ui.toolbutton_add,
                     ui.tablewidget_regions):
            item.setEnabled(True)

        if get_selected_row(ui.tablewidget_regions) is not None:
            for item in (ui.detail_pane,
                         ui.toolbutton_delete):
                item.setEnabled(True)


    def iss_add_regions(self):
        # Interactively add regions to define ISs
        ui = self.ui.internal_surfaces
        rp = self.regions_popup
        self.iss_cancel_add() # Reenable input widgets
        selections = rp.get_selection_list()
        is_type = rp.combobox.currentIndex()
        if not selections:
            return
        self.iss_add_regions_1(selections, is_type) # Indices will be assigned
        self.iss_setup_current_tab() # Update the widgets


    def iss_add_regions_1(self, selections,
                          is_type=DEFAULT_IS_TYPE, indices=None):
        # Used by both interactive and load-time add-region handlers
        ui = self.ui.internal_surfaces

        if self.iss_region_dict is None:
            self.iss_region_dict = self.ui.regions.get_region_dict()

        tw = ui.tablewidget_regions
        nrows = tw.rowCount()
        tw.setRowCount(nrows+1)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        item = make_item('+'.join(selections))

        if indices is None: # interactive
            indices = [None] * len(selections)
            autoselect = True
        else: # loading file
            assert len(selections) == len(indices)
            autoselect = False

        for (i, region_name) in enumerate(selections):
            idx = indices[i]
            if idx is None:
                idx = self.iss_find_index()
                indices[i] = idx
            self.iss[idx] = {'region': region_name}
            region_data = self.iss_region_dict.get(region_name)
            if region_data is None: # ?
                self.warn("no data for region %s" % region_name)
                continue
            self.iss_set_region_keys(region_name, idx, region_data, is_type)
            self.iss_region_dict[region_name]['available'] = False # Mark as in-use

        item.setData(UserRole, (tuple(indices), tuple(selections)))
        tw.setItem(nrows, 0, item)

        item = make_item(IS_NAMES[is_type])
        tw.setItem(nrows, 1, item)

        self.fixup_iss_table(tw)

        if autoselect:
            tw.setCurrentCell(nrows, 0)


    def iss_find_index(self):
        n = 1
        while n in self.iss:
            n += 1
        return n


    def iss_delete_regions(self):
        ui = self.ui.internal_surfaces
        tw = ui.tablewidget_regions
        row = get_selected_row(tw)
        if row is None: # No selection
            return

        # Unset keywords
        kwlist = list(self.project.keywordItems())
        for kw in kwlist:
            key, args = kw.key, kw.args
            if key.startswith('is_') and args and args[0] in self.iss_current_indices:
                self.unset_keyword(key, args=args)

        # TODO: fix any resulting holes in index sequence!

        for r in self.iss_current_regions:
            if r in self.iss_region_dict:
                self.iss_region_dict[r]['available'] = True

        for i in self.iss_current_indices:
            del self.iss[i]

        self.iss_current_regions = []
        self.iss_current_indices = []

        tw.removeRow(row)
        self.fixup_iss_table(tw)
        self.iss_setup_current_tab()


    def handle_iss_region_selection(self):
        ui = self.ui.internal_surfaces
        table = ui.tablewidget_regions
        row = get_selected_row(table)
        if row is None:
            indices = []
            regions = []
        else:
            (indices, regions) = table.item(row,0).data(UserRole)
        self.iss_current_indices, self.iss_current_regions = indices, regions
        enabled = (row is not None)
        ui.toolbutton_delete.setEnabled(enabled)
        ui.detail_pane.setEnabled(enabled)
        if not enabled:
            # Clear
            for widget in widget_iter(ui.detail_pane):
                if isinstance(widget, LineEdit):
                    widget.setText('')
            return
        self.setup_iss() # reinitialize all widgets in current tab


    def fixup_iss_table(self, tw, stretch_column=0):
        ui = self.ui.internal_surfaces
        hv = QtWidgets.QHeaderView
        if PYQT5:
            resize = tw.horizontalHeader().setSectionResizeMode
        else:
            resize = tw.horizontalHeader().setResizeMode
        ncols = tw.columnCount()
        for n in range(0, ncols):
            resize(n, hv.Stretch if n==stretch_column else hv.ResizeToContents)

        # trim excess vertical space - can't figure out how to do this in designer
        header_height = tw.horizontalHeader().height()

        # TODO FIXME scrollbar handling is not right - scrollbar status can change
        # outside of this function.  We need to call this everytime window geometry changes
        scrollbar_height = tw.horizontalScrollBar().isVisible() * (4+tw.horizontalScrollBar().height())
        nrows = tw.rowCount()
        if nrows==0:
            height = header_height+scrollbar_height
        else:
            height =  (header_height+scrollbar_height
                       + nrows*tw.rowHeight(0) + 4) # extra to avoid unneeded scrollbar

        if tw == ui.tablewidget_regions: # main table, adjust top splitter
            ui.top_frame.setMaximumHeight(height+40)
            ui.top_frame.setMinimumHeight(header_height+40)
            ui.top_frame.updateGeometry()
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(header_height)
        else: # mass fraction tables
            tw.setMaximumHeight(height) # Works for tablewidget inside groupbox
            tw.setMinimumHeight(height) #? needed? should we allow scrollbar?
        tw.updateGeometry() #? needed?

    def iss_set_region_keys(self, name, idx, data, is_type=None):
        # Update the keys which define the region the IS applies to
        if is_type is not None:
            val = IS_TYPES[is_type]
            self.update_keyword('is_type', val, args=[idx])

        no_k = self.project.get_value('no_k')
        for (key, val) in zip(('x_w', 'y_s', 'z_b',
                               'x_e', 'y_n', 'z_t'),
                              data['from']+data['to']):
            # is_z_t and is_z_b keywords should not be added when no_k=True
            if no_k and key in ('z_t', 'z_b'):
                continue
            self.update_keyword('is_'+key, val, args=[idx])


    def iss_change_region_name(self, old_name, new_name):
        ui = self.ui.internal_surfaces
        for (key, val) in self.iss.items():
            if val.get('region') == old_name:
                self.iss[key]['region'] = new_name
                tw = ui.tablewidget_regions
                for i in range(tw.rowCount()):
                    data = tw.item(i,0).data(UserRole)
                    indices, names = data
                    if key in indices:
                        item = tw.item(i,0)
                        names = (new_name if n==old_name else n for n in names)
                        item.setData(UserRole, (indices, names))
                        item.setText('+'.join(names))
                        break
                break


    def reset_iss(self):
        self.iss.clear()
        self.iss_current_indices = []
        self.iss_current_regions = []
        self.iss_region_dict = None
        self.iss_current_solid = self.P = None
        ui = self.ui.internal_surfaces
        ui.tablewidget_regions.clearContents()
        ui.tablewidget_regions.setRowCount(0)
        # anything else to do here?


    def iss_to_str(self):
        ui = self.ui.internal_surfaces
        tw = ui.tablewidget_regions
        data = [tw.item(i,0).data(UserRole)
                for i in range(tw.rowCount())]
        return JSONEncoder().encode(data)


    def iss_regions_from_str(self, s):
        if not s:
            return
        data = JSONDecoder().decode(s)
        for (indices, regions) in data:
            self.iss_add_regions_1(regions, indices)


    def iss_extract_regions(self):
        if self.iss:
            # We assume that ic regions have been initialized correctly
            # from mfix_gui_comments.
            # TODO: verify that there is an IS region for each IS
            return

        if self.iss_region_dict is None:
            self.iss_region_dict = self.ui.regions.get_region_dict()

        for is_ in self.project.iss:
            d = is_.keyword_dict
            extent = [d.get('is_'+k,None) for k in ('x_w', 'y_s', 'z_b',
                                                    'x_e', 'y_n', 'z_t')]
            extent = [0 if x is None else x.value for x in extent]
            #if any (x is None for x in extent):
            #    self.warn("internal surface %s: invalid extents %s" %
            #               (is_.ind, extent))
            #    continue
            for (region_name, data) in self.iss_region_dict.items():
                ext2 = [0 if x is None else x for x in
                        (data.get('from',[]) + data.get('to',[]))]
                if ext2 == extent:
                    if data.get('available', True):
                        self.iss_add_regions_1([region_name], [is_.ind])
                        break
            else:
                self.warn("internal surface %s: could not match defined region %s" %
                          (is_.ind, extent))


    def setup_iss(self):
        ui = self.ui.internal_surfaces

        # Grab a fresh copy, may have been updated
        self.iss_region_dict = self.ui.regions.get_region_dict()

        # Mark regions which are in use (this gets reset each time we get here)
        for (i, data) in self.iss.items(): # Should we also consider BCs/ICs?
            region = data['region']
            if region in self.iss_region_dict:
                self.iss_region_dict[region]['available'] = False

        self.fixup_iss_table(ui.tablewidget_regions)
        row = get_selected_row(ui.tablewidget_regions)
        enabled = (row is not None)
        ui.toolbutton_delete.setEnabled(enabled)
        ui.detail_pane.setEnabled(enabled)
        self.iss_setup_current_tab()


    def iss_setup_current_tab(self):
        #Input is only needed for semi-permeable surfaces.
        #Gas permeability:
        # Specification only available for semipermeable regions
        # DEFAULT value 1.0d32
        # Sets keyword IS_PC(#,1)
        #Internal resistance coefficient:
        # Specification only available for semipermeable regions
        # DEFAULT value 0.0
        # Sets keyword IS_PC(#,2)
        #Solids velocity through surface:
        # Specification only available for semipermeable regions
        # DEFAULT value 0.0
        # Sets keyword IS_VEL_s(#,PHASE)
        # There should be a line for each solids phase. Use the user provided solids name.
