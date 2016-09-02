# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtWidgets import (QLabel, QLineEdit, QPushButton, QGridLayout,
                            QHBoxLayout, QWidget, QGroupBox)

from qtpy.QtGui import QPixmap # QPicture doesn't work with Qt4

UserRole = QtCore.Qt.UserRole

from constants import *
from widgets.regions_popup import RegionsPopup
from widgets.base import LineEdit, ComboBox

from project import Equation, FloatExp

from tools.general import (set_item_noedit, set_item_enabled,
                           get_combobox_item, get_selected_row,
                           widget_iter)

# We don't need extended JSON here
from json import JSONDecoder, JSONEncoder

def safe_float(val):
    try:
        return float(val)
    except ValueError:
        return 0.0

FLUID_TAB = 0
SOLIDS_TAB_DUMMY_L = 1
SOLIDS_TAB = 2
SOLIDS_TAB_DUMMY_R = 3
#SCALAR_TAB = 4

class PSS(object):
    #Point Source Task Pane Window: This section allows a user to define point sources for the
    #described model. This section relies on regions named in the Regions section.

    def init_pss(self):
        ui = self.ui.point_sources
        #The top of the task pane is where users define/select PS regions
        #Icons to add/remove/duplicate boundary conditions are given at the top
        #Clicking the 'add' and 'duplicate' buttons triggers a popup window where the user must select
        #a region to apply the point source.
        ui.toolbutton_add.clicked.connect(self.pss_show_regions_popup)
        ui.toolbutton_delete.clicked.connect(self.pss_delete_regions)
        ui.toolbutton_delete.setEnabled(False) # Need a selection

        ui.tablewidget_regions.itemSelectionChanged.connect(self.handle_pss_region_selection)

        self.pss_current_tab = FLUID_TAB # If fluid is disabled, we will switch
        self.pss_current_solid = self.P = None
        ui.pushbutton_fluid.pressed.connect(lambda: self.pss_change_tab(FLUID_TAB,None))

        # Trim width of "Fluid" button, like we do for
        # dynamically-created "Solid #" buttons
        b =  ui.pushbutton_fluid
        w = b.fontMetrics().boundingRect(b.text()).width() + 20
        b.setMaximumWidth(w)


    def pss_show_regions_popup(self):
        # Users cannot select inapplicable regions.
        # PS regions can be points, planes, or volumes (not STLs)
        # No region can define more than one point source.
        pass

    def pss_cancel_add(self):
        ui = self.ui.point_sources

        for item in (ui.toolbutton_add,
                     ui.tablewidget_regions):
            item.setEnabled(True)

        if get_selected_row(ui.tablewidget_regions) is not None:
            for item in (ui.detail_pane,
                         ui.toolbutton_delete):
                item.setEnabled(True)


    def pss_add_regions(self):
        # Interactively add regions to define PSs
        ui = self.ui.point_sources
        rp = self.regions_popup
        self.pss_cancel_add() # Reenable input widgets
        selections = rp.get_selection_list()
        if not selections:
            return
        self.pss_add_regions_1(selections) # Indices will be assigned
        self.pss_setup_current_tab() # Update the widgets

    def pss_add_regions_1(self, selections, indices=None):
        # Used by both interactive and load-time add-region handlers

        if self.pss_region_dict is None:
            self.pss_region_dict = self.ui.regions.get_region_dict()

        ui = self.ui.initial_conditions
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
                idx = self.pss_find_index()
                indices[i] = idx
            self.pss[idx] = {'region': region_name}
            region_data = self.pss_region_dict.get(region_name)
            if region_data is None: # ?
                self.warn("no data for region %s" % region_name)
                continue
            self.pss_set_region_keys(region_name, idx, region_data)
            self.pss_region_dict[region_name]['available'] = False # Mark as in-use

        item.setData(UserRole, (tuple(indices), tuple(selections)))
        tw.setItem(nrows, 0, item)

        self.fixup_pss_table(tw)

        if autoselect:
            tw.setCurrentCell(nrows, 0)


    def pss_find_index(self):
        n = 1
        while n in self.pss:
            n += 1
        return n


    def pss_delete_regions(self):
        ui = self.ui.initial_conditions
        tw = ui.tablewidget_regions
        row = get_selected_row(tw)
        if row is None: # No selection
            return

        # Unset keywords
        kwlist = list(self.project.keywordItems())
        for kw in kwlist:
            key, args = kw.key, kw.args
            # TODO use keyword_args here instead of startswith
            if key.startswith('ic_') and args and args[0] in self.pss_current_indices:
                self.unset_keyword(key, args=args)

        # TODO: fix any resulting holes in index sequence!

        for r in self.pss_current_regions:
            if r in self.pss_region_dict:
                self.pss_region_dict[r]['available'] = True

        for i in self.pss_current_indices:
            del self.pss[i]

        self.pss_current_regions = []
        self.pss_current_indices = []

        tw.removeRow(row)
        self.pss_setup_current_tab()


    def handle_pss_region_selection(self):
        ui = self.ui.initial_conditions
        table = ui.tablewidget_regions
        row = get_selected_row(table)
        if row is None:
            indices = []
            regions = []
        else:
            (indices, regions) = table.item(row,0).data(UserRole)
        self.pss_current_indices, self.pss_current_regions = indices, regions
        enabled = (row is not None)
        ui.toolbutton_delete.setEnabled(enabled)
        ui.detail_pane.setEnabled(enabled)
        if not enabled:
            # Clear
            for widget in widget_iter(ui.detail_pane):
                if isinstance(widget, LineEdit):
                    widget.setText('')
            return
        self.pss_setup_current_tab() # reinitialize all widgets in current tab


    def fixup_pss_table(self, tw, stretch_column=0):
        ui = self.ui.initial_conditions
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







    def setup_pss(self):
        #Tabs group point source parameters for phases. Tabs are unavailable if no input
        #is required from the user.
        #    Fluid tab - Unavailable if the fluid phase was disabled.
        b = ui.pushbutton_fluid
        b.setText(self.fluid_phase_name)
        b.setEnabled(not self.fluid_solver_disabled)
        if self.fluid_solver_disabled:
            if self.bcs_current_tab == 0: # Don't stay on disabled tab
                self.bcs_change_tab(SOLIDS_TAB, 1) # what if no solids?  we shouldn't be here
        font = b.font()
        font.setBold(self.bcs_current_tab == 0)
        b.setFont(font)

        #    Each solid phase will have its own tab. The tab name should be the name of the solid



    def pss_setup_fluid_tab(self):
        #Fluid (tab)

        #    Define mass flowrate:
        # Sets keyword PS_MASSFLOW_G(#)
        # DEFAULT value of 0.0

        #    Define temperature
        # Specification only available when solving energy equations
        # DEFAULT value 293.15
        # Sets keyword PS_T_G(#)

        #    Select species and set mass fractions (table format)
        # Specification only available when solving species equations
        # DEFAULT value 1.0 of last defined species
        # Sets keyword PS_X_G(#,#)
        # Error check: if specified, mass fractions must sum to one

        #    Define X-axial velocity:
        # Specification always
        # No DEFAULT value (blank)
        # Sets keyword PS_U_G(#)
        #    Define Y-axial velocity:
        # Specification always
        # No DEFAULT value (blank)
        # Sets keyword PS_V_G(#)
        #    Define Z-axial velocity:
        # Specification always
        # No DEFAULT value (blank)
        # Sets keyword PS_W_G(#)
        pass


    def setup_pss_solids_tab(self, P):
        #Solids-# (tab) - At this time, only TFM solids can be defined with point sources.
        #At some point in the future, this could be extended to PIC solids, but never DEM.

        #Define mass flowrate:
        # Select mass inflow specification type:
        # Sets keyword PS_MASSFLOW_S(#,#)
        # DEFAULT value of 0.0

        #Define temperature
        # Specification only available when solving energy equations
        # DEFAULT value 293.15
        # Sets keyword PS_T_S(#,#)

        #Select species and set mass fractions (table format)
        # Specification only available when solving species equations
        # DEFAULT value 1.0 of last defined species
        # Sets keyword PS_X_S(#,#,#)
        # Error check: if specified, mass fractions must sum to one

        #Define X-axial velocity:
        # Specification always
        # No DEFAULT value (blank)
        # Sets keyword PS_U_S(#,#)

        #Define Y-axial velocity:
        # Specification always
        # No DEFAULT value (blank)
        # Sets keyword PS_V_S(#,#)

        #Define Z-axial velocity:
        # Specification always
        # No DEFAULT value (blank)
        # Sets keyword PS_W_S(#,#)
        pass


    def setup_pss_scalar_tab():
        # No scalar tab for point sources
        pass
