# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict
import numpy as np

from qtpy import QtCore, QtWidgets, PYQT5


MESH_EXTENT_KEYS = ['xmin', 'xlength', 'ymin', 'ylength', 'zmin', 'zlength']
MESH_CELL_KEYS = ['imax', 'jmax', 'kmax']

def safe_float(value):
    """try to convert the value to a float, if ValueError, send error to gui
    and rturn 0"""
    try:
        return float(value)
    except ValueError:
        return 0.0


def safe_int(value):
    """try to convert the value to a int, if ValueError, send error to gui
    and rturn 0"""
    try:
        return int(value)
    except ValueError:
        return 0

class KeyHandler(QtCore.QObject):
    value_updated = QtCore.Signal(object, object, object)
    update_mesh_extents = QtCore.Signal(object, object, object)
    update_mesh_cells = QtCore.Signal(object, object, object)
    def updateValue(self, key, val, args):
        if key in MESH_EXTENT_KEYS:
            self.update_mesh_extents.emit(key, val, args)
        elif key in MESH_CELL_KEYS:
            self.update_mesh_cells.emit(key, val, args)


class Mesh(object):
    # Mesh Task Pane Window: This section allows a user to define the mesh
    # Methods that need vtk are in the vtkwidget
    def init_mesh(self):
        self.mesh_extents = []
        self.mesh_cells = []
        ui = self.ui.mesh
        self.cell_spacing_widgets = [ui.lineedit_mesh_cells_size_x, ui.lineedit_mesh_cells_size_y,
            ui.lineedit_mesh_cells_size_z]

        # key hadler
        self.mesh_key_handler = KeyHandler()
        self.project.register_widget(self.mesh_key_handler, MESH_EXTENT_KEYS + MESH_CELL_KEYS, [])
        self.mesh_key_handler.update_mesh_extents.connect(self.update_background_mesh_extents)
        self.mesh_key_handler.update_mesh_cells.connect(self.update_background_mesh_cells)

        # connect mesh tab btns
        for i, btn in enumerate([ui.pushbutton_mesh_uniform, ui.pushbutton_mesh_mesher]):
            btn.clicked.connect(lambda ignore, i=i, btn=btn:self.change_mesh_tab(i, btn))

        ui.combobox_mesher.currentIndexChanged.connect(
            self.change_mesher_options)
        ui.checkbox_internal_external_flow.value_updated.connect(self.toggle_out_stl_value)

        # setup tables
        for table in [ui.table_mesh_control_points_x, ui.table_mesh_control_points_y, ui.table_mesh_control_points_z]:
            table.dtype = OrderedDict
            table._setModel() # FIXME: Should be in __init__
            table.set_selection_model()
            table.set_value(OrderedDict())
            table.set_columns(['position', 'cells', 'stretch'])
            table.show_vertical_header(True)
            table.auto_update_rows(True)
#            table.new_selection.connect(self.update_region_parameters)
#            table.clicked.connect(self.cell_clicked)
            table.default_value = OrderedDict()
#            table.value_changed.connect(self.table_value_changed)

    def toggle_out_stl_value(self, wid, value, args):
        val = -1.0
        if list(value.values())[0]:
            val = 1.0
        self.update_keyword('out_stl_value', val)

    def change_mesh_tab(self, tabnum, btn):
        """switch mesh stacked widget based on selected"""
        ui = self.ui.mesh
        self.animate_stacked_widget(
            ui.stackedwidget_mesh,
            ui.stackedwidget_mesh.currentIndex(),
            tabnum,
            direction='horizontal',
            line=ui.line_mesh,
            to_btn=btn,
            btn_layout=ui.gridlayout_mesh_tab_btns,
            )

    def change_mesher_options(self):
        """switch the mesh options stacked widget"""
        ui = self.ui.mesh
        mesher = str(ui.combobox_mesher.currentText()).lower()

        current_index = 0
        for i in range(ui.stackedwidget_mesher_options.count()):
            widget = ui.stackedwidget_mesher_options.widget(i)
            if mesher == str(widget.objectName()).lower():
                current_index = i
                break

        self.animate_stacked_widget(
            ui.stackedwidget_mesher_options,
            ui.stackedwidget_mesher_options.currentIndex(),
            current_index,
            direction='horizontal',
            )

        enable = mesher != 'cutcell'
        ui.pushbutton_generate_mesh.setEnabled(enable)
        ui.pushbutton_remove_mesh.setEnabled(enable)

    def init_background_mesh(self):
        self.mesh_extents = []
        for key in MESH_EXTENT_KEYS:
            self.mesh_extents.append(safe_float(self.project.get_value(key, default=0.0)))

        self.mesh_cells = []
        for key in MESH_CELL_KEYS:
            self.mesh_cells.append(safe_int(self.project.get_value(key, default=0)) + 1)

        self.update_background_mesh()

    def update_background_mesh_cells(self, key, val, args):
        """collect cells changes, check if value is different"""
        if not self.mesh_cells: return
        val = safe_int(val)
        ind = MESH_CELL_KEYS.index(key)
        old_val = self.mesh_cells[ind]
        if old_val != val:
            self.mesh_cells[ind] = val
            self.update_background_mesh()

    def update_background_mesh_extents(self, key, val, args):
        """collect extents changes, check if value is different"""
        if not self.mesh_extents: return
        val = safe_float(val)
        ind = MESH_EXTENT_KEYS.index(key)
        old_val = self.mesh_extents[ind]
        if old_val != val:
            self.mesh_extents[ind] = val
            self.update_background_mesh()

    def update_background_mesh(self):
        """update the background mesh"""
        extents = self.mesh_extents
        cells = self.mesh_cells

        # average cell width
        for (f, t), c, wid in zip(zip(extents[::2], extents[1::2]), cells, self.cell_spacing_widgets):
            if c-1 > 0:
                w = (t-f)/(c-1)
            else:
                w = t-f
            wid.setText('{0:.2e}'.format(w))

        # determine cell spacing
        spacing = [
            np.linspace(extents[0], extents[1], cells[0]),
            np.linspace(extents[2], extents[3], cells[1]),
            np.linspace(extents[4], extents[5], cells[2])
            ]

        self.vtkwidget.update_background_mesh(spacing)


