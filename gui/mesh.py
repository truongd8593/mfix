# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict
import numpy as np

from qtpy import QtCore, QtWidgets, PYQT5

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
    update_mesh = QtCore.Signal(bool)
    def updateValue(self, key, val, args):
        self.update_mesh.emit(True)

class Mesh(object):
    # Mesh Task Pane Window: This section allows a user to define the mesh
    # Methods that need vtk are in the vtkwidget
    def init_mesh(self):
        ui = self.ui.mesh

        self.cell_spacing_widgets = [
            ui.lineedit_mesh_cells_size_x,
            ui.lineedit_mesh_cells_size_y,
            ui.lineedit_mesh_cells_size_z
            ]

        self._meshhandler = KeyHandler()
        self.project.register_widget(
            self._meshhandler, ['xmin', 'xlength', 'ymin', 'ylength', 'zmin',
            'zlength', 'imax', 'jmax', 'kmax'], [])
        self._meshhandler.update_mesh.connect(self.update_background_mesh)

        # connect mesh tab btns
        for i, btn in enumerate([ui.pushbutton_mesh_uniform,
                                 ui.pushbutton_mesh_mesher]):
            btn.clicked.connect(lambda ignore, i=i, btn=btn:self.change_mesh_tab(i, btn))

        ui.combobox_mesher.currentIndexChanged.connect(
            self.change_mesher_options)

        ui.checkbox_internal_external_flow.value_updated.connect(self.toggle_out_stl_value)

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

    def update_background_mesh(self):
        """update the background mesh"""
        extents = []
        for key in ['xmin', 'xlength', 'ymin', 'ylength', 'zmin', 'zlength']:
            extents.append(safe_float(self.project.get_value(key, default=0.0)))

        cells = []
        for key in ['imax', 'jmax', 'kmax']:
            cells.append(safe_int(self.project.get_value(key, default=0)) + 1)

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
