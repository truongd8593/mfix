# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict
from itertools import groupby

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


def ctrl_pts_to_loc(ctrl):
    """given control points, convert to location"""
    loc = [0]
    last = 0
    for pt in ctrl.values():
        er = pt['stretch']
        cp = pt['position']
        nc = pt['cells']

        dx = (cp-last)/nc
        for i in range(nc):
            # no expansion
            if er == 1:
                loc.append(loc[-1] + dx)
        last = loc[-1]
    return loc


def linspace(f, t, c):
    """copy of numpy's linspace,"""
    if c == 1:
        return [f, t]
    dx = (t-f)/float(c-1)
    l = [0]
    for i in range(c-1):
        l.append(l[-1]+dx)
    return l


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
        self.mesh_spacing = [[],[],[]]
        ui = self.ui.mesh
        self.cell_spacing_widgets = [ui.lineedit_mesh_cells_size_x, ui.lineedit_mesh_cells_size_y,
            ui.lineedit_mesh_cells_size_z]
        self.mesh_spacing_tables = [ui.table_mesh_control_points_x, ui.table_mesh_control_points_y, ui.table_mesh_control_points_z]

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
        for table in self.mesh_spacing_tables:
            table.dtype = OrderedDict
            table._setModel() # FIXME: Should be in __init__
            table.set_selection_model()
            table.set_value(OrderedDict())
            table.set_columns(['position', 'cells', 'stretch', 'first', 'last'])
            table.show_vertical_header(True)
            table.auto_update_rows(True)
#            table.new_selection.connect(self.update_region_parameters)
#            table.clicked.connect(self.cell_clicked)
            table.default_value = OrderedDict()
#            table.value_changed.connect(self.table_value_changed)
#            table.fit_to_contents()

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
        prj = self.project
        self.mesh_extents = []
        self.mesh_spacing = [[],[],[]]
        for key in MESH_EXTENT_KEYS:
            self.mesh_extents.append(safe_float(prj.get_value(key, default=0.0)))

        self.mesh_cells = []
        for key in MESH_CELL_KEYS:
            self.mesh_cells.append(safe_int(prj.get_value(key, default=0)) + 1)

        # collect dx, dy, dx
        for i, (s, c, e) in enumerate(zip(['dx', 'dy', 'dz'], MESH_CELL_KEYS, MESH_EXTENT_KEYS[1::2])):
            d = [prj.get_value(s, args=args) for args in sorted(prj.get_key_indices(s))]
            l = len(d)

            # if there are spacing, update the keywords.
            if l>0:
                self.mesh_cells[i] = l
                self.update_keyword(c, l)
                m = sum(d)
                self.mesh_extents[i*2+1] = m
                self.update_keyword(e, m)
                self.mesh_spacing[i] = d
                self.extract_mesh_spacing(i, d)

        # collect variable grid spacing keywords
        for i, k in enumerate(['x', 'y', 'z']):
            indices = prj.get_key_indices('cp' + k)
            if indices:
                table_dic = OrderedDict()
                for j, ind in enumerate(sorted(indices)):
                    table_dic[j] = {
                        'position': prj.get_value('cp' + k, 0, args=ind),
                        'cells': prj.get_value('nc' + k, 1, args=ind),
                        'stretch': prj.get_value('er' + k, 1, args=ind),
                        'first': prj.get_value('first_d' + k, 0, args=ind),
                        'last': prj.get_value('last_d' + k, 0, args=ind)}
                self.mesh_spacing_tables[i].set_value(table_dic)
                self.mesh_spacing_tables[i].fit_to_contents()

        self.update_background_mesh()

    def extract_mesh_spacing(self, index, spacing):
        """given a list of cell spacing, convert to control points"""

        start = 0
        table_dic = OrderedDict()
        d = ['x','y','z'][index]

        for i, (val, count)  in enumerate([(k, sum(1 for i in g)) for k, g in groupby(spacing)]):
            loc = count*val + start
            start = loc
            self.update_keyword('cp' + d, loc, args=i)
            self.update_keyword('nc' + d, count, args=i)
            table_dic[i] = {'position': loc, 'cells': count, 'stretch': 1,
                            'first': 0, 'last': 0}
        for i in range(len(spacing)):
            self.update_keyword('d' + d, None, args=i)

        self.mesh_spacing_tables[index].set_value(table_dic)
        self.mesh_spacing_tables[index].fit_to_contents()

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
        spacing = []
        for i, table in enumerate(self.mesh_spacing_tables):
            val = table.value
            if val:
                spacing.append(ctrl_pts_to_loc(val))
            else:
                spacing.append(linspace(extents[i*2], extents[i*2+1], cells[i]))

        self.vtkwidget.update_background_mesh(spacing)
