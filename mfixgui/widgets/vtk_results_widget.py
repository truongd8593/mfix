# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict, Mapping
import glob
import os
import copy
from bisect import  bisect_left
from xml.etree import ElementTree

# graphics libraries
try:
    import pyqtgraph as pg
    PYQTGRAPH_AVAILABLE = True
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    DEFAULT_PEN = pg.mkPen(color='#64B5F6', width=2)
except ImportError:
    pg = None
    PYQTGRAPH_AVAILABLE = False
    DEFAULT_PEN = 'b'
except RuntimeError:
    pg = None
    PYQTGRAPH_AVAILABLE = False
    DEFAULT_PEN = 'b'

try:
    import vtk
    VTK_AVAILABLE = True
    from mfixgui.colormaps.color_maps import build_vtk_lookup_tables, build_qicons
    LOOKUP_TABLES = build_vtk_lookup_tables()
    GLYPHS = {
        'sphere':   vtk.vtkSphereSource,
        'cube':      vtk.vtkCubeSource,
        'cylinder': vtk.vtkCylinderSource,
        'cone':     vtk.vtkConeSource}
except:
    vtk = None
    VTK_AVAILABLE = False
    build_vtk_lookup_tables, get_color_map_pngs = None, None
    LOOKUP_TABLES = {}

from qtpy import QtCore, QtWidgets, QtGui, uic
from mfixgui.widgets.base_vtk import BaseVtkWidget, vtk, VTK_AVAILABLE
from mfixgui.widgets.base import CustomPopUp
from mfixgui.tools.general import get_icon

PLOT_ITEMS = OrderedDict([
    ['Select an item', {}],
    ['dt', {'left':'dt', 'bottom':'Time Step', 'var':'dt'}],
    ['nit', {'left':'Number of Iterations', 'bottom':'Time Step', 'var':'nit'}],
    ['time', {'left':'Simulation Time [s]', 'bottom':'Ellapsed Wall Time [s]', 'var':'time', 'var2':'walltime_elapsed'}],
    ])

SETTINGS = QtCore.QSettings('MFIX', 'MFIX')
DEFAULT_PLAYBACK_SPEED = 100
DEFAULT_MAXIMUM_POINTS = 10000
DEFAULT_GEO_COLOR = QtGui.QColor(224, 224, 224)

UI_FILE_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), os.pardir, 'uifiles')


def parse_pvd_file(fname):
    '''given a pvd file, return a dict of time (float) : file_name'''
    f_dict = OrderedDict()
    if not os.path.exists(fname): return f_dict
    try:
        tree = ElementTree.parse(fname)
    except:
        return f_dict
    root = tree.getroot()
    for data_set in root.iter('DataSet'):
        f_dict[float(data_set.attrib['timestep'])] = data_set.attrib['file']

    return f_dict

def build_time_dict(search_str):
    '''given a glob search string, return a dictionary of
    time (float) : file_name'''
    f_dict = OrderedDict()

    files = glob.glob(search_str)
    for f in sorted(files):
        time = None
        with open(f) as xml:
            for i, line in enumerate(xml):
                if '<!-- Time =' in line:
                    try:
                        time = float(line.replace('<!-- Time =', '').replace('sec. -->', ''))
                    except:
                        pass
                    break

                if i > 4:
                    break
        if time is not None:
            f_dict[time] = os.path.basename(f)
    return f_dict

def update(d, u):
    for k, v in u.items():
        if isinstance(v, Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d

class ColorMapPopUp(QtWidgets.QDialog):
    applyEvent = QtCore.Signal(object, object, object)
    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self, parent)

        self.color = None
        ui = self.ui = uic.loadUi(os.path.join(UI_FILE_DIR, 'color_map.ui'), self)

        self.setWindowTitle('Color Map')

        ui.toolbutton_select_color.clicked.connect(self.select_color)
        ui.toolbutton_select_color.setStyleSheet("QToolButton{{ background: {};}}".format(DEFAULT_GEO_COLOR.name()))
        ui.lineedit_from.dtype = float
        ui.lineedit_to.dtype = float

        for name, icons in build_qicons().items():
            if not name.endswith('_reversed'):
                ui.combobox_color_map.addItem(icons.get('bar', QtGui.QIcon()), name)
        self.set_color_map('viridis')

        btn = ui.buttonBox.button(QtWidgets.QDialogButtonBox.Apply)
        btn.clicked.connect(self.emit_apply_event)

        ui.toolbutton_auto_scale.clicked.connect(self.auto_scale)

    def emit_apply_event(self):
        self.applyEvent.emit(self.geo, self.button, self.array)

    def set_(self, array):
        self.array = array
        self.color = self.array.get('color', DEFAULT_GEO_COLOR)
        self.set_color(self.color)
        self.set_color_map(self.array.get('color_map', 'viridis'))
        self.ui.checkbox_reversed.setChecked(self.array.get('reversed', False))

        d_range = [0, 1]
        self.ui.lineedit_from.updateValue(None,
            self.array.get('from', self.array.get('range', d_range)[0]))
        self.ui.lineedit_to.updateValue(None,
            self.array.get('to', self.array.get('range', d_range)[1]))

        single_color = self.array.get('single_color', False)
        self.ui.checkbox_single_color.setChecked(single_color)
        self.ui.widget_color_map.setEnabled(not single_color)

    def get(self):
        color_map = self.ui.combobox_color_map.currentText()
        reverse = self.ui.checkbox_reversed.isChecked()
        if reverse:
            color_map += '_reversed'
        return {
            'color':        self.color,
            'single_color': self.ui.checkbox_single_color.isChecked(),
            'color_map':    color_map,
            'reversed':     reverse,
            'from':         self.ui.lineedit_from.value,
            'to':           self.ui.lineedit_to.value,
            }

    def select_color(self):
        col = QtWidgets.QColorDialog.getColor()
        if not col.isValid():
            return
        self.color = col
        self.set_color(col)

    def set_color(self, color):
        if isinstance(color, QtGui.QColor):
            self.ui.toolbutton_select_color.setStyleSheet("QToolButton{{ background: {};}}".format(
                color.name()))

    def set_color_map(self, color_map):
        color_map = color_map.replace('_reversed', '')
        index = self.ui.combobox_color_map.findText(color_map)
        self.ui.combobox_color_map.setCurrentIndex(index)

    def auto_scale(self):
        d_range = [0, 1]
        self.ui.lineedit_from.updateValue(None,
            self.array.get('range', d_range)[0])
        self.ui.lineedit_to.updateValue(None,
            self.array.get('range', d_range)[1])


class ParticleOptions(QtWidgets.QDialog):
    applyEvent = QtCore.Signal()
    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self, parent)

        self.color = None
        ui = self.ui = uic.loadUi(os.path.join(UI_FILE_DIR, 'particle_options.ui'), self)

        ui.lineedit_maximum_particles.updateValue(None, DEFAULT_MAXIMUM_POINTS)
        ui.lineedit_maximum_particles.dtype = int
        self.setWindowTitle('Particle Options')

        btn = self.ui.buttonBox.button(QtWidgets.QDialogButtonBox.Apply)
        btn.clicked.connect(self.emit_apply_event)

    def emit_apply_event(self):
        self.applyEvent.emit()

    def get(self):
        """return options"""
        return {
            'max_points': self.ui.lineedit_maximum_particles.value,
            'mapper': self.ui.combobox_mapper.currentText(),
            'glyph': self.ui.combobox_glyph.currentText(),
            }


class GraphicsVtkWidget(BaseVtkWidget):
    """vtk widget for showing results"""
    def __init__(self, parent=None):
        BaseVtkWidget.__init__(self, parent)

        # find MfixGui
        def find_gui(parent):
            if parent.objectName() == 'mfixgui':
                return parent
            else:
                return find_gui(parent.parent())

        self.mfixgui = find_gui(self)

        self.cell_arrays = {}
        self.node_arrays = {}
        self.point_arrays = {}
        self.frame_index = -1

        self.play_timer = QtCore.QTimer()
        self.play_timer.timeout.connect(self.forward)

        # look for vtu files
        self.project_file = SETTINGS.value('project_file')
        self.project_name = os.path.splitext(os.path.basename(self.project_file))[0]
        self.project_dir = os.path.dirname(self.project_file)

        self.file_watcher = QtCore.QFileSystemWatcher()
        self.file_watcher.addPath(self.project_dir)
        self.file_watcher.directoryChanged.connect(self.look_for_files)

        # look for files
        self.look_for_files()

        # dialogs
        self.color_dialog  = ColorMapPopUp(self)
        self.color_dialog.applyEvent.connect(self.change_color)

        self.particle_option_dialog = ParticleOptions(self)
        self.particle_option_dialog.applyEvent.connect(self.change_particle_options)

        self.init_toolbar()
        self.init_vtk()
        self.init_geometry()

        self.change_frame(0)
        self.reset_view()

    def init_vtk(self):

        self.actors = {}
        self.mappers = {}
        self.lookuptables = {}

        self.ugrid_cell_mapper = None
        self.ugrid_node_mapper = None
        self.particle_mapper = None

    def init_ugrid(self):
        '''setup the cell/point vtk stuff'''
        # cells
        self.enable_toolbar_geo('cells', visible=False)
        self.ugrid_reader = vtk.vtkXMLUnstructuredGridReader()

        self.ugrid_cell_mapper = self.mappers['cells'] = vtk.vtkDataSetMapper()
        self.ugrid_cell_mapper.SetInputConnection(self.ugrid_reader.GetOutputPort())
        self.ugrid_cell_mapper.SetScalarVisibility(True)
        self.ugrid_cell_mapper.SetScalarModeToUseCellFieldData()
        self.change_color_bar('cells', 'viridis')

        actor = self.actors['cells'] = vtk.vtkActor()
        actor.SetMapper(self.ugrid_cell_mapper)

        self.vtkrenderer.AddActor(actor)
        self.change_visibility('cells', False) # hide cells by default

        # points
        self.enable_toolbar_geo('nodes')

        self.ugrid_cell_to_points = vtk.vtkCellDataToPointData()
        self.ugrid_cell_to_points.SetInputConnection(self.ugrid_reader.GetOutputPort())

        self.ugrid_node_mapper = self.mappers['nodes'] = vtk.vtkDataSetMapper()
        self.ugrid_node_mapper.SetInputConnection(self.ugrid_cell_to_points.GetOutputPort())
        self.ugrid_node_mapper.SetScalarVisibility(True)
        self.ugrid_node_mapper.SetScalarModeToUsePointFieldData()
        self.change_color_bar('nodes', 'viridis')

        actor = self.actors['nodes'] = vtk.vtkActor()
        actor.SetMapper(self.ugrid_node_mapper)

        self.vtkrenderer.AddActor(actor)

    def init_particles(self):
        '''setup the particle vtk stuff'''
        self.enable_toolbar_geo('points')
        self.particle_reader = vtk.vtkXMLPolyDataReader()

        ss = vtk.vtkSphereSource()

        self.glyph_mask = vtk.vtkMaskPoints()
        self.glyph_mask.SetInputConnection(self.particle_reader.GetOutputPort())
        self.glyph_mask.RandomModeOn()
        self.glyph_mask.SetRandomModeType(1)
        self.glyph_mask.SetMaximumNumberOfPoints(DEFAULT_MAXIMUM_POINTS)

        self.glyph = vtk.vtkGlyph3D()
        self.glyph.SetInputConnection(self.glyph_mask.GetOutputPort())
        self.glyph.SetSourceConnection(ss.GetOutputPort())
        self.glyph.SetColorModeToColorByVector()
#        self.glyph.SetColorModeToColorByScalar()

        self.particle_mapper = self.mappers['points'] = vtk.vtkPolyDataMapper()
        self.particle_mapper.SetInputConnection(self.glyph.GetOutputPort())
#        self.particle_mapper.SetScalarModeToUsePointFieldData()
        self.change_color_bar('points', 'viridis')

        actor = self.actors['points'] = vtk.vtkActor()
        actor.SetMapper(self.particle_mapper)

        self.vtkrenderer.AddActor(actor)

    def init_geometry(self):
        self.enable_toolbar_geo('geometry')
        poly_data = self.mfixgui.vtkwidget.collect_toplevel_geometry()

        # Create a mapper
        mapper = self.mappers['geometry'] = vtk.vtkPolyDataMapper()
        if poly_data is not None:
            mapper.SetInputConnection(poly_data.GetOutputPort())

        # Create an actor
        actor = self.actors['geometry'] = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(DEFAULT_GEO_COLOR.getRgbF()[:3])
        actor.GetProperty().SetOpacity(0.4)

        self.vtkrenderer.AddActor(actor)

    def enable_toolbar_geo(self, geo, visible=True):
        for name, wid in self.visual_btns[geo].items():
            if name == 'component':
                continue
            if visible and name == 'visible':
                wid.setChecked(True)
            wid.setEnabled(True)

    def init_toolbar(self):
        self.init_base_toolbar()

        # more buttons
        self.toolbutton_visible = QtWidgets.QToolButton()
        ##self.toolbutton_visible.setCheckable(True)
        self.toolbutton_visible.setIcon(get_icon('visibility.png'))
        self.visible_menu = CustomPopUp(self, self.toolbutton_visible)
        self.visible_menu.finished.connect(lambda ignore: self.toolbutton_visible.setDown(False))
        self.toolbutton_visible.pressed.connect(self.toggle_visible_menu)

        # --- visual representation menu ---
        layout = self.visible_menu.layout
        layout.setContentsMargins(0, 5, 5, 5)
        self.visual_btns = {}
        for i, geo in enumerate(['Cells', 'Nodes', 'Points', 'Geometry']):
            geo_name = geo
            geo = geo.lower().replace(' ', '_')
            btns = self.visual_btns[geo] = {}
            # tool button
            toolbutton = QtWidgets.QToolButton(self.visible_menu)
            toolbutton.clicked.connect(lambda down, g=geo: self.change_visibility(g, down))
            toolbutton.setCheckable(True)
            toolbutton.setChecked(False)
            toolbutton.setAutoRaise(True)
            toolbutton.setIcon(get_icon('visibility.png'))
            toolbutton.setEnabled(False)
            layout.addWidget(toolbutton, i, 0)
            btns['visible'] = toolbutton

            combo = None
            if not geo == 'geometry':
                # color by
                combo = QtWidgets.QComboBox(self.visible_menu)
                combo.activated.connect(lambda item, g=geo, c=combo: self.change_color_by(g, c))
                combo.setEnabled(False)
                layout.addWidget(combo, i, 1)
                btns['color_by'] = combo

                # component
                combo2 = QtWidgets.QComboBox(self.visible_menu)
                combo2.activated.connect(lambda item, g=geo, c=combo, c2=combo2: self.change_color_by(g, c, c2))
                combo2.addItems(['mag', 'x', 'y', 'z'])
                combo2.setEnabled(False)
                layout.addWidget(combo2, i, 2)
                btns['component'] = combo2


            toolbutton = QtWidgets.QToolButton(self.visible_menu)
            size = QtCore.QSize(25, 25)
            toolbutton.setMinimumSize(size)
            toolbutton.setMaximumSize(size)
            toolbutton.setIconSize(size)
            if not geo == 'geometry':
                toolbutton.clicked.connect(lambda ignore, g=geo, t=toolbutton, c=combo: self.handle_change_color(g, t, c))
                toolbutton.setIcon(build_qicons().get('viridis', {}).get('icon', QtGui.QIcon))
            else:
                toolbutton.clicked.connect(lambda ignore, t=toolbutton: self.change_geo_color(t))
                toolbutton.setStyleSheet("QToolButton{{ background: {};}}".format(DEFAULT_GEO_COLOR.name()))
            toolbutton.setAutoRaise(True)
            toolbutton.setEnabled(False)
            layout.addWidget(toolbutton, i, 3)
            btns['color'] = toolbutton

            opacity = QtWidgets.QDoubleSpinBox(self.visible_menu)
            opacity.setRange(0, 1)
            if geo == 'geometry':
                opacity.setValue(0.4)
            else:
                opacity.setValue(1.0)
            opacity.setSingleStep(0.1)
            opacity.valueChanged.connect(lambda o, g=geo: self.change_opacity(g, o))
            opacity.setEnabled(False)
            layout.addWidget(opacity, i, 4)
            btns['opacity'] = opacity

            # label
            label = QtWidgets.QLabel(geo_name, self.visible_menu)
            if geo == 'points':
                layout.addWidget(label, i, 10)
            else:
                layout.addWidget(label, i, 10, 1, 2)

            # more
            if geo in ['points']:
                toolbutton = QtWidgets.QToolButton()
                toolbutton.setAutoRaise(True)
                toolbutton.setIcon(get_icon('right.png'))
                toolbutton.setEnabled(False)
                toolbutton.clicked.connect(self.handle_particle_options)
                layout.addWidget(toolbutton, i, 11)
                btns['more'] = toolbutton

        self.toolbutton_back = QtWidgets.QToolButton(self.visible_menu)
        self.toolbutton_back.clicked.connect(self.handle_begining)
        self.toolbutton_back.setIcon(get_icon('previous.png'))

        self.toolbutton_play = QtWidgets.QToolButton(self.visible_menu)
        self.toolbutton_play.clicked.connect(self.handle_play_stop)
        self.toolbutton_play.setIcon(get_icon('play.png'))

        self.toolbutton_forward = QtWidgets.QToolButton(self.visible_menu)
        self.toolbutton_forward.clicked.connect(self.handle_end)
        self.toolbutton_forward.setIcon(get_icon('next.png'))

        self.frame_spinbox = QtWidgets.QSpinBox(self.visible_menu)
        self.frame_spinbox.valueChanged.connect(self.change_frame)
        self.frame_spinbox.setMaximum(9999999)

        self.toolbutton_play_speed = QtWidgets.QToolButton(self.visible_menu)
        self.toolbutton_play_speed.setIcon(get_icon('speed.png'))

        self.speed_menu = CustomPopUp(self, self.toolbutton_play_speed)
        self.speed_menu.finished.connect(lambda ignore: self.toolbutton_play_speed.setDown(False))
        self.speed_slider = QtWidgets.QSlider(self.speed_menu)
        self.speed_slider.setRange(0, 1000)
        self.speed_slider.setValue(DEFAULT_PLAYBACK_SPEED)
        self.speed_slider.setOrientation(QtCore.Qt.Horizontal)
        self.speed_slider.valueChanged.connect(self.handle_speed_changed)
        self.speed_menu.layout.addWidget(self.speed_slider)
        self.toolbutton_play_speed.pressed.connect(self.speed_menu.popup)

        self.checkbox_snap = QtWidgets.QCheckBox('Save Snapshots', self.visible_menu)

        for btn in [self.toolbutton_visible, self.toolbutton_back,
                    self.toolbutton_play, self.toolbutton_forward,
                    self.frame_spinbox, self.toolbutton_play_speed,
                    self.checkbox_snap]:
            self.button_bar_layout.addWidget(btn)
            if isinstance(btn, QtWidgets.QToolButton):
                btn.setAutoRaise(True)

        self.button_bar_layout.addStretch()

    def showEvent(self, event):
        # has to be called after the widget is visible
        self.vtkiren.Initialize()

    def hideEvent(self, event):
        self.stop()

    def close(self):
        BaseVtkWidget.close(self)

        # clean up timer
        self.play_timer.stop()


    def toggle_visible_menu(self):
        if self.visible_menu.isVisible():
            self.visible_menu.hide()
        else:
            self.visible_menu.popup()



    def show_visible_menu(self):
        # update comboboxes based on avaliable arrays
        for type_, array in [('cells', self.cell_arrays),
                             ('nodes', self.node_arrays),
                             ('points', self.point_arrays)]:
            btns = self.visual_btns[type_]
            combo = btns['color_by']
            text = combo.currentText()
            combo.clear()
            if array:
                combo.addItems(array.keys())
                index = combo.findText(text)
                combo.setCurrentIndex(index)
            combo.setEnabled(bool(array))

        self.visible_menu.popup()

    def handle_speed_changed(self):
        if self.play_timer.isActive():
            self.play_timer.stop()
            self.handle_play_stop()

    def stop(self):
        self.toolbutton_play.setIcon(get_icon('play.png'))
        self.play_timer.stop()

    def handle_play_stop(self):
        if self.play_timer.isActive():
            self.stop()
        else:
            self.toolbutton_play.setIcon(get_icon('stop.png'))
            self.play_timer.start(self.speed_slider.value())

    def handle_begining(self):
        self.change_frame(0)

    def handle_end(self):
        self.change_frame(max(len(self.vtu_files), len(self.vtp_files)))

    def forward(self):
        self.change_frame(self.frame_index + 1)

    def change_frame(self, index):

        # assume that whatever one is bigger has the smaller time step
        n_vtp = len(self.vtp_files)
        n_vtu = len(self.vtu_files)
        n_max = max(n_vtp, n_vtu)

        if n_max > 0:
            if index >= n_max:
                index = n_max-1
            elif index < 0:
                index = 0

            self.frame_spinbox.setValue(index)

            if index == self.frame_index:
                return
            else:
                self.frame_index = index

            if n_vtp > n_vtu:
                time = list(self.vtp_files.keys())[index]
                self.read_vtp(self.vtp_files[time])
                if n_vtu:
                    self.read_vtu(self.vtu_files.values()[bisect_left(self.vtu_files.keys(), time)-1])
            else:
                time = list(self.vtu_files.keys())[index]
                self.read_vtu(self.vtu_files[time])
                if n_vtp:
                    self.read_vtp(self.vtp_files.values()[bisect_left(self.vtp_files.keys(), time)-1])
            self.render()

            if self.checkbox_snap.isChecked():
                self.screenshot(True, fname=os.path.join(self.project_dir, self.project_name+'_'+str(index).zfill(4)+'.png'))
        else:
            self.frame_spinbox.setValue(0)

    def look_for_files(self):
        base_name = os.path.join(self.project_dir, self.project_name.upper())
        vtu_pvd = base_name + '.pvd'
        if os.path.exists(vtu_pvd):
            self.vtu_files = parse_pvd_file(vtu_pvd)
        else:
            self.vtu_files = build_time_dict(os.path.join(self.project_dir, '*.vtu'))
        self.vtp_files = parse_pvd_file(base_name + '_DES.pvd')

    # --- vtk functions ---
    def read_vtu(self, path):
        init = False
        if self.ugrid_cell_mapper is None:
            self.init_ugrid()
            init = True

        path = os.path.join(self.project_dir, path)
        self.ugrid_reader.SetFileName(path)
        self.ugrid_reader.Update()

        # TODO: Build Once
        data = self.ugrid_reader.GetOutput()
        cell_data = data.GetCellData()
        new_array_info = {}
        for i in range(cell_data.GetNumberOfArrays()):
            array = cell_data.GetArray(i)
            new_array_info[cell_data.GetArrayName(i)] = {
                'i':i,
                'components':array.GetNumberOfComponents(),
                'range': array.GetRange(),}

        self.cell_arrays = update(self.cell_arrays, copy.deepcopy(new_array_info))
        self.node_arrays = update(self.node_arrays, copy.deepcopy(new_array_info))

        if init:
            name = cell_data.GetArrayName(0)
            for t, m in [('cells', self.ugrid_cell_mapper), ('nodes', self.ugrid_node_mapper)]:
                m.SelectColorArray(name)
                combo = self.visual_btns[t]['color_by']
                combo.addItems(self.cell_arrays.keys())
                combo.setCurrentIndex(combo.findText(name))

    def read_vtp(self, path):
        init = False
        if self.particle_mapper is None:
            self.init_particles()
            init = True

        path = os.path.join(self.project_dir, path)
        self.particle_reader.SetFileName(path)
        self.particle_reader.Update()

        data = self.particle_reader.GetOutput()
        point_data = data.GetPointData()
        new_array_info = {}
        for i in range(point_data.GetNumberOfArrays()):
            array = point_data.GetArray(i)
            new_array_info[point_data.GetArrayName(i)] = {
                'i':i,
                'components':array.GetNumberOfComponents(),
                'range': array.GetRange()}

        self.point_arrays = update(self.point_arrays, copy.deepcopy(new_array_info))

        if 'Diameter' in self.point_arrays:
            self.glyph.SetScaleModeToScaleByScalar()
            self.glyph.SetInputArrayToProcess(0, 0, 0, 0, 'Diameter')
        if init:
            name = point_data.GetArrayName(0)
            self.glyph.SetInputArrayToProcess(1, 0, 0, 0, name)
            combo = self.visual_btns['points']['color_by']
            combo.addItems(self.point_arrays.keys())
            combo.setCurrentIndex(combo.findText(name))

    def change_visibility(self, geo, visible):
        if geo in self.actors:
            self.actors[geo].SetVisibility(visible)
            self.render()

    def change_color_by(self, geo, colorby, component=None):
        array_name = colorby.currentText()

        index = None
        if component is not None:
            comp = component.currentText()
            index = {'x':0, 'y':1, 'z':2, 'mag':None}[comp]

        if geo == 'points':
            self.glyph.SetInputArrayToProcess(1, 0, 0, 0, array_name)
#            self.particle_mapper.SelectColorArray(array_name)
            array = self.point_arrays[array_name]
#            self.glyph.SetRange(array.get('from', 0), array.get('to', 1))
        elif geo == 'cells':
            self.ugrid_cell_mapper.SelectColorArray(array_name)
            array = self.cell_arrays[array_name]
        elif geo == 'nodes':
            self.ugrid_node_mapper.SelectColorArray(array_name)
            array = self.node_arrays[array_name]
        else:
            return

        mapper = self.mappers.get(geo)
        mapper.SetScalarRange(array.get('from', 0), array.get('to', 1))

        single_color = array.get('single_color', False)
        if single_color:
            mapper.ScalarVisibilityOff()
            color = array.get('color', QtCore.Qt.white)
            btn = self.visual_btns[geo]['color']
            btn.setStyleSheet("QToolButton{{ background: {};}}".format(color.name()))
            btn.setIcon(QtGui.QIcon())
        else:
            self.change_color_bar(geo, array.get('color_map', 'viridis'), index)
            mapper.ScalarVisibilityOn()

        is_comp = array['components'] == 3 and not geo == 'points'
        self.visual_btns[geo]['component'].setEnabled(is_comp)

        self.render()

    def change_color_bar(self, geo, colormap, component=None):
        """change the color map"""
        mapper = self.mappers.get(geo, None)
        if mapper is None: return
        lut = self.lookuptables.get(geo, None)

        if colormap is not None:
            new_lut = vtk.vtkLookupTable()
            new_lut.DeepCopy(LOOKUP_TABLES.get(colormap, 'viridis'))
            # check component in old lut
            if lut is not None and component is None and lut.GetVectorMode() != 0:
                component = lut.GetVectorComponent()
        else:
            new_lut = lut

        if isinstance(component, int):
            new_lut.SetVectorModeToComponent()
            new_lut.SetVectorComponent(component)
        else:
            new_lut.SetVectorModeToMagnitude()

        mapper.SetLookupTable(new_lut)
        self.lookuptables[geo] = new_lut

        if colormap is not None:
            btn = self.visual_btns[geo]['color']
            btn.setIcon(build_qicons().get(colormap).get('icon', QtGui.QIcon()))
            btn.setStyleSheet("QToolButton{{ background: {};}}".format(None))

    def handle_change_color(self, geo, button, colorby):
        """popup the color bar dialog"""

        array_name = colorby.currentText()
        if len(array_name) == 0: return
        if geo == 'points':
            array = self.point_arrays[array_name]
        elif geo == 'cells':
            array = self.cell_arrays[array_name]
        elif geo == 'nodes':
            array = self.node_arrays[array_name]
        else:
            return

        self.color_dialog.set_(array)
        self.color_dialog.geo = geo
        self.color_dialog.button = button
        result = self.color_dialog.exec_()
        if result == QtWidgets.QDialog.Rejected:
            return

        self.change_color(geo, button, array)

    def change_color(self, geo, button, array):
        """change the color or color map of an actor"""
        mapper = self.mappers.get(geo)
        actor = self.actors.get(geo)

        params = self.color_dialog.get()

        array.update(params)

        color = params.get('color', QtCore.Qt.white)
        color_map = params.get('color_map', 'viridis')

        self.change_color_bar(geo, color_map)

        single_color = params.get('single_color', False)
        if single_color:
            button.setStyleSheet("QToolButton{{ background: {};}}".format(color.name()))
            actor.GetProperty().SetColor(color.getRgbF()[:3])
            button.setIcon(QtGui.QIcon())
            mapper.ScalarVisibilityOff()
        else:
            button.setStyleSheet("QToolButton{{ background: {};}}".format(None))
            mapper.ScalarVisibilityOn()

        mapper.SetScalarRange(array.get('from', 0), array.get('to', 1))

        self.render()

    def change_geo_color(self, button):
        """Change the color of the geometry actor"""
        col = QtWidgets.QColorDialog.getColor(parent=self)
        if not col.isValid():
            return

        button.setStyleSheet("QToolButton{{ background: {};}}".format(
            col.name()))

        self.actors['geometry'].GetProperty().SetColor(col.getRgbF()[:3])

        self.render()

    def change_opacity(self, geo, opacity):
        """change the opactiy of an actor"""
        if geo in self.actors:
            self.actors[geo].GetProperty().SetOpacity(opacity)
            self.render()


    def handle_particle_options(self):
        result = self.particle_option_dialog.exec_()
        if result == QtWidgets.QDialog.Rejected:
            return

        self.change_particle_options()

    def change_particle_options(self):
        data = self.particle_option_dialog.get()
        self.glyph_mask.SetMaximumNumberOfPoints(data.get('max_points', DEFAULT_MAXIMUM_POINTS))
        self.set_glyph_source(data.get('glyph', 'sphere'))
        self.render()

    def set_glyph_source(self, name):
        gs = GLYPHS[name]()
        self.glyph.SetSourceConnection(gs.GetOutputPort())