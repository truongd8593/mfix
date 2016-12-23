# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict
from distutils.version import LooseVersion
import glob
import os
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
    from mfixgui.colormaps.color_maps import build_vtk_lookup_tables, get_color_map_pngs
    LOOKUP_TABLES = build_vtk_lookup_tables()
    print(LOOKUP_TABLES)
except:
    vtk = None
    VTK_AVAILABLE = False
    build_vtk_lookup_tables, get_color_map_pngs = None, None
    LOOKUP_TABLES = {}

from qtpy import QtCore, QtWidgets, QtGui, uic
from mfixgui.widgets.base_vtk import BaseVtkWidget, vtk, VTK_AVAILABLE
from mfixgui.widgets.base import CustomPopUp
from mfixgui.tools.general import get_icon, clear_layout, get_unique_string

PLOT_ITEMS = OrderedDict([
    ['Select an item', {}],
    ['dt', {'left':'dt', 'bottom':'Time Step', 'var':'dt'}],
    ['nit', {'left':'Number of Iterations', 'bottom':'Time Step', 'var':'nit'}],
    ['time', {'left':'Simulation Time [s]', 'bottom':'Ellapsed Wall Time [s]', 'var':'time', 'var2':'walltime_elapsed'}],
    ])

SETTINGS = QtCore.QSettings('MFIX', 'MFIX')
DEFAULT_PLAYBACK_SPEED = 100
DEFAULT_MAXIMUM_POINTS = 10000


def parse_pvd_file(fname):
    '''given a pvd file, return a dict of time (float) : file_name'''
    f_dict = OrderedDict()
    if not os.path.exists(fname): return f_dict
    tree = ElementTree.parse(fname)
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


class ColorMapPopUp(QtWidgets.QDialog):
    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self, parent)

        self.color = None

        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(thisdir, 'uifiles')
        ui = self.ui = uic.loadUi(os.path.join(uidir, 'color_map.ui'), self)

        self.setWindowTitle('Color Map')

        ui.toolbutton_select_color.clicked.connect(self.select_color)
        ui.lineedit_from.dtype = float
        ui.lineedit_to.dtype = float

        for name, path in get_color_map_pngs().items():
            ui.combobox_color_map.addItem(QtGui.QIcon(QtGui.QPixmap(path).scaled(100, 25, transformMode=QtCore.Qt.SmoothTransformation)), name)

        self.set_color_map('viridis')

    def set_(self, array):
        self.array = array
        self.set_color(self.array.get('color', QtCore.Qt.white))
        self.set_color_map(self.array.get('color_map', 'viridis'))

        d_range = [0, 1]
        self.ui.lineedit_from.updateValue(None,
            self.array.get('from', self.array.get('range', d_range)[0]))
        self.ui.lineedit_to.updateValue(None,
            self.array.get('to', self.array.get('range', d_range)[1]))

    def get(self):
        return {
            'color':        self.color,
            'single_color': self.ui.checkbox_single_color.isChecked(),
            'color_map':    self.ui.combobox_color_map.currentText(),
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
        index = self.ui.combobox_color_map.findText(color_map)
        self.ui.combobox_color_map.setCurrentIndex(index)


class GraphicsVtkWidget(BaseVtkWidget):
    """vtk widget for showing results"""
    def __init__(self, parent=None):
        BaseVtkWidget.__init__(self, parent)

        self.cell_arrays = {}
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

        self.init_vtk()
        self.init_toolbar()

        self.change_frame(0)
        self.reset_view()

    def init_vtk(self):

        self.actors = {}

        # unstructured grid
        self.ugrid_reader = vtk.vtkXMLUnstructuredGridReader()

        self.ugrid_mapper = vtk.vtkDataSetMapper()
        self.ugrid_mapper.SetInputConnection(self.ugrid_reader.GetOutputPort())
        self.ugrid_mapper.SetScalarVisibility(True)
        self.ugrid_mapper.SetScalarModeToUseCellFieldData()
        self.ugrid_mapper.SetLookupTable(LOOKUP_TABLES['viridis'])

        self.ugrid_actor = vtk.vtkActor()
        self.ugrid_actor.SetMapper(self.ugrid_mapper)
        self.actors['cells'] = self.ugrid_actor

        self.vtkrenderer.AddActor(self.ugrid_actor)

        # particles
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

        self.particle_mapper = vtk.vtkPolyDataMapper()
        self.particle_mapper.SetInputConnection(self.glyph.GetOutputPort())
        self.particle_mapper.SetLookupTable(LOOKUP_TABLES['viridis'])

        self.particle_actor = vtk.vtkActor()
        self.particle_actor.SetMapper(self.particle_mapper)
        self.actors['points'] = self.particle_actor

        self.vtkrenderer.AddActor(self.particle_actor)

    def init_toolbar(self):
        self.init_base_toolbar()

        # more buttons
        self.toolbutton_visible = QtWidgets.QToolButton()
        self.toolbutton_visible.setIcon(get_icon('visibility.png'))
        self.visible_menu = CustomPopUp(self, self.toolbutton_visible)
        self.visible_menu.finished.connect(lambda ignore: self.toolbutton_visible.setDown(False))
        self.toolbutton_visible.pressed.connect(self.show_visible_menu)

        # --- visual representation menu ---
        layout = self.visible_menu.layout
        self.visual_btns = {}
        for i, geo in enumerate(['Cells', 'Nodes', 'Points', 'Geometry']):
            geo_name = geo
            geo = geo.lower().replace(' ', '_')
            btns = self.visual_btns[geo] = {}
            # tool button
            toolbutton = QtWidgets.QToolButton(self.visible_menu)
            toolbutton.clicked.connect(lambda down, g=geo: self.change_visibility(g, down))
            toolbutton.setCheckable(True)
            toolbutton.setChecked(True)
            toolbutton.setAutoRaise(True)
            toolbutton.setIcon(get_icon('visibility.png'))
            layout.addWidget(toolbutton, i, 0)
            btns['visible'] = toolbutton

            combo = None
            if not geo == 'geometry':
                # color by
                combo = QtWidgets.QComboBox(self.visible_menu)
                combo.activated.connect(lambda item, g=geo, c=combo: self.change_color_by(g, c))
                layout.addWidget(combo, i, 1)
                btns['color_by'] = combo

                # component
                combo2 = QtWidgets.QComboBox(self.visible_menu)
                combo2.activated.connect(lambda item, g=geo, c=combo, c2=combo2: self.change_color_by(g, c, c2))
                combo2.addItems(['mag', 'x', 'y', 'z'])
                combo2.setEnabled(False)
                layout.addWidget(combo2, i, 2)
                btns['component'] = combo2


            toolbutton = QtWidgets.QToolButton()
            toolbutton.clicked.connect(lambda ignore, g=geo, t=toolbutton, c=combo: self.change_color(g, t, c))
            toolbutton.setAutoRaise(True)
            layout.addWidget(toolbutton, i, 3)
            btns['color'] = toolbutton

            opacity = QtWidgets.QDoubleSpinBox(self.visible_menu)
            opacity.setRange(0, 1)
            opacity.setValue(1.0)
            opacity.setSingleStep(0.1)
            opacity.valueChanged.connect(lambda o, g=geo: self.change_opacity(g, o))
            layout.addWidget(opacity, i, 4)
            btns['opacity'] = opacity

            # label
            label = QtWidgets.QLabel(geo_name, self.visible_menu)
            layout.addWidget(label, i, 10)

        self.toolbutton_back = QtWidgets.QToolButton()
        self.toolbutton_back.clicked.connect(self.handle_begining)
        self.toolbutton_back.setIcon(get_icon('previous.png'))

        self.toolbutton_play = QtWidgets.QToolButton()
        self.toolbutton_play.clicked.connect(self.handle_play_stop)
        self.toolbutton_play.setIcon(get_icon('play.png'))

        self.toolbutton_forward = QtWidgets.QToolButton()
        self.toolbutton_forward.clicked.connect(self.handle_end)
        self.toolbutton_forward.setIcon(get_icon('next.png'))

        self.frame_spinbox = QtWidgets.QSpinBox()
        self.frame_spinbox.valueChanged.connect(self.change_frame)
        self.frame_spinbox.setMaximum(99999)

        self.toolbutton_play_speed = QtWidgets.QToolButton()
        self.toolbutton_play_speed.setIcon(get_icon('speed.png'))

        self.speed_menu = CustomPopUp(self, self.toolbutton_play_speed)
        self.speed_menu.finished.connect(lambda ignore: self.toolbutton_play_speed.setDown(False))
        self.speed_slider = QtWidgets.QSlider()
        self.speed_slider.setRange(0, 1000)
        self.speed_slider.setValue(DEFAULT_PLAYBACK_SPEED)
        self.speed_slider.setOrientation(QtCore.Qt.Horizontal)
        self.speed_slider.valueChanged.connect(self.handle_speed_changed)
        self.speed_menu.layout.addWidget(self.speed_slider)
        self.toolbutton_play_speed.pressed.connect(self.speed_menu.popup)

        for btn in [self.toolbutton_visible, self.toolbutton_back,
                    self.toolbutton_play, self.toolbutton_forward,
                    self.frame_spinbox, self.toolbutton_play_speed]:
            self.button_bar_layout.addWidget(btn)
            if isinstance(btn, QtWidgets.QToolButton):
                btn.setAutoRaise(True)

        self.button_bar_layout.addStretch()

        self.color_dialog  = ColorMapPopUp()

    def showEvent(self, event):
        # has to be called after the widget is visible
        self.vtkiren.Initialize()

    def hideEvent(self, event):
        self.stop()

    def close(self):
        BaseVtkWidget.close(self)

        # clean up timer
        self.play_timer.stop()

    def show_visible_menu(self):
        # update comboboxes based on avaliable arrays

        for type_, array in [('cells', self.cell_arrays),
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
                time = self.vtp_files.keys()[index]
                self.read_vtp(self.vtp_files[time])
                if n_vtu:
                    self.read_vtu(self.vtu_files.values()[bisect_left(self.vtu_files.keys(), time)-1])
            else:
                time = self.vtu_files.keys()[index]
                self.read_vtu(self.vtu_files[time])
                if n_vtp:
                    self.read_vtp(self.vtp_files.values()[bisect_left(self.vtp_files.keys(), time)-1])
            self.render()
        else:
            self.frame_spinbox.setValue(0)

    def look_for_files(self):
        base_name = os.path.join(self.project_dir, self.project_name)
        self.vtu_files = parse_pvd_file(base_name + '.pvd')
        if not self.vtu_files:
            self.vtu_files = build_time_dict(os.path.join(self.project_dir, '*.vtu'))
        self.vtp_files = parse_pvd_file(base_name + '_DES.pvd')

    # --- vtk functions ---
    def read_vtu(self, path):

        path = os.path.join(self.project_dir, path)
        self.ugrid_reader.SetFileName(path)
        self.ugrid_reader.Update()

        # TODO: Build Once
        data = self.ugrid_reader.GetOutput()
        cell_data = data.GetCellData()
        self.cell_arrays = {}
        for i in range(cell_data.GetNumberOfArrays()):
            array = cell_data.GetArray(i)
            self.cell_arrays[cell_data.GetArrayName(i)] = {
                'i':i,
                'components':array.GetNumberOfComponents(),
                'range': array.GetRange(),}

    def read_vtp(self, path):
        path = os.path.join(self.project_dir, path)
        self.particle_reader.SetFileName(path)
        self.particle_reader.Update()

        data = self.particle_reader.GetOutput()
        point_data = data.GetPointData()
        self.point_arrays = {}
        for i in range(point_data.GetNumberOfArrays()):
            array = point_data.GetArray(i)
            self.point_arrays[point_data.GetArrayName(i)] = {
                'i':i,
                'components':array.GetNumberOfComponents(),
                'range': array.GetRange()}

        if 'Diameter' in self.point_arrays:
            self.glyph.SetScaleModeToScaleByScalar()
            self.glyph.SetInputArrayToProcess(0, 0, 0, 0, 'Diameter')


        self.particle_mapper.SetScalarRange(0, 100)

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
            array = self.point_arrays[array_name]
        elif geo == 'cells':
            self.ugrid_mapper.SelectColorArray(array_name)
            if index is not None:
                self.ugrid_mapper.SetArrayComponent(index)
            array = self.cell_arrays[array_name]
        else:
            return

        is_comp = array['components'] == 3
        self.visual_btns[geo]['component'].setEnabled(is_comp)

        self.render()

    def change_color(self, geo, button, colorby):

        array_name = colorby.currentText()
        if geo == 'points':
            array = self.point_arrays[array_name]
            mapper = self.particle_mapper
        elif geo == 'cells':
            array = self.cell_arrays[array_name]
            mapper = self.ugrid_mapper
        else:
            return

        self.color_dialog.set_(array)
        self.color_dialog.exec_()
        params = self.color_dialog.get()

        array.update(params)

        color = params['color']
        if isinstance(color, QtGui.QColor):
            button.setStyleSheet("QToolButton{{ background: {};}}".format(
                color.name()))

        mapper.SetLookupTable(LOOKUP_TABLES[params.get('color_map', 'viridis')])

        actor = None
        if geo == 'points':
            actor = self.particle_actor
        elif geo == 'cells':
            actor = self.ugrid_actor
        if actor is not None and color is not None:
            actor.GetProperty().SetColor(color.getRgbF()[:3])

        self.render()

    def change_opacity(self, geo, opacity):
        if geo in self.actors:
            self.actors[geo].GetProperty().SetOpacity(opacity)
            self.render()


class BaseGraphicTab(QtWidgets.QWidget):
    """Base graphics plot widget, provides options to redirect to a specific
    graphic widget, i.e. plot"""
    def __init__(self, plot_dict, name, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

        self.tabWidget = parent

        self.x = []
        self.y = []
        self.plot_dict = plot_dict
        self.name = name

        # close button
        self.close_btn = QtWidgets.QToolButton()
        self.close_btn.setAutoRaise(True)
        self.close_btn.setIcon(get_icon('close.png'))
        self.close_btn.clicked.connect(self.close)
        self.close_btn.setIconSize(QtCore.QSize(10, 10))

        # plot selection btns
        self.layout = QtWidgets.QGridLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.build_option_widgets()


    def build_option_widgets(self):

        # spacers
        # left
        spacer = QtWidgets.QSpacerItem(1000, 0, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum,)
        self.layout.addItem(spacer, 0, 0)
        # right
        spacer = QtWidgets.QSpacerItem(1000, 0, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum,)
        self.layout.addItem(spacer, 0, 100)
        # top
        spacer = QtWidgets.QSpacerItem(0, 1000, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Expanding,)
        self.layout.addItem(spacer, 0, 0)
        # bottom
        spacer = QtWidgets.QSpacerItem(0, 1000, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Expanding,)
        self.layout.addItem(spacer, 100, 0)

        # pyqtgraph
        combobox = QtWidgets.QComboBox()
        # TODO: add more/build dynamically
        combobox.addItems(PLOT_ITEMS.keys())
        # disable already used plots
        model = combobox.model()
        for i in range(combobox.count()):
            item = model.item(i)
            if item.text() in self.plot_dict:
                item.setFlags(item.flags() & ~(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled))


        self.layout.addWidget(combobox, 1, 1)

        plotbtn = QtWidgets.QToolButton()
        plotbtn.setText('Plot')
        plotbtn.setSizePolicy(QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding))
        plotbtn.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon)
        plotbtn.setIcon(get_icon('timeline.png'))
        plotbtn.setIconSize(QtCore.QSize(24, 24))
        plotbtn.clicked.connect(lambda: self.create_plot_widget(combobox))
        self.layout.addWidget(plotbtn, 1, 2)

        combobox.setEnabled(PYQTGRAPH_AVAILABLE)
        plotbtn.setEnabled(PYQTGRAPH_AVAILABLE)

        # vtk
        plotbtn = QtWidgets.QToolButton()
        plotbtn.setText('VTK')
        plotbtn.setSizePolicy(QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding))
        plotbtn.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon)
        plotbtn.setIcon(get_icon('geometry.png'))
        plotbtn.setIconSize(QtCore.QSize(24, 24))
        plotbtn.clicked.connect(self.create_vtk_widget)
        self.layout.addWidget(plotbtn, 2, 2)

        plotbtn.setEnabled(VTK_AVAILABLE)


    def create_plot_widget(self, combobox):
        name = combobox.currentText()

        if name == 'Select an item': return

        clear_layout(self.layout)
        self.change_name(name)
        props = PLOT_ITEMS[name]

        plot = pg.PlotWidget()
        plot.addLegend()
        plot.setLabel('left', props['left'])
        plot.setLabel('bottom', props['bottom'])
        plot.getPlotItem().showGrid(True, True, 0.5)
        plot.setDownsampling(ds=True, auto=True, mode='subsample')

        self.curve = plot.plot([], pen=DEFAULT_PEN, name=name)
        self.layout.addWidget(plot)

    def create_vtk_widget(self):
        clear_layout(self.layout)
        self.change_name('VTK')
        self.vtk_widget = GraphicsVtkWidget()
        self.layout.addWidget(self.vtk_widget)

    def plot(self, x=None, y=None, append=False):
        if append:
            if x is not None:
                self.x.append(x)
            if y is not None:
                self.y.append(y)
        if len(self.y)>1:
            if len(self.x) == 0:
                self.curve.setData(self.y)
            else:
                self.curve.setData(self.x, self.y)

    def get_index(self):
        return self.tabWidget.indexOf(self)

    def add_close_btn(self):
        self.tabWidget.tabBar().setTabButton(self.get_index(), QtWidgets.QTabBar.RightSide, self.close_btn)

    def close(self):
        self.tabWidget.removeTab(self.get_index())
        if self.name in self.plot_dict:
            self.plot_dict.pop(self.name)

        if hasattr(self, 'vtk_widget'):
            self.vtk_widget.close()

    def change_name(self, name):
        self.tabWidget.tabBar().setTabText(self.get_index(), name)

        if self.name in self.plot_dict:
            self.plot_dict[name] = self.plot_dict.pop(self.name)

        self.name = name

class GraphicTabs(object):
    """mixin to the gui.MfixGui class to handle plots and graphics"""
    def init_graphic_tabs(self):

        self.plot_dict = {}

        # Add corner widget to tabs
        corner_widget = QtWidgets.QWidget()
        corner_layout = QtWidgets.QHBoxLayout(corner_widget)
        corner_layout.setContentsMargins(0,0,0,0)
        toolbutton_add_plot = QtWidgets.QToolButton()
        toolbutton_add_plot.setAutoRaise(True)
        toolbutton_add_plot.setIcon(get_icon('add.png'))
        toolbutton_add_plot.clicked.connect(self.handle_new_tab)
        corner_layout.addWidget(toolbutton_add_plot)
        corner_layout.addStretch() #TODO: this doesn't work

        # graphics tab widget
        self.ui.tabWidgetGraphics.setCornerWidget(corner_widget)

    def handle_new_tab(self, checked, name='New'):
        """callback to add a new tab to tabWidgetGraphics"""

        name = get_unique_string(name, self.plot_dict.keys())

        tab_w = self.ui.tabWidgetGraphics
        new_tab = BaseGraphicTab(self.plot_dict, name, tab_w)
        self.plot_dict[name] = new_tab
        index = tab_w.addTab(new_tab, name)
        new_tab.add_close_btn()
        self.ui.tabWidgetGraphics.setCurrentIndex(index)

    def update_plots(self, status):
        for k, plot in self.plot_dict.items():
            if k in PLOT_ITEMS:
                if k not in status: continue
                props = PLOT_ITEMS[k]
                if 'var2' in props:
                    plot.plot(x=status[props['var2']], y=status[props['var']], append=True)
                else:
                    plot.plot(y=status[props['var']], append=True)
