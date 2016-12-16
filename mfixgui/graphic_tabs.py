# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict
from distutils.version import LooseVersion
import glob
import os

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
except:
    vtk = None
    VTK_AVAILABLE = False

from qtpy import QtCore, QtWidgets, QtGui
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

class GraphicsVtkWidget(BaseVtkWidget):
    """vtk widget for showing results"""
    def __init__(self, parent=None):
        BaseVtkWidget.__init__(self, parent)

        self.cell_arrays = {}
        self.frame_index = -1

        self.play_timer = QtCore.QTimer()
        self.play_timer.timeout.connect(self.forward)

        # look for vtu files
        self.project_dir = os.path.dirname(SETTINGS.value('project_file'))
        self.vtu_files = sorted(glob.glob(os.path.join(self.project_dir, '*.vtu')))


        self.init_vtk()
        self.init_toolbar()

        self.change_frame(0)
        self.reset_view()

    def init_vtk(self):

        # unstructured grid
        self.ugrid_reader = vtk.vtkXMLUnstructuredGridReader()

        self.ugrid_mapper = vtk.vtkDataSetMapper()
        self.ugrid_mapper.SetInputConnection(self.ugrid_reader.GetOutputPort())
        self.ugrid_mapper.SetScalarVisibility(True)
        self.ugrid_mapper.SetScalarModeToUseCellFieldData()
        self.ugrid_mapper.CreateDefaultLookupTable()

        self.ugrid_actor = vtk.vtkActor()
        self.ugrid_actor.SetMapper(self.ugrid_mapper)

        self.vtkrenderer.AddActor(self.ugrid_actor)

        # particles
        self.particle_reader = vtk.vtkXMLPolyDataReader()

        ss = vtk.vtkSphereSource()

        self.glyph = vtk.vtkGlyph3D()
        self.glyph.SetInputConnection(self.particle_reader.GetOutputPort())
        self.glyph.SetSourceConnection(ss.GetOutputPort())

        self.particle_mapper = vtk.vtkPolyDataMapper()
        self.particle_mapper.SetInputConnection(self.glyph.GetOutputPort())
        self.particle_mapper.CreateDefaultLookupTable()

        self.particle_actor = vtk.vtkActor()
        self.particle_actor.SetMapper(self.particle_mapper)

        self.vtkrenderer.AddActor(self.particle_actor)

    def init_toolbar(self):
        self.init_base_toolbar()

        # more buttons
        self.toolbutton_visible = QtWidgets.QToolButton()
        self.toolbutton_visible.setIcon(get_icon('visibility.png'))
        self.visible_menu = CustomPopUp(self, self.toolbutton_visible)
        self.visible_menu.finished.connect(lambda ignore: self.toolbutton_visible.setDown(False))
        self.toolbutton_visible.pressed.connect(self.visible_menu.popup)

        # --- visual representation menu ---
        layout = self.visible_menu.layout
        self.visual_btns = {}
        for i, geo in enumerate(['Cells', 'Nodes', 'Points']):
            geo_name = geo
            geo = geo.lower().replace(' ', '_')
            btns = self.visual_btns[geo] = {}
            # tool button
            toolbutton = QtWidgets.QToolButton(self.visible_menu)
            toolbutton.clicked.connect(lambda down, g=geo: self.change_visibility(geo, down))
            toolbutton.setCheckable(True)
            toolbutton.setChecked(True)
            toolbutton.setAutoRaise(True)
            toolbutton.setIcon(get_icon('visibility.png'))
            layout.addWidget(toolbutton, i, 0)
            btns['visible'] = toolbutton

            # color by
            combo = QtWidgets.QComboBox(self.visible_menu)
            combo.activated.connect(lambda item, g=geo: self.change_color_by(g, item))
            layout.addWidget(combo, i, 1)

            # color bar
            combo = QtWidgets.QComboBox(self.visible_menu)
            combo.activated.connect(lambda item, g=geo: self.change_color_bar(g, item))
            layout.addWidget(combo, i, 2)

            opacity = QtWidgets.QDoubleSpinBox(self.visible_menu)
            opacity.setRange(0, 1)
            opacity.setValue(1.0)
            opacity.setSingleStep(0.1)
            opacity.valueChanged.connect(lambda o, g=geo: self.change_opacity(geo, o))
            layout.addWidget(opacity, i, 3)
            btns['opacity'] = opacity

            # label
            label = QtWidgets.QLabel(geo_name, self.visible_menu)
            layout.addWidget(label, i, 4)

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

    def showEvent(self, event):
        # has to be called after the widget is visible
        self.vtkiren.Initialize()

    def hideEvent(self, event):
        self.stop()

    def close(self):
        BaseVtkWidget.close(self)

        # clean up timer
        self.play_timer.stop()

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
        self.change_frame(len(self.vtu_files))

    def forward(self):
        self.change_frame(self.frame_index + 1)

    def change_frame(self, index):
        # look for more files
        self.vtu_files = sorted(glob.glob(os.path.join(self.project_dir, '*.vtu')))
        n_files = len(self.vtu_files)
        if n_files> 0:
            if index >= n_files:
                index = n_files-1
            elif index < 0:
                index = 0

            if index == self.frame_index:
                self.frame_spinbox.setValue(index)
                return
            else:
                self.frame_index = index

            self.frame_spinbox.setValue(index)
            self.read_vtu(self.vtu_files[index])

        self.vtp_files = sorted(glob.glob(os.path.join(self.project_dir, '*.vtp')))
        n_files = len(self.vtp_files)
        if n_files> 0:
            if index >= n_files:
                index = n_files-1
            elif index < 0:
                index = 0

            if index == self.frame_index:
                self.frame_spinbox.setValue(index)
                return
            else:
                self.frame_index = index

            self.frame_spinbox.setValue(index)
            self.read_vtp(self.vtp_files[index])

        self.render()

    # --- vtk functions ---
    def read_vtu(self, path):

        self.ugrid_reader.SetFileName(path)
        self.ugrid_reader.Update()

        # TODO: Build Once
        data = self.ugrid_reader.GetOutput()
        cell_data = data.GetCellData()
        self.cell_arrays = {}
        for i in range(cell_data.GetNumberOfArrays()):
            self.cell_arrays[cell_data.GetArrayName(i)] = {'i':i, 'components':cell_data.GetArray(i).GetNumberOfComponents()}

        self.ugrid_mapper.ColorByArrayComponent(0, 0)

    def read_vtp(self, path):
        self.particle_reader.SetFileName(path)
        self.ugrid_reader.Update()

        data = self.ugrid_reader.GetOutput()
        point_data = data.GetPointData()
        self.point_arrays = {}
        for i in range(point_data.GetNumberOfArrays()):
            self.cell_arrays[point_data.GetArrayName(i)] = {'i':i, 'components':point_data.GetArray(i).GetNumberOfComponents()}

        self.glyph.SetScaleModeToScaleByScalar()
        self.glyph.SetColorModeToColorByVector()
        self.glyph.SetInputArrayToProcess(0, 0, 0, 0, 'Diameter')
        self.glyph.SetInputArrayToProcess(1, 0, 0, 0, 'Velocity')

        self.particle_mapper.SetScalarRange(0, 100)

    def change_visibility(self, geo, visible):
        print(geo, visible)

    def change_color_by(self, geo, item):
        print(geo, item)

    def change_color_bar(self, geo, item):
        print(geo, item)

    def change_opacity(self, geo, opacity):
        print(geo, opacity)


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
