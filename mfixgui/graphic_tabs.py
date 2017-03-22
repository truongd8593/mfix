# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict
import os

# graphics libraries
try:
    import pyqtgraph as pg
    PYQTGRAPH_AVAILABLE = True
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
except (ImportError, RuntimeError):
    pg = None
    PYQTGRAPH_AVAILABLE = False

from qtpy import QtCore, QtWidgets
from mfixgui.tools.general import get_icon, clear_layout, get_unique_string

PLOT_ITEMS = OrderedDict([
    ['Select an item', {}],
    ['DT vs Simulation Time', {
        'left':'Time Step [s]',
        'bottom':'Simulation Time [s]',
        'x_var':'time',
        'y_var':'dt'}],
    ['NIT vs Simulation time', {
        'left':'Number of Iterations [-]',
        'bottom':'Simulation Time [s]',
        'x_var':'time',
        'y_var':'nit'}],
    ['Elapsed time vs Simulation time', {
        'left':'Simulation Time [s]',
        'bottom':'Elapsed Wall Time [s]',
        'x_var':'walltime_elapsed',
        'y_var':'time'}],
#    ['Residuals vs Simulation time', {
#        'left':'Residual',
#        'bottom':'Simulation Time [s]',
#        'x_var':'time',
#        'y_var':'residuals'}],
    ])

TABLEAU20 = [(31, 119, 180),
             (255, 127, 14),
             (44, 160, 44),
             (214, 39, 40),
             (148, 103, 189),
             (140, 86, 75),
             (227, 119, 194),
             (127, 127, 127),
             (188, 189, 34),
             (219, 219, 141),
             (23, 190, 207),
             (158, 218, 229),
             (174, 199, 232),
             (255, 187, 120),
             (152, 223, 138),
             (255, 152, 150),
             (197, 176, 213),
             (196, 156, 148),
             (247, 182, 210),
             (199, 199, 199)
             ]

if PYQTGRAPH_AVAILABLE:
    DEFAULT_PENS = []
    for color in TABLEAU20:
        DEFAULT_PENS.append(pg.mkPen(color=color, width=2))
else:
    DEFAULT_PENS = []

class BaseGraphicTab(QtWidgets.QWidget):
    """Base graphics plot widget, provides options to redirect to a specific
    graphic widget, i.e. plot"""
    def __init__(self, plot_dict, name, parent=None, loadvtk=True):
        QtWidgets.QWidget.__init__(self, parent)

        self.tabWidget = parent
        self.loadvtk = loadvtk

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

        if self.loadvtk and 'MFIX_NO_VTK' not in os.environ:
            from mfixgui.widgets.base_vtk import VTK_AVAILABLE
            plotbtn.setEnabled(VTK_AVAILABLE)
        else:
            plotbtn.setEnabled(False)

    def create_plot_widget(self, combobox):
        name = combobox.currentText()

        if name == 'Select an item': return

        clear_layout(self.layout)
        self.change_name(name)
        props = PLOT_ITEMS.get(name, None)
        if props is None:
            return

        plot = self.plot_widget = pg.PlotWidget()
        plot.addLegend()
        plot.setLabel('left', props['left'])
        plot.setLabel('bottom', props['bottom'])
        plot.getPlotItem().showGrid(True, True, 0.5)
        plot.setDownsampling(ds=True, auto=True, mode='subsample')

        if 'residuals' in name.lower():
            plot.setLogMode(y=True)
            self.curve = {}
        else:
            self.curve = plot.plot([], pen=DEFAULT_PENS[0], name=name)
        self.layout.addWidget(plot)

    def create_vtk_widget(self):
        clear_layout(self.layout)
        self.change_name('VTK')
        from mfixgui.widgets.vtk_results_widget import GraphicsVtkWidget
        self.vtk_widget = GraphicsVtkWidget(self)
        self.layout.addWidget(self.vtk_widget)

    def plot(self, x=None, y=None, append=False):
        """update the plot"""
        if y is None: return

        # append the data to the internal data structure
        if append:
            if x is not None:
                self.x.append(x)
            if isinstance(y, (list, tuple)):
                if self.y == []:
                    self.y = {}
                for item in y:
                    if len(item) == 2:
                        key = item[0].strip()
                        if key:
                            self.y.setdefault(key,[]).append(float(item[1]))
            elif y is not None:
                self.y.append(y)

        # handle residual ploting
        if isinstance(self.y, dict):
            for i, (key, y) in enumerate(self.y.items()):
                curve = self.curve.get(key, None)

                # no curve, create one
                if curve is None:
                    curve = self.curve[key] = self.plot_widget.plot([], pen=DEFAULT_PENS[i], name=key)
                if len(y) > 1:
                    if len(self.x) == 0:
                        curve.setData(y)
                    elif len(self.x) > 1 and self.x[-1] - self.x[0] > 0:
                        curve.setData(self.x, y)
        # normal plotting
        elif len(self.y)>1:
            if len(self.x) == 0:
                self.curve.setData(self.y)
            elif len(self.x) > 1 and self.x[-1] - self.x[0] > 0:
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
    def init_graphic_tabs(self, loadvtk=True):
        self.loadvtk=loadvtk
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
        new_tab = BaseGraphicTab(self.plot_dict, name, tab_w, self.loadvtk)
        self.plot_dict[name] = new_tab
        index = tab_w.addTab(new_tab, name)
        new_tab.add_close_btn()
        self.ui.tabWidgetGraphics.setCurrentIndex(index)

    def update_plots(self, status):
        # loop through plots
        for k, plot in self.plot_dict.items():
            if k in PLOT_ITEMS:
                props = PLOT_ITEMS.get(k, None)
                if props is None: continue
                x_var = props.get('x_var', None)
                y_var = props.get('y_var', None)
                if y_var is None: continue
                plot.plot(x=status.get(x_var, None), y=status.get(y_var, None), append=True)
