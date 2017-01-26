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
    DEFAULT_PEN = pg.mkPen(color='#64B5F6', width=2)
except ImportError:
    pg = None
    PYQTGRAPH_AVAILABLE = False
    DEFAULT_PEN = 'b'
except RuntimeError:
    pg = None
    PYQTGRAPH_AVAILABLE = False
    DEFAULT_PEN = 'b'

from qtpy import QtCore, QtWidgets
from mfixgui.tools.general import get_icon, clear_layout, get_unique_string

PLOT_ITEMS = OrderedDict([
    ['Select an item', {}],
    ['dt', {'left':'dt', 'bottom':'Time Step', 'var':'dt'}],
    ['nit', {'left':'Number of Iterations', 'bottom':'Time Step', 'var':'nit'}],
    ['time', {'left':'Simulation Time [s]', 'bottom':'Elapsed Wall Time [s]', 'var':'time', 'var2':'walltime_elapsed'}],
    ])

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

        if self.loadvtk and not 'MFIX_NO_VTK' in os.environ:
            from mfixgui.widgets.base_vtk import VTK_AVAILABLE
            plotbtn.setEnabled(VTK_AVAILABLE)
        else:
            plotbtn.setEnabled(False)

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
        from mfixgui.widgets.vtk_results_widget import GraphicsVtkWidget
        self.vtk_widget = GraphicsVtkWidget(self)
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
        for k, plot in self.plot_dict.items():
            if k in PLOT_ITEMS:
                if k not in status: continue
                props = PLOT_ITEMS[k]
                if 'var2' in props:
                    plot.plot(x=status[props['var2']], y=status[props['var']], append=True)
                else:
                    plot.plot(y=status[props['var']], append=True)
