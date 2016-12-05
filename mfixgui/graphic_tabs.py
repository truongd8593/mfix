# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division

# graphics libraries
try:
    import pyqtgraph as pg
    PYQTGRAPH_AVAILABLE = True
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
except ImportError:
    pg = None
    PYQTGRAPH_AVAILABLE = False

from qtpy import QtCore, QtWidgets

from mfixgui.tools.general import get_icon, clear_layout

class BaseGraphicTab(QtWidgets.QWidget):
    """Base graphics plot widget, provides options to redirect to a specific
    graphic widget, i.e. plot"""
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

        self.tabWidget = parent

        self.x = []
        self.y = []
        self.z = []

        # close button
        self.close_btn = QtWidgets.QToolButton()
        self.close_btn.setAutoRaise(True)
        self.close_btn.setIcon(get_icon('close.png'))
        self.close_btn.pressed.connect(self.close)
        self.close_btn.setIconSize(QtCore.QSize(10, 10))

        # plot selection btns
        self.layout = QtWidgets.QGridLayout(self)

        self.build_option_widgets()


    def build_option_widgets(self):
        if PYQTGRAPH_AVAILABLE:
            combobox = QtWidgets.QComboBox()
            combobox.addItems(['dt'])
            self.layout.addWidget(combobox, 1, 1)

            plotbtn = QtWidgets.QToolButton()
            plotbtn.setIcon(get_icon('timeline.png'))
            plotbtn.pressed.connect(lambda: self.create_plot_widget(combobox))
            self.layout.addWidget(plotbtn, 1, 2)

    def create_plot_widget(self, combobox):
        clear_layout(self.layout)
        plot = pg.PlotWidget()
        plot.addLegend()
        plot.setLabel('left', 'dt')
        plot.setLabel('bottom', 'Time Step')
        plot.getPlotItem().showGrid(True, True, 0.5)
        plot.setDownsampling(ds=True, auto=True, mode='subsample')
        self.curve = plot.plot([], pen='b', name='dt')
        self.layout.addWidget(plot)

    def plot(self, x=None, y=None, z=None, append=False):
        if append:
            self.x.append(x)
            self.y.append(y)
            self.z.append(z)
        self.curve.setData(self.x)

    def get_index(self):
        return self.tabWidget.indexOf(self)

    def add_close_btn(self):
        self.tabWidget.tabBar().setTabButton(self.get_index(), QtWidgets.QTabBar.RightSide, self.close_btn)

    def close(self):
        self.tabWidget.removeTab(self.get_index())


class GraphicTabs(object):
    """mixin to the gui.MfixGui class to handle plots etc."""
    def init_graphic_tabs(self):

        self.plot_dict = {}

        # Add corner widget to tabs
        corner_widget = QtWidgets.QWidget()
        corner_layout = QtWidgets.QHBoxLayout(corner_widget)
        corner_layout.setSpacing(0)
        corner_layout.setMargin(0)
        toolbutton_add_plot = QtWidgets.QToolButton()
        toolbutton_add_plot.setAutoRaise(True)
        toolbutton_add_plot.setIcon(get_icon('add.png'))
        toolbutton_add_plot.pressed.connect(self.handle_new_tab)
        corner_layout.addWidget(toolbutton_add_plot)

        # graphics tab widget
        self.ui.tabWidgetGraphics.setCornerWidget(corner_widget)

    def handle_new_tab(self, name='Plot'):
        """callback to add a new tab to tabWidgetGraphics"""
        tab_w = self.ui.tabWidgetGraphics
        new_tab = BaseGraphicTab(tab_w)
        self.plot_dict[name] = new_tab
        index = tab_w.addTab(new_tab, name)
        new_tab.add_close_btn()
        self.ui.tabWidgetGraphics.setCurrentIndex(index)
