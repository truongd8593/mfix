""" MFIX GUI """


# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals
import glob
import logging
import os
import shutil
import signal
import subprocess
import sys
import time
from collections import OrderedDict
import copy

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen

# import qt
from qtpy import QtCore, QtWidgets, QtGui
from qtpy.QtCore import (QObject, QThread, pyqtSignal, QUrl, QTimer, QSettings,
                         Qt)

# TODO: add pyside?
try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic


# local imports
from widgets.vtkwidget import VtkWidget
from widgets.base import (LineEdit, CheckBox, ComboBox, SpinBox, DoubleSpinBox,
                          Table)
from tools.mfixproject import Project, Keyword
from tools.general import (make_callback, get_icon, get_unique_string,
                           widget_iter, set_script_directory)
from tools.namelistparser import buildKeywordDoc

# look for pyqtnode
try:
    from pyqtnode import NodeWidget
except ImportError:
    NodeWidget = None

logging.basicConfig(stream=sys.stdout,
                    filemode='w', level=logging.INFO,
                    format='%(name)s - %(levelname)s - %(message)s')

SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), ))
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))
LOG = logging.getLogger(__name__)
LOG.debug(SCRIPT_DIRECTORY)
set_script_directory(SCRIPT_DIRECTORY)  # should this be in an __init__.py?

def get_mfix_home():
    " return the top level directory (since __file__ is mfix_home/gui/gui.py) "
    return os.path.dirname(
        os.path.dirname(os.path.realpath(__file__)))

# --- Main Gui ---
class MfixGui(QtWidgets.QMainWindow):
    '''
    Main window class handling all gui interactions
    '''
    def __init__(self, app, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)

        # reference to qapp instance
        self.app = app

        # load ui file
        self.customWidgets = {'LineEdit':      LineEdit,
                              'CheckBox':      CheckBox,
                              'ComboBox':      ComboBox,
                              'DoubleSpinBox': DoubleSpinBox,
                              'SpinBox':       SpinBox,
                              'Table':         Table,
                              }

        self.ui = uic.loadUi(os.path.join('uifiles', 'gui.ui'), self)

        self.ui.general = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'general.ui'), self.ui.general)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.general)

        self.ui.geometry = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'geometry.ui'), self.ui.geometry)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.geometry)

        self.ui.mesh = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'mesh.ui'), self.ui.mesh)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.mesh)

        self.ui.regions = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'regions.ui'), self.ui.regions)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.regions)

        self.ui.model_setup = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'model_setup.ui'), self.ui.model_setup)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.model_setup)

        self.ui.numerics = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'numerics.ui'), self.ui.numerics)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.numerics)

        self.ui.output = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'output.ui'), self.ui.output)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.output)

        self.ui.monitors = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'monitors.ui'), self.ui.monitors)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.monitors)

        self.ui.run = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'run.ui'), self.ui.run)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.run)

        self.ui.interact = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'interact.ui'), self.ui.interact)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.interact)

        self.ui.post_processing = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'post_processing.ui'), self.ui.post_processing)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.post_processing)

        # load settings
        self.settings = QSettings('MFIX', 'MFIX')

        # set title and icon
        self.setWindowTitle('MFIX')
        self.setWindowIcon(get_icon('mfix.png'))

        # build keyword documentation from namelist docstrings
        self.keyword_doc = buildKeywordDoc(os.path.join(SCRIPT_DIRECTORY,
                                                        os.pardir, 'model'))

        # create project manager
        self.project = ProjectManager(self, self.keyword_doc)

        # --- data ---
        self.modebuttondict = {'modeler':   self.ui.pushButtonModeler,
                               'workflow':  self.ui.pushButtonWorkflow,
                               'developer': self.ui.pushButtonDeveloper,
                               }

        self.animation_speed = 400
        self.animating = False
        self.stack_animation = None

        # --- icons ---
        # loop through all widgets, because I am lazy
        for widget in widget_iter(self.ui):
            if isinstance(widget, QtWidgets.QToolButton):
                name = str(widget.objectName())
                if 'add' in name:
                    widget.setIcon(get_icon('add.png'))
                elif 'delete' in name or 'remove' in name:
                    widget.setIcon(get_icon('remove.png'))
                elif 'copy' in name:
                    widget.setIcon(get_icon('copy.png'))

        self.ui.toolbutton_new.setIcon(get_icon('newfolder.png'))
        self.ui.toolbutton_open.setIcon(get_icon('openfolder.png'))
        self.ui.toolbutton_save.setIcon(get_icon('save.png'))

        self.ui.toolbutton_compile.setIcon(get_icon('build.png'))
        self.ui.toolbutton_run.setIcon(get_icon('play.png'))
        self.ui.toolbutton_restart.setIcon(get_icon('restart.png'))
        self.ui.toolbutton_interact.setIcon(get_icon('flash.png'))

        self.ui.geometry.toolbutton_add_geometry.setIcon(get_icon('geometry.png'))
        self.ui.geometry.toolbutton_add_filter.setIcon(get_icon('filter.png'))
        self.ui.geometry.toolbutton_geometry_union.setIcon(get_icon('union.png'))
        self.ui.geometry.toolbutton_geometry_intersect.setIcon(
            get_icon('intersect.png'))
        self.ui.geometry.toolbutton_geometry_difference.setIcon(
            get_icon('difference.png'))

        self.ui.toolButtonTFMSolidsDatabase.setIcon(get_icon('download.png'))

        # --- Connect Signals to Slots---
        # open/save/new project
        self.ui.toolbutton_open.clicked.connect(self.open_project)
        self.ui.toolbutton_save.clicked.connect(self.save_project)

        # mode (modeler, workflow, developer)
        for mode, btn in self.modebuttondict.items():
            btn.released.connect(make_callback(self.mode_changed, mode))

        # navigation tree
        self.ui.treewidget_model_navigation.itemSelectionChanged.connect(
            self.navigation_changed)

        # build/run/connect MFIX
        self.ui.run.run_mfix_button.clicked.connect(self.run_mfix)
        self.ui.run.connect_mfix_button.clicked.connect(self.connect_mfix)
        self.ui.run.clear_output_button.clicked.connect(self.clear_output)
        self.ui.run.mfix_executables.activated.connect(self.update_run)

        # --- Threads ---
        self.run_thread = MfixThread(self)
        self.clear_thread = MfixThread(self)
        self.monitor_executables_thread = MonitorExecutablesThread(self)

        def make_handler(qtextbrowser):
            " make a closure to read output from external process "

            def handle_line(line, color=None):
                " closure to read output from external process "

                log = logging.getLogger(__name__)
                log.debug(str(line).strip())
                cursor = qtextbrowser.textCursor()
                cursor.movePosition(cursor.End)
                char_format = QtGui.QTextCharFormat()
                if color:
                    char_format.setForeground(QtGui.QColor(color))
                cursor.setCharFormat(char_format)
                cursor.insertText(line)
                scrollbar = qtextbrowser.verticalScrollBar()
                scrollbar.setValue(scrollbar.maximum())

            return handle_line

        def update_run_executables(*args):
            "Updates the list of executables available"
            self.ui.run.mfix_executables.clear()
            output = self.monitor_executables_thread.get_output()
            for executable in output:
                self.ui.run.mfix_executables.addItem(executable)
            mfix_available = bool(output)
            self.ui.run.mfix_executables.setVisible(mfix_available)
            self.ui.run.mfix_executables_warning.setVisible(not mfix_available)

        self.run_thread.line_printed.connect(
            make_handler(self.ui.command_output))
        self.clear_thread.line_printed.connect(
            make_handler(self.ui.command_output))
        update_run_executables()

        self.monitor_executables_thread.sig.connect(update_run_executables)
        self.monitor_executables_thread.start()

        # --- setup widgets ---
        self.__setup_simple_keyword_widgets()
        self.__setup_other_widgets()

        self.__setup_regions()

        # --- vtk setup ---
        self.__setup_vtk_widget()

        # --- workflow setup ---
        if NodeWidget is not None:
            self.__setup_workflow_widget()

        # --- default ---
        self.mode_changed('modeler')
        self.change_pane('geometry')

        # autoload last project
        if self.get_project_dir():
            self.open_project(self.get_project_dir())

        # print number of keywords
        self.print_internal('Registered {} keywords'.format(
            len(self.project.registered_keywords)))

    def set_navigation_item_state(self, item_name, state):
        on = Qt.ItemIsSelectable | Qt.ItemIsEnabled
        off = Qt.ItemFlags(0)
        tree = self.ui.treewidget_model_navigation
        flags =  Qt.MatchFixedString | Qt.MatchRecursive
        items = tree.findItems(item_name, flags, 0)
        assert len(items) == 1, "Multiple menu items matching %s"%item_name
        items[0].setFlags(on if state else off)

    def set_solver(self, index):
        """ handler for "Solver" combobox in Model Setup """
        # NB.  This depends on strings matching items in the ui file
        model_setup = self.ui.model_setup
        solver_name = model_setup.combobox_solver.currentText()
        self.print_internal("set solver to %d: %s" % (index, solver_name))

        item_names =  ("Solids", "Continuum Solids Model",
                       "Discrete Element Model", "Particle in Cell Model")

        states = {"Single phase": (False, False, False, False),
                  "MFIX-TFM": (True, True, False, False),
                  "MFIX-DEM": (True, False, True, False),
                  "MFIX-PIC": (True, False, False, True),
                  "MFIX-Hybrid": (True, True, True, False)}

        for item_name, state in zip(item_names, states[solver_name]):
            self.set_navigation_item_state(item_name, state)

        # Options which require TFM, DEM, or PIC
        enabled = 0 < index < 4
        interphase = model_setup.groupbox_interphase
        interphase.setEnabled(enabled)

        # TFM only
        # use a groupbox here, instead of accessing combox + label?
        enabled = (index == 1)
        model_setup.combobox_subgrid_model.setEnabled(enabled)
        model_setup.label_subgrid_model.setEnabled(enabled)
        model_setup.groupbox_subgrid_params.setEnabled(enabled and
                                                       self.subgrid_model > 0)



    def disable_fluid_solver(self, state):
        self.set_navigation_item_state("Fluid", not state)

    def set_subgrid_model(self, index):
        self.subgrid_model = index
        groupbox_subgrid_params = self.ui.model_setup.groupbox_subgrid_params
        groupbox_subgrid_params.setEnabled(index > 0)

    def __setup_other_widgets(self):
        """ setup widgets which are not tied to a simple keyword
        """
        # move to another file - cgw
        model_setup = self.ui.model_setup
        combobox_solver = model_setup.combobox_solver
        combobox_solver.currentIndexChanged.connect(self.set_solver)
        self.set_solver(0) # Default - Single Phase - (?)

        checkbox_disable_fluid_solver = model_setup.checkbox_disable_fluid_solver
        checkbox_disable_fluid_solver.stateChanged.connect(self.disable_fluid_solver)
        self.disable_fluid_solver(False)

        combobox_subgrid_model = model_setup.combobox_subgrid_model
        combobox_subgrid_model.currentIndexChanged.connect(self.set_subgrid_model)
        self.set_subgrid_model(0)

    def __setup_simple_keyword_widgets(self):
        """
        Look for and connect simple keyword widgets to the project manager.
        Keyword information from the namelist doc strings is added to each
        keyword widget. The widget must be named: *_keyword_<keyword> where
        <keyword> is the actual keyword. For example:
        lineedit_keyword_run_name
        """

        # loop through all widgets looking for *_keyword_<keyword>
        for widget in widget_iter(self.ui):
            name_list = str(widget.objectName()).lower().split('_')

            if 'keyword' in name_list: # perf - cgw
                keyword = '_'.join(name_list[name_list.index('keyword')+1:])

                # set the key attribute to the keyword
                widget.key = keyword

                # add info from keyword documentation
                if keyword in self.keyword_doc:
                    doc = self.keyword_doc[keyword]
                    widget.setdtype(doc['dtype'])
                    if 'required' in doc:
                        widget.setValInfo(
                            req=doc['required'] == 'true')
                    if 'validrange' in doc:
                        if 'max' in doc['validrange']:
                            widget.setValInfo(
                                _max = doc['validrange']['max'])
                        if 'min' in doc['validrange']:
                            widget.setValInfo(
                                _min = doc['validrange']['min'])

                    if 'initpython' in doc:
                        widget.default(
                            self.keyword_doc[keyword]['initpython'])

                    if isinstance(widget, QtWidgets.QComboBox) and widget.count() < 1:
                            widget.addItems(list(doc['valids'].keys()))

                # register the widget with the project manager
                self.project.register_widget(widget, [keyword])

                # connect to unsaved method
                widget.value_updated.connect(self.unsaved)

    def __setup_vtk_widget(self):
        " setup the vtk widget "

        self.vtkwidget = VtkWidget(self.project, parent=self)
        self.ui.horizontalLayoutModelGraphics.addWidget(self.vtkwidget)

        # register with project manager
        self.project.register_widget(self.vtkwidget,
                                     ['xmin', 'xlength', 'ymin', 'ylength',
                                      'zmin', 'zlength', 'imax', 'jmax',
                                      'kmax'])

    def __setup_workflow_widget(self):

        self.nodeChart = NodeWidget(showtoolbar=False)
        # Build default node library
        self.nodeChart.nodeLibrary.buildDefaultLibrary()
        self.ui.horizontalLayoutPyqtnode.addWidget(self.nodeChart)

    def __setup_regions(self):
        " setup the region connections etc."

        self.ui.regions.combobox_regions_shape.addItems(['box', 'sphere',
                                                         'point'])

        self.ui.regions.toolbutton_region_add.pressed.connect(
            self.new_region)
        self.ui.regions.toolbutton_region_delete.pressed.connect(
            self.delete_region)
        self.ui.regions.toolbutton_region_copy.pressed.connect(
            self.copy_region)

        self.ui.regions.tablewidget_regions.dtype = OrderedDict
        self.ui.regions.tablewidget_regions._setModel()
        self.ui.regions.tablewidget_regions.set_columns(['shape', 'from',
                                                         'to'])
        self.ui.regions.tablewidget_regions.show_vertical_header(True)
        self.ui.regions.tablewidget_regions.set_value(OrderedDict())
        self.ui.regions.tablewidget_regions.auto_update_rows(True)
        self.ui.regions.tablewidget_regions.set_selection_model('cell',
                                                                multi=False)

        self.ui.regions.tablewidget_regions.new_selection.connect(
            self.update_region_parameters)

        for widget in widget_iter(self.ui.regions.groupbox_region_parameters):
            if hasattr(widget, 'value_updated'):
                widget.value_updated.connect(self.region_value_changed)

                # lineedit_regions_to_x
                name = str(widget.objectName())
                if '_to_' in name:
                    widget.key = '_'.join(name.split('_')[-2:])
                    widget.dtype = float
                elif '_from_' in str(widget.objectName()):
                    widget.key = '_'.join(name.split('_')[-2:])
                    widget.dtype = float
                elif 'name' in name:
                    widget.key = 'name'
                    widget.dtype = str
                elif 'shape' in name:
                    widget.key = 'shape'
                    widget.dtype = str

    def get_project_dir(self):
        " get the current project directory "

        last_dir = self.settings.value('project_dir')
        if last_dir:
            return last_dir
        else:
            return None

    def mode_changed(self, mode):
        " change the Modeler, Workflow, Developer tab"

        current_index = 0
        for i in range(self.ui.stackedwidget_mode.count()):
            widget = self.ui.stackedwidget_mode.widget(i)
            if mode == str(widget.objectName()):
                current_index = i
                break

        for key, btn in self.modebuttondict.items():
            btn.setChecked(mode == key)

        self.animate_stacked_widget(self.ui.stackedwidget_mode,
                                    self.ui.stackedwidget_mode.currentIndex(),
                                    current_index,
                                    'horizontal')

    # --- modeler pane navigation ---
    def change_pane(self, name):
        """ change to the specified pane """

        clist = self.ui.treewidget_model_navigation.findItems(
                    name,
                    Qt.MatchContains | Qt.MatchRecursive,
                    0)

        for item in clist:
            if str(item.text(0)).lower() == name.lower():
                break

        self.ui.treewidget_model_navigation.setCurrentItem(item)
        self.navigation_changed()

    def navigation_changed(self):
        """ an item in the tree was selected, change panes """
        current_selection = self.ui.treewidget_model_navigation.selectedItems()

        if current_selection:

            text = str(current_selection[-1].text(0))
            text = '_'.join(text.lower().split(' '))

            current_index = 0
            for i in range(self.ui.stackedWidgetTaskPane.count()):
                widget = self.ui.stackedWidgetTaskPane.widget(i)
                if text == str(widget.objectName()):
                    current_index = i
                    break

            self.animate_stacked_widget(
                self.ui.stackedWidgetTaskPane,
                self.ui.stackedWidgetTaskPane.currentIndex(),
                current_index)

    # --- animation methods ---
    def animate_stacked_widget(self, stackedwidget, from_, to,
                               direction='vertical', line=None, to_btn=None,
                               btn_layout=None):
        """ animate changing of qstackedwidget """

        # check to see if already animating
        if self.animating and self.stack_animation is not None:
            self.stack_animation.stop()

        from_widget = stackedwidget.widget(from_)
        to_widget = stackedwidget.widget(to)

        # get from geometry
        width = from_widget.frameGeometry().width()
        height = from_widget.frameGeometry().height()

        # offset
        # bottom to top
        if direction == 'vertical' and from_ < to:
            offsetx = 0
            offsety = height
        # top to bottom
        elif direction == 'vertical' and from_ > to:
            offsetx = 0
            offsety = -height
        elif direction == 'horizontal' and from_ < to:
            offsetx = width
            offsety = 0
        elif direction == 'horizontal' and from_ > to:
            offsetx = -width
            offsety = 0
        else:
            return

        # move to widget and show
        # set the geometry of the next widget
        to_widget.setGeometry(0 + offsetx, 0 + offsety, width, height)
        to_widget.show()
        to_widget.raise_()

        # animate
        # from widget
        animnow = QtCore.QPropertyAnimation(from_widget, "pos")
        animnow.setDuration(self.animation_speed)
        animnow.setEasingCurve(QtCore.QEasingCurve.InOutQuint)
        animnow.setStartValue(
            QtCore.QPoint(0,
                          0))
        animnow.setEndValue(
            QtCore.QPoint(0 - offsetx,
                          0 - offsety))

        # to widget
        animnext = QtCore.QPropertyAnimation(to_widget, "pos")
        animnext.setDuration(self.animation_speed)
        animnext.setEasingCurve(QtCore.QEasingCurve.InOutQuint)
        animnext.setStartValue(
            QtCore.QPoint(0 + offsetx,
                          0 + offsety))
        animnext.setEndValue(
            QtCore.QPoint(0,
                          0))

        # line
        animline = None
        if line is not None and to_btn is not None:
            animline = QtCore.QPropertyAnimation(line, "pos")
            animline.setDuration(self.animation_speed)
            animline.setEasingCurve(QtCore.QEasingCurve.InOutQuint)
            animline.setStartValue(
                QtCore.QPoint(line.geometry().x(),
                              line.geometry().y()))
            animline.setEndValue(
                QtCore.QPoint(to_btn.geometry().x(),
                              line.geometry().y()))

        # animation group
        self.stack_animation = QtCore.QParallelAnimationGroup()
        self.stack_animation.addAnimation(animnow)
        self.stack_animation.addAnimation(animnext)
        if animline is not None:
            self.stack_animation.addAnimation(animline)
        self.stack_animation.finished.connect(
            make_callback(self.animate_stacked_widget_finished,
                          stackedwidget, from_, to, btn_layout, line)
            )
        self.stack_animation.stateChanged.connect(
            make_callback(self.animate_stacked_widget_finished,
                          stackedwidget, from_, to, btn_layout, line)
            )

        self.animating = True
        self.stack_animation.start()

    def animate_stacked_widget_finished(self, widget, from_, to,
                                        btn_layout=None, line=None):
        """ cleanup after animation """
        if self.stack_animation.state() == QtCore.QAbstractAnimation.Stopped:
            widget.setCurrentIndex(to)
            from_widget = widget.widget(from_)
            from_widget.hide()
            from_widget.move(0, 0)

            if btn_layout is not None and line is not None:
                btn_layout.addItem(btn_layout.takeAt(
                    btn_layout.indexOf(line)), 1, to)

            self.animating = False

    # --- helper methods ---
    def message(self,
                title='Warning',
                icon='warning',
                text='This is a warning.',
                buttons=['ok'],
                default='ok',
                infoText=None,
                detailedtext=None,
                ):
        '''
        Create a message box:
        title = 'title'
        icon = 'warning' or 'info'
        text = 'test to show'
        buttons = ['ok',...] where value is 'ok', 'yes', 'no', 'cancel',
            'discard'
        default = 'ok' the default selected button
        infotext = 'extended information text'
        detailedtext = 'Some details'
        '''
        msgBox = QtWidgets.QMessageBox(self)
        msgBox.setWindowTitle(title)

        # Icon
        if icon == 'warning':
            icon = QtWidgets.QMessageBox.Warning
        else:
            icon = QtWidgets.QMessageBox.Information

        msgBox.setIcon(icon)

        # Text
        msgBox.setText(text)

        if infoText:
            msgBox.setInformativeText(infoText)

        if detailedtext:
            msgBox.setDetailedText(detailedtext)

        # buttons
        qbuttonDict = {'ok':      QtWidgets.QMessageBox.Ok,
                       'yes':     QtWidgets.QMessageBox.Yes,
                       'no':      QtWidgets.QMessageBox.No,
                       'cancel':  QtWidgets.QMessageBox.Cancel,
                       'discard': QtWidgets.QMessageBox.Discard,
                       }
        for button in buttons:
            msgBox.addButton(qbuttonDict[button])

            if button == default:
                msgBox.setDefaultButton(qbuttonDict[button])

        ret = msgBox.exec_()

        for key, value in qbuttonDict.items():
            if value == ret:
                result = key
                break

        return result

    def print_internal(self, text):
        if not text.endswith('\n'):
            text += '\n'
        LOG.info(str(text).strip())
        cursor = self.ui.command_output.textCursor()
        cursor.movePosition(cursor.End)
        cursor.insertText(text)
        cursor.movePosition(cursor.End)
        vbar = self.ui.command_output.verticalScrollBar()
        vbar.triggerAction(QtWidgets.QAbstractSlider.SliderToMaximum)
        self.ui.command_output.ensureCursorVisible()

    def update_run(self):
        """Enable/disable run options based on selected executable """
        mfix_exe = self.ui.run.mfix_executables.currentText()
        if not mfix_exe:
            return
        config = self.monitor_executables_thread.output[mfix_exe]
        smp_enabled = 'smp' in config
        dmp_enabled = 'dmp' in config
        pymfix_enabled = 'pymfix' in mfix_exe
        self.ui.run.openmp_threads.setEnabled(smp_enabled)
        if not dmp_enabled:
            self.ui.run.spinbox_keyword_nodesi.setValue(1)
            self.ui.run.spinbox_keyword_nodesj.setValue(1)
            self.ui.run.spinbox_keyword_nodesk.setValue(1)
        self.ui.run.spinbox_keyword_nodesi.setEnabled(dmp_enabled)
        self.ui.run.spinbox_keyword_nodesj.setEnabled(dmp_enabled)
        self.ui.run.spinbox_keyword_nodesk.setEnabled(dmp_enabled)
        self.ui.run.mfix_host.setEnabled(pymfix_enabled)
        self.ui.run.mfix_port.setEnabled(pymfix_enabled)
        self.ui.run.connect_mfix_button.setEnabled(pymfix_enabled)

    def clear_output(self):
        self.run_thread.start_command(
            cmd='rm -v -f *.LOG *.OUT *.RES *.SP? *.pvd *vtp',
            cwd=self.get_project_dir())

    def run_mfix(self):
        mfix_exe = self.ui.run.mfix_executables.currentText()
        config = self.monitor_executables_thread.output[mfix_exe]

        if mfix_exe.endswith('pymfix'):
            # run with python or python3, depending on sys.executable
            mfix_exe = '{} {}'.format(sys.executable, mfix_exe)

        if 'dmp' in config:
            nodesi = int(self.ui.nodesi.text()) # FIXME this is not right
            nodesj = int(self.ui.nodesj.text())
            nodesk = int(self.ui.nodesk.text())
            total = nodesi * nodesj * nodesk
            # FIXME: maybe we should save NODES* keywords to runname.mfx instead of passing them on command line?
            mfix_exe = 'mpirun -np {} {} NODESI={} NODESJ={} NODESK={}'.format(
                total, mfix_exe, nodesi, nodesj, nodesk)

        run_cmd = '{} -f {}'.format(mfix_exe, self.get_mfix_dat())
        self.run_thread.start_command(cmd=run_cmd, cwd=self.get_project_dir())

    def update_residuals(self):
        self.ui.residuals.setText(str(self.update_residuals_thread.residuals))
        if self.update_residuals_thread.job_done:
            self.ui.mfix_browser.setHTML('')

    def connect_mfix(self):
        """ connect to running instance of mfix """
        url = "http://{}:{}".format(self.ui.run.mfix_host.text(),
                                    self.ui.run.mfix_port.text())
        log = logging.getLogger(__name__)
        log.debug("trying to connect to {}".format(url))
        qurl = QUrl(url)
        self.ui.interact.mfix_browser.load(qurl)

        self.update_residuals_thread = UpdateResidualsThread()
        self.update_residuals_thread.sig.connect(self.update_residuals)
        self.update_residuals_thread.start()

        self.change_pane('interact')


    def get_mfix_dat(self):
        name = self.project.run_name.value
        for char in ('.', '"', "'", '/', '\\', ':'):
            name = name.replace(char, '_')
        mfix_dat = name + '.mfx'
        return os.path.join(self.get_project_dir(), mfix_dat)

    # --- open/save/new ---
    def save_project(self):
        project_dir = self.settings.value('project_dir')

        # export geometry
        self.vtkwidget.export_stl(os.path.join(project_dir, 'geometry.stl'))

        self.setWindowTitle('MFIX - %s' % project_dir)
        self.project.writeDatFile(self.get_mfix_dat())

    def unsaved(self):
        project_dir = self.settings.value('project_dir')
        self.setWindowTitle('MFIX - %s *' % project_dir)

    def new_project(self, project_dir=None):
        if not project_dir:
            project_dir = str(
                QtWidgets.QFileDialog.getExistingDirectory(
                    self, 'Create Project in Directory',
                    "",
                    QtWidgets.QFileDialog.ShowDirsOnly))
        if len(project_dir) < 1:
            return
        try:
            shutil.copyfile('mfix.dat.template',
                            os.path.join(project_dir, 'mfix.dat'))
        except IOError:
            self.message(title='Warning',
                         icon='warning',
                         text=('You do not have write access to this '
                               'directory.\n'),
                         buttons=['ok'],
                         default='ok',
                         )
            return
        self.open_project(project_dir)

    def open_project(self, project_dir=None):
        """
        Open MFiX Project
        """
        if not project_dir:
            project_dir = str(
                 QtWidgets.QFileDialog.getExistingDirectory(
                     self, 'Open Project Directory',
                     "",
                     QtWidgets.QFileDialog.ShowDirsOnly))

        if len(project_dir) < 1:
            return

        writable = True
        try:
            import tempfile
            testfile = tempfile.TemporaryFile(dir=project_dir)
            testfile.close()
        except IOError:
            writable = False

        if not writable:
            self.message(title='Warning',
                         icon='warning',
                         text=('You do not have write access to this '
                               'directory.\n'),
                         buttons=['ok'],
                         default='ok',
                         )
            return

        mfix_dat = os.path.abspath(os.path.join(project_dir, 'mfix.dat'))

        if not os.path.exists(mfix_dat):
            self.message(title='Warning',
                         icon='warning',
                         text=('mfix.dat file does not exist in this '
                               'directory.\n'),
                         buttons=['ok'],
                         default='ok',
                         )
            return

        self.settings.setValue('project_dir', project_dir)
        self.setWindowTitle('MFIX - %s' % project_dir)

        # read the file
        with open(mfix_dat, 'r') as mfix_dat_file:
            src = mfix_dat_file.read()
        self.ui.mfix_dat_source.setPlainText(src)
        # self.mode_changed('developer')

        self.project.load_mfix_dat(mfix_dat)

        self.ui.model_setup.energy_eq.setChecked(self.project['energy_eq'])

        # cgw - lots more model setup todo here

    # --- region methods ---
    def new_region(self):
        'create a new region'

        data = self.ui.regions.tablewidget_regions.value

        name = get_unique_string('new', list(data.keys()))

        data[name] = {'shape': 'box', 'from': [0, 0, 0], 'to': [0, 0, 0]}

        self.vtkwidget.new_region(name, data[name])

        self.ui.regions.tablewidget_regions.set_value(data)
        self.ui.regions.tablewidget_regions.fit_to_contents()

    def delete_region(self):
        'remove the currently selected region'

        row = self.ui.regions.tablewidget_regions.current_row()

        if row >= 0:
            data = self.ui.regions.tablewidget_regions.value
            name = list(data.keys())[row]
            data.pop(name)
            self.ui.regions.tablewidget_regions.set_value(data)
            self.vtkwidget.delete_region(name)

            if data:
                self.ui.regions.groupbox_region_parameters.setEnabled(True)
            else:
                self.ui.regions.groupbox_region_parameters.setEnabled(False)

    def copy_region(self):
        'copy the currently selected region'

        row = self.ui.regions.tablewidget_regions.current_row()

        if row >= 0:
            data = self.ui.regions.tablewidget_regions.value
            name = list(data.keys())[row]
            new_region = copy.deepcopy(data[name])

            new_name = get_unique_string(name, list(data.keys()))
            data[new_name] = new_region

            self.ui.regions.tablewidget_regions.set_value(data)

            self.vtkwidget.new_region(new_name, data[new_name])

    def update_region_parameters(self):
        'a new region was selected, update region widgets'

        row = self.ui.regions.tablewidget_regions.current_row()

        if row >= 0:
            self.ui.regions.groupbox_region_parameters.setEnabled(True)

            data = self.ui.regions.tablewidget_regions.value
            name = list(data.keys())[row]

            self.ui.regions.lineedit_regions_name.updateValue(None, name)
            self.ui.regions.combobox_regions_shape.updateValue(
                None,
                data[name]['shape'])

            for widget, value in zip([self.ui.regions.lineedit_regions_from_x,
                                      self.ui.regions.lineedit_regions_from_y,
                                      self.ui.regions.lineedit_regions_from_z],
                                     data[name]['from']
                                     ):
                widget.updateValue(None, value)

            for widget, value in zip([self.ui.regions.lineedit_regions_to_x,
                                      self.ui.regions.lineedit_regions_to_y,
                                      self.ui.regions.lineedit_regions_to_z],
                                     data[name]['to']
                                     ):
                widget.updateValue(None, value)

            self.ui.regions.tablewidget_regions.fit_to_contents()

        else:
            self.ui.regions.groupbox_region_parameters.setEnabled(False)

    def region_value_changed(self, widget, value, args):

        row = self.ui.regions.tablewidget_regions.current_row()
        data = self.ui.regions.tablewidget_regions.value
        name = list(data.keys())[row]
        key = value.keys()[0]

        if 'to' in key or 'from' in key:
            item = key.split('_')
            index = ['x', 'y', 'z'].index(item[1])
            data[name][item[0]][index] = value.values()[0]

        elif 'name' in key:

            new_name = get_unique_string(value.values()[0], list(data.keys()))
            data = OrderedDict([(new_name, v) if k == name else (k, v) for
                                k, v in data.items()])

            self.vtkwidget.change_region_name(name, new_name)

        elif 'shape' in key:
            data[name]['shape'] = value.values()[0]

            self.vtkwidget.change_region_shape(name, data[name])

        self.vtkwidget.update_region(name, data[name])

        self.ui.regions.tablewidget_regions.set_value(data)


# --- Threads ---
class MfixThread(QThread):
    line_printed = pyqtSignal(str, str)

    def __init__(self, parent):
        super(MfixThread, self).__init__(parent)
        self.cmd = None
        self.cwd = None

    def start_command(self, cmd, cwd):
        log = logging.getLogger(__name__)
        log.debug("Running in %s : %s", cwd, cmd)
        self.cmd = cmd
        self.cwd = cwd
        self.start()

    def run(self):
        if self.cmd:
            popen = subprocess.Popen(self.cmd,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=True, cwd=self.cwd)
            lines_iterator = iter(popen.stdout.readline, b"")
            for line in lines_iterator:
                self.line_printed.emit(str(line), None)
            lines_iterator = iter(popen.stderr.readline, b"")
            for line in lines_iterator:
                self.line_prvinted.emit(str(line), "red")


class MonitorExecutablesThread(QThread):

    sig = QtCore.pyqtSignal()

    def __init__(self, parent):
        QThread.__init__(self)
        self.parent = parent
        self.mfix_home = get_mfix_home()
        self.output = self.get_output()

    def get_output(self):
        """ returns a dict mapping full (mfix,pymfix) paths
        to configuration options."""

        config_options = {}
        # TODO: get config options from mfix itself, via 'mfix --print-config',
        # rather than inferring from build dir names

        # Check system PATH dirs first
        PATH = os.environ.get("PATH")
        if PATH:
            for dir in PATH.split(os.pathsep):
                for name in 'mfix', 'mfix.exe', 'pymfix', 'pymfix.exe':
                    path = os.path.join(dir, name)
                    if os.path.isfile(path):
                        config_options[path] = '' # Unknown configuration options


        # Look in subdirs of build dir
        build = os.path.join(self.mfix_home,'build')
        all_builddirs = os.path.join(build,'*')
        for exe in 'mfix', 'mfix.exe', 'pymfix', 'pymfix.exe':
            for path in glob.glob(os.path.join(all_builddirs, exe)):
                config_options[path] = '' # Unknown configuration options

        # Look for named symlinks in run_dir/.build
        run_dir = self.parent.get_project_dir()
        if run_dir is not None:
            dot_build = os.path.join(run_dir,'.build')
            if os.path.isdir(dot_build):
                for subdir in os.listdir(dot_build):
                    path = os.path.join(dot_build, subdir)
                    if os.path.islink(path): # Symlink name indicates config options
                        for name in 'mfix', 'mfix.exe', 'pymfix', 'pymfix.exe':
                            exe = os.path.join(run_dir, name)
                            if os.path.isfile(exe):
                                config_options[exe] = subdir

        return config_options

    def run(self):
        while True:
            output = self.get_output()
            if output != self.output:
                self.output = output
                self.sig.emit()
            time.sleep(1)

class MonitorOutputFilesThread(QThread):

    sig = QtCore.pyqtSignal(object)

    def run(self):
        while True:
            self.job_done = False
            output = glob.glob(run_name+'*')
            if output:
                self.sig.emit('output')
            else:
                self.sig.emit('no_output')
            time.sleep(1)


class UpdateResidualsThread(QThread):

    sig = QtCore.pyqtSignal(object)

    def run(self):
        while True:
            self.job_done = False
            try:
                self.residuals = urlopen('http://localhost:5000/residuals').read()
            except Exception:
                log = logging.getLogger(__name__)
                log.debug("cannot retrieve residuals; pymfix process must have terminated.")
                self.job_done = True
                return
            time.sleep(1)
            self.sig.emit('update')


# --- Project Manager ---
class ProjectManager(Project):
    '''
    Manages the editing of the MFiX.dat file
    '''
    def __init__(self, parent=None, keyword_doc=None):
        Project.__init__(self, keyword_doc=keyword_doc)

        self.parent = parent
        self.widgetList = []
        self.registered_keywords = set()

    def submit_change(self, widget, newValueDict, args=None,
                      forceUpdate=False):
        '''
        Submit a value change

        Examples:
        submitChange(lineEdit, {'run_name':'new run name'}, args)
        '''

        if isinstance(args, int):
            args = [args]

        for key, newValue in newValueDict.items():

            if isinstance(newValue, dict):
                if args:
                    for ind, value in newValue.items():
                        self._change(widget, key, value, args=args+[ind],
                                     forceUpdate=forceUpdate)
                else:
                    for ind, value in newValue.items():
                        self._change(widget, key, value, args=[ind],
                                     forceUpdate=forceUpdate)
            else:
                self._change(widget, key, newValue, args=args,
                             forceUpdate=forceUpdate)

        self._cleanDeletedItems()

    def _change(self, widget, key, newValue, args=None, forceUpdate=False):

        key = key.lower()
        updatedValue = None

        if isinstance(newValue, Keyword):
            keyword = newValue
            newValue = keyword.value
        else:
            keyword = None

        if args is None:
            args = []

        updatedValue = self.addKeyword(key, newValue, args)

        if args:
            keystring = '{}({}) = {}'.format(
                key,
                ','.join([str(arg) for arg in args]),
                str(updatedValue))
        else:
            keystring = '{} = {}'.format(key, str(updatedValue))
        self.parent.print_internal(keystring)

        if updatedValue is not None:
            for wid, keys in self.widgetList:
                if not wid == widget:
                    if key in keys or 'all' in keys:
                        wid.updateValue(key, updatedValue, args)

    def load_mfix_dat(self, mfixDat):

        self.parsemfixdat(fname=mfixDat)

        # emit loaded keys
        for keyword in self.keywordItems():
            self.submit_change('Loading', {keyword.key: keyword.value},
                               args=keyword.args, forceUpdate=True)

    def register_widget(self, widget, keys=None):
        '''
        Register a widget with the project manager. The widget must have a
        valueUpdated signal to connect to.

        Example:
        registerWidget(widget, ['list', 'of', 'keys'])
        '''

        LOG.debug(
            'ProjectManager: Registering {} with keys {}'.format(
                widget.objectName(),
                keys)
            )
        self.widgetList.append([widget, keys])
        widget.value_updated.connect(self.submit_change)

        self.registered_keywords = self.registered_keywords.union(set(keys))

    def objectName(self):
        return 'Project Manager'


if __name__ == '__main__':
    qapp = QtWidgets.QApplication(sys.argv)

    mfix = MfixGui(qapp)

    mfix.show()

    # have to initialize vtk after the widget is visible!
    mfix.vtkwidget.vtkiren.Initialize()

    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    qapp.exec_()

    qapp.deleteLater()
    sys.exit()
