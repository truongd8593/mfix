""" MFIX GUI """


# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals
import logging
import os
import shutil
import signal
import sys
import subprocess
import time
import urllib2
import glob


# import qt
from qtpy import QtCore, QtGui
from qtpy.QtCore import QObject, QThread, pyqtSignal, QUrl, QTimer, QSettings

# TODO: add pyside?
try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic


# local imports
try:
    from pyqtnode import NodeWidget
except ImportError:
    NodeWidget = None
from widgets.vtkwidget import VtkWidget
from widgets.base import (LineEdit, CheckBox, ComboBox, SpinBox, DoubleSpinBox)
from tools.mfixproject import Project, KeyWord
from tools.general import (get_image_path, make_callback, get_icon,
                           widget_iter, set_script_directory)
from tools.namelistparser import buildKeywordDoc

logging.basicConfig(stream=sys.stdout,
                    filemode='w', level=logging.DEBUG,
                    format='%(name)s - %(levelname)s - %(message)s')

SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), ))
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))
LOG = logging.getLogger(__name__)
LOG.debug(SCRIPT_DIRECTORY)
set_script_directory(SCRIPT_DIRECTORY)  # should this be in an __init__.py?


# --- Main Gui ---
class MfixGui(QtGui.QMainWindow):
    '''
    Main window class handling all gui interactions
    '''
    def __init__(self, app, parent=None):
        QtGui.QMainWindow.__init__(self, parent)

        # reference to qapp instance
        self.app = app

        # load ui file
        self.customWidgets={'LineEdit':      LineEdit,
                            'CheckBox':      CheckBox,
                            'ComboBox':      ComboBox,
                            'DoubleSpinBox': DoubleSpinBox,
                            'SpinBox':       SpinBox,
                            }
        self.ui = uic.loadUi(os.path.join('uifiles', 'gui.ui'), self)

        # load settings
        self.settings = QSettings('MFIX', 'MFIX')

        # set title and icon
        self.setWindowTitle('MFIX')
        self.setWindowIcon(get_icon('mfix.png'))

        self.project = ProjectManager(self)

        # --- data ---
        self.modebuttondict = {'modeler':   self.ui.pushButtonModeler,
                               'workflow':  self.ui.pushButtonWorkflow,
                               'developer': self.ui.pushButtonDeveloper,
                               }

        self.booleanbtndict = {
            'union':      self.ui.toolbutton_geometry_union,
            'intersection':  self.ui.toolbutton_geometry_intersect,
            'difference': self.ui.toolbutton_geometry_difference,
            }
        self.animation_speed = 400
        self.animating = False
        self.stack_animation = None

        # --- icons ---
        # loop through all widgets, because I am lazy
        for widget in widget_iter(self.ui):
            if isinstance(widget, QtGui.QToolButton):
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

        self.ui.toolbutton_add_geometry.setIcon(get_icon('geometry.png'))
        self.ui.toolbutton_add_filter.setIcon(get_icon('filter.png'))
        self.ui.toolbutton_geometry_union.setIcon(get_icon('union.png'))
        self.ui.toolbutton_geometry_intersect.setIcon(
            get_icon('intersect.png'))
        self.ui.toolbutton_geometry_difference.setIcon(
            get_icon('difference.png'))

        self.ui.toolButtonTFMSolidsDatabase.setIcon(get_icon('download.png'))

        # --- Connect Signals to Slots---
        # open/save/new project
        self.ui.toolbutton_open.pressed.connect(self.open_project)
        self.ui.toolbutton_save.pressed.connect(self.save_project)

        # mode (modeler, workflow, developer)
        for mode, btn in self.modebuttondict.items():
            btn.released.connect(make_callback(self.mode_changed, mode))

        # navigation tree
        self.ui.treewidget_model_navigation.itemSelectionChanged.connect(
            self.navigation_changed)

        # build/run/connect MFIX
        self.ui.build_mfix_button.pressed.connect(self.build_mfix)
        self.ui.run_mfix_button.pressed.connect(self.run_mfix)
        self.ui.connect_mfix_button.pressed.connect(self.connect_mfix)
        self.ui.clear_output_button.pressed.connect(self.clear_output)

        # --- Threads ---
        self.build_thread = BuildThread(self)
        self.run_thread = RunThread(self)
        self.clear_thread = ClearThread(self)

        def make_handler(qtextbrowser):
            " make a closure to read stdout from external process "

            def handle_line(line):
                " closure to read stdout from external process "

                log = logging.getLogger(__name__)
                log.debug(str(line).strip())
                cursor = qtextbrowser.textCursor()
                cursor.movePosition(cursor.End)
                cursor.insertText(line)
                qtextbrowser.ensureCursorVisible()

            return handle_line

        self.build_thread.line_printed.connect(
            make_handler(self.ui.command_output))
        self.run_thread.line_printed.connect(
            make_handler(self.ui.command_output))
        self.clear_thread.line_printed.connect(
            make_handler(self.ui.command_output))

        # --- setup simple widgets ---
        self.__setup_simple_keyword_widgets()

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

    def __setup_simple_keyword_widgets(self):
        """
        Look for and connect simple keyword widgets to the project manager.
        Keyword informtation from the namelist doc strings is added to each
        keyword widget. The widget must be named: *_keyword_<keyword> where
        <keyword> is the actual keyword. For example:
        lineedit_keyword_run_name
        """

        # build keyword documentation from namelist docstrings
        self.keyword_doc = buildKeywordDoc(os.path.join(SCRIPT_DIRECTORY,
                                                        os.pardir, 'model'))

        # loop through all widgets looking for *_keyword_<keyword>
        for widget in widget_iter(self.ui):
            name_list = str(widget.objectName()).lower().split('_')

            if 'keyword' in name_list:
                keyword = '_'.join(name_list[name_list.index('keyword')+1:])

                # set the key attribute to the keyword
                widget.key = keyword

                # add info from keyword documentation
                if keyword in self.keyword_doc:
                    widget.setdtype(self.keyword_doc[keyword]['dtype'])

                    if 'required' in self.keyword_doc[keyword]:
                        widget.setValInfo(
                            req=self.keyword_doc[keyword]['required'] == 'true')
                    if 'validrange' in self.keyword_doc[keyword]:
                        if 'max' in self.keyword_doc[keyword]['validrange']:
                            widget.setValInfo(
                                _max = self.keyword_doc[keyword]['validrange']['max'])
                        if 'min' in self.keyword_doc[keyword]['validrange']:
                            widget.setValInfo(
                                _min = self.keyword_doc[keyword]['validrange']['min'])

                    if 'initpython' in self.keyword_doc[keyword]:
                        widget.default(
                            self.keyword_doc[keyword]['initpython'])

                    if isinstance(widget, QtGui.QComboBox) and widget.count() < 1:
                            widget.addItems(list(self.mfixKeyWordDoc[key]['valids'].keys()))

                # register the widget with the project manager
                self.project.register_widget(widget, [keyword])
                
                # connect to unsaved method
                widget.value_updated.connect(self.unsaved)

    def __setup_vtk_widget(self):
        " setup the vtk widget "

        self.vtkwidget = VtkWidget(parent=self)
        self.ui.horizontalLayoutModelGraphics.addWidget(self.vtkwidget)

        # --- geometry button ---
        self.add_geometry_menu = QtGui.QMenu(self)
        self.ui.toolbutton_add_geometry.setMenu(self.add_geometry_menu)

        action = QtGui.QAction('STL File',  self.add_geometry_menu)
        action.triggered.connect(self.vtkwidget.add_stl)
        self.add_geometry_menu.addAction(action)

        self.add_geometry_menu.addSeparator()

        for geo in self.vtkwidget.primitivedict.keys():
            action = QtGui.QAction(geo, self.add_geometry_menu)
            action.triggered.connect(
                make_callback(self.vtkwidget.add_primitive, geo))
            self.add_geometry_menu.addAction(action)

        self.add_geometry_menu.addSeparator()

        for geo in self.vtkwidget.parametricdict.keys():
            action = QtGui.QAction(geo.replace('_', ' '),
                                   self.add_geometry_menu)
            action.triggered.connect(
                make_callback(self.vtkwidget.add_parametric, geo))
            self.add_geometry_menu.addAction(action)

        # --- filter button ---
        self.add_filter_menu = QtGui.QMenu(self)
        self.ui.toolbutton_add_filter.setMenu(self.add_filter_menu)

        for geo in self.vtkwidget.filterdict.keys():
            action = QtGui.QAction(geo.replace('_', ' '),
                                   self.add_filter_menu)
            action.triggered.connect(
                make_callback(self.vtkwidget.add_filter, geo))
            self.add_filter_menu.addAction(action)

        # tree widget icons
        self.ui.treeWidgetGeometry.setStyleSheet(
            "QTreeView::indicator:unchecked {image: url(%s);}"
            "QTreeView::indicator:checked {image: url(%s);}"
            % (get_image_path('visibilityofftransparent.png'),
               get_image_path('visibility.png'))
            )

        # setup signals
        self.ui.toolbutton_remove_geometry.pressed.connect(
            self.vtkwidget.remove_geometry)
        self.ui.toolbutton_copy_geometry.pressed.connect(
            self.vtkwidget.copy_geometry)

        # connect boolean
        for key, btn in self.booleanbtndict.items():
            btn.pressed.connect(
                make_callback(self.vtkwidget.boolean_operation, key))

        # connect parameter widgets
        for widget in widget_iter(self.ui.stackedWidgetGeometryDetails):
            if isinstance(widget, QtGui.QLineEdit):
                widget.editingFinished.connect(
                    make_callback(self.vtkwidget.parameter_edited, widget))
            elif isinstance(widget, QtGui.QCheckBox):
                widget.stateChanged.connect(
                    make_callback(self.vtkwidget.parameter_edited, widget))

        # --- mesh ---
        # connect mesh tab btns
        for i, btn in enumerate([self.ui.pushbutton_mesh_uniform,
                                 self.ui.pushbutton_mesh_controlpoints,
                                 self.ui.pushbutton_mesh_mesher]):
            btn.pressed.connect(
                make_callback(self.vtkwidget.change_mesh_tab, i, btn))

        self.ui.pushbutton_mesh_autosize.pressed.connect(
            self.vtkwidget.auto_size_mesh_extents)

        # connect unifrom mesh
        for widget in [self.ui.lineedit_mesh_min_x,
                       self.ui.lineedit_mesh_max_x,
                       self.ui.lineedit_mesh_min_y,
                       self.ui.lineedit_mesh_max_y,
                       self.ui.lineedit_mesh_min_z,
                       self.ui.lineedit_mesh_max_z,
                       self.ui.lineedit_mesh_cells_x,
                       self.ui.lineedit_mesh_cells_y,
                       self.ui.lineedit_mesh_cells_z
                       ]:
            widget.editingFinished.connect(self.vtkwidget.update_mesh)

        # connect mesher
        self.ui.pushbutton_generate_mesh.pressed.connect(self.vtkwidget.mesher)

    def __setup_workflow_widget(self):

        self.nodeChart = NodeWidget(showtoolbar=False)
        # Build defualt node library
        self.nodeChart.nodeLibrary.buildDefualtLibrary()
        self.ui.horizontalLayoutPyqtnode.addWidget(self.nodeChart)

    def get_project_dir(self):
        " get the current project directory"

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
                    QtCore.Qt.MatchContains | QtCore.Qt.MatchRecursive,
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
        msgBox = QtGui.QMessageBox(self)
        msgBox.setWindowTitle(title)

        # Icon
        if icon == 'warning':
            icon = QtGui.QMessageBox.Warning
        else:
            icon = QtGui.QMessageBox.Information

        msgBox.setIcon(icon)

        # Text
        msgBox.setText(text)

        if infoText:
            msgBox.setInformativeText(infoText)

        if detailedtext:
            msgBox.setDetailedText(detailedtext)

        # buttons
        qbuttonDict = {'ok':      QtGui.QMessageBox.Ok,
                       'yes':     QtGui.QMessageBox.Yes,
                       'no':      QtGui.QMessageBox.No,
                       'cancel':  QtGui.QMessageBox.Cancel,
                       'discard': QtGui.QMessageBox.Discard,
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
        LOG.debug(str(text).strip())
        cursor = self.ui.command_output.textCursor()
        cursor.movePosition(cursor.End)
        cursor.insertText(text)
        cursor.movePosition(cursor.End)
        self.ui.command_output.ensureCursorVisible()

    # --- mfix methods ---
    def set_mpi(self):
        " set MPI checkbox if more than one node selected "
        nodesi = int(self.ui.nodes_i.text())
        nodesj = int(self.ui.nodes_j.text())
        nodesk = int(self.ui.nodes_k.text())
        if max(nodesi, nodesj, nodesk) > 1:
            self.ui.dmp_button.setChecked(True)

    def set_nodes(self):
        " set all nodes to one if MPI checkbox is unchecked "
        if self.ui.dmp_button.isChecked():
            self.ui.nodes_i.setValue(1)
            self.ui.nodes_j.setValue(1)
            self.ui.nodes_k.setValue(1)

    def make_build_cmd(self):
        """ build mfix """
        mfix_home = os.path.dirname(
            os.path.dirname(os.path.realpath(__file__)))
        dmp = '--dmp' if self.ui.dmp_button.isChecked() else ''
        smp = '--smp' if self.ui.smp_button.isChecked() else ''
        return os.path.join(mfix_home, 'configure_mfix --python %s %s && make -j pymfix' % (smp, dmp))

    def build_mfix(self):
        """ build mfix """
        self.build_thread.start_command(self.make_build_cmd(),
                                        self.get_project_dir())

    def clear_output(self):
        """ build mfix """
        self.run_thread.start_command('echo "Removing:";'
                                      ' ls *.LOG *.OUT *.RES *.SP? *.pvd *vtp;'
                                      ' rm -f *.LOG *.OUT *.RES *.SP? *.pvd *vtp',
                                      self.get_project_dir())

    def run_mfix(self):
        """ build mfix """
        if not self.ui.dmp_button.isChecked():
            pymfix_exe = os.path.join(self.get_project_dir(), 'pymfix')
        else:
            nodesi = int(self.ui.nodes_i.text())
            nodesj = int(self.ui.nodes_j.text())
            nodesk = int(self.ui.nodes_k.text())
            total = nodesi*nodesj*nodesk
            pymfix_exe = 'mpirun -np {} ./pymfix NODESI={} NODESJ={} NODESK={}'.format(total, nodesi, nodesj, nodesk)

        build_and_run_cmd = '{} && {}'.format(self.make_build_cmd(),
                                              pymfix_exe)
        self.run_thread.start_command(build_and_run_cmd,
                                      self.get_project_dir())

    def update_residuals(self):
        self.ui.residuals.setText(self.updater.residuals)
        if self.updater.job_done:
            self.ui.mfix_browser.setHTML('')


    def connect_mfix(self):
        """ connect to running instance of mfix """
        url = "http://{}:{}".format(self.ui.mfix_host.text(),
                                    self.ui.mfix_port.text())
        log = logging.getLogger(__name__)
        log.debug("trying to connect to {}".format(url))
        qurl = QUrl(url)
        self.ui.mfix_browser.load(qurl)

        self.updater = UpdateResidualsThread()
        self.updater.sig.connect(self.update_residuals)
        self.updater.start()

        self.change_pane('interact')

    # --- open/save/new ---
    def save_project(self):
        project_dir = self.settings.value('project_dir')

        # export geometry
        self.vtkwidget.export_stl(os.path.join(project_dir, 'geometry.stl'))

        # save project
        self.setWindowTitle('MFIX - %s' % project_dir)
        self.project.writeDatFile(os.path.join(project_dir, 'mfix.dat'))

    def unsaved(self):
        project_dir = self.settings.value('project_dir')
        self.setWindowTitle('MFIX - %s *' % project_dir)

    def new_project(self, project_dir=None):
        if not project_dir:
            project_dir = str(
                QtGui.QFileDialog.getExistingDirectory(
                    self, 'Create Project in Directory',
                    "",
                    QtGui.QFileDialog.ShowDirsOnly))
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
                 QtGui.QFileDialog.getExistingDirectory(
                     self, 'Open Project Directory',
                     "",
                     QtGui.QFileDialog.ShowDirsOnly))

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

        self.ui.energy_eq.setChecked(self.project['energy_eq'])


# --- Threads ---
class MfixThread(QThread):

    line_printed = pyqtSignal(str)

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
                                     shell=True, cwd=self.cwd)
            lines_iterator = iter(popen.stdout.readline, b"")
            for line in lines_iterator:
                self.line_printed.emit(line)


class RunThread(MfixThread):
    line_printed = pyqtSignal(str)


class BuildThread(MfixThread):
    line_printed = pyqtSignal(str)


class ClearThread(MfixThread):
    line_printed = pyqtSignal(str)


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
                self.residuals = urllib2.urlopen('http://localhost:5000/residuals').read()
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
    def __init__(self, parent = None):
        Project.__init__(self)

        self.parent = parent
        self.widgetList = []

    def submit_change(self, widget, newValueDict, args=None, forceUpdate=False):
        '''
        Submit a value change

        Examples:
        submitChange(lineEdit, {'run_name':'new run name'}, args)
        '''

        if isinstance(args, int):
            args = [args]

        for key, newValue in newValueDict.items():
            if hasattr(widget, 'objectName'):
                name = widget.objectName()
            else:
                name = widget

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

        if isinstance(newValue, KeyWord):
            keyword = newValue
            newValue = keyword.value
        else:
            keyword = None

        if args is None:
            args = []

        updatedValue = self.addKeyword(key, newValue, args)

        if args:
            keystring = '{}({}) = {}'.format(key,
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

    def objectName(self):
        return 'Project Manager'


if __name__ == '__main__':
    qapp = QtGui.QApplication(sys.argv)

    mfix = MfixGui(qapp)

    mfix.show()

    # have to initialize vtk after the widget is visible!
    mfix.vtkwidget.vtkiren.Initialize()

    # exit with Ctrl-C at the terminal
    timer = QTimer()
    timer.start(500)  # You may change this if you wish.
    timer.timeout.connect(lambda: None)  # Let the interpreter run each 500 ms.
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    qapp.exec_()

    qapp.deleteLater()
    sys.exit()
