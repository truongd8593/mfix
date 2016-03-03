""" MFIX GUI """


# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals
import logging
import os
import signal
import sys
import subprocess


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
from tools.mfixproject import Project
from tools.general import (get_image_path, make_callback, get_icon,
                           widget_iter, set_script_directory)

logging.basicConfig(stream=sys.stdout,
                    filemode='w', level=logging.DEBUG,
                    format='%(name)s - %(levelname)s - %(message)s')

SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), ))
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))
LOG = logging.getLogger(__name__)
LOG.debug(SCRIPT_DIRECTORY)
set_script_directory(SCRIPT_DIRECTORY)  # should this be in an __init__.py


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
        self.ui = uic.loadUi(os.path.join('uifiles', 'gui.ui'), self)

        # load settings
        self.settings = QSettings('MFIX', 'MFIX')

        # set title and icon
        self.setWindowTitle('MFIX')
        self.setWindowIcon(get_icon('mfix.png'))

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

        # --- icons ---
        self.ui.toolbutton_new.setIcon(get_icon('newfolder.png'))
        self.ui.toolbutton_open.setIcon(get_icon('openfolder.png'))
        self.ui.toolbutton_save.setIcon(get_icon('save.png'))

        self.ui.toolbutton_add_geometry.setIcon(get_icon('geometry.png'))
        self.ui.toolbutton_add_filter.setIcon(get_icon('filter.png'))
        self.ui.toolbutton_remove_geometry.setIcon(get_icon('remove.png'))
        self.ui.toolbutton_copy_geometry.setIcon(get_icon('copy.png'))
        self.ui.toolbutton_geometry_union.setIcon(get_icon('union.png'))
        self.ui.toolbutton_geometry_intersect.setIcon(
            get_icon('intersect.png'))
        self.ui.toolbutton_geometry_difference.setIcon(
            get_icon('difference.png'))

        self.ui.toolButtonRegionAdd.setIcon(get_icon('add.png'))
        self.ui.toolButtonRegionDelete.setIcon(get_icon('remove.png'))
        self.ui.toolButtonRegionCopy.setIcon(get_icon('copy.png'))

        self.ui.toolButtonFluidSpeciesAdd.setIcon(get_icon('add.png'))
        self.ui.toolButtonFluidSpeciesDelete.setIcon(get_icon('remove.png'))
        self.ui.toolButtonFluidSpeciesCopy.setIcon(get_icon('copy.png'))

        self.ui.toolButtonTFMSolidsAdd.setIcon(get_icon('add.png'))
        self.ui.toolButtonTFMSolidsDelete.setIcon(get_icon('remove.png'))
        self.ui.toolButtonTFMSolidsCopy.setIcon(get_icon('copy.png'))
        self.ui.toolButtonTFMSolidsDatabase.setIcon(get_icon('download.png'))

        self.ui.toolButtonTFMSolidsSpeciesAdd.setIcon(get_icon('add.png'))
        self.ui.toolButtonTFMSolidsSpeciesDelete.setIcon(
            get_icon('remove.png'))
        self.ui.toolButtonTFMSolidsSpeciesCopy.setIcon(get_icon('copy.png'))

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

        self.ui.run_name.textChanged.connect(self.unsaved)
        self.ui.description.textChanged.connect(self.unsaved)

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
            make_handler(self.ui.build_output))
        self.run_thread.line_printed.connect(make_handler(self.ui.run_output))
        self.clear_thread.line_printed.connect(
            make_handler(self.ui.run_output))

        # --- vtk setup ---
        self.__setup_vtk_widget()

        # --- workflow setup ---
        if NodeWidget is not None:
            self.__setup_workflow_widget()

        # --- default ---
        self.mode_changed('modeler')
        top = self.ui.treewidget_model_navigation.topLevelItem(0)
        self.ui.treewidget_model_navigation.setCurrentItem(top)

        # autoload last project
        if self.get_project_dir():
            self.open_project(self.get_project_dir())

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

        # connect boolean
        for key, btn in self.booleanbtndict.items():
            btn.pressed.connect(
                make_callback(self.vtkwidget.boolean_operation, key))

        # connect primitive lineedits
        for widget in widget_iter(self.ui.stackedWidgetGeometryDetails):
            if isinstance(widget, QtGui.QLineEdit):
                widget.editingFinished.connect(
                    make_callback(self.vtkwidget.primitive_edited, widget))

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

        self.ui.stackedwidget_mode.setCurrentIndex(current_index)

        for key, btn in self.modebuttondict.items():
            btn.setChecked(mode == key)

    def navigation_changed(self):
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

            self.ui.stackedWidgetTaskPane.setCurrentIndex(current_index)

    def build_mfix(self):
        """ build mfix """
        mfix_home = os.path.dirname(
            os.path.dirname(os.path.realpath(__file__)))
        dmp = '--dmp' if self.ui.dmp_button.isChecked() else ''
        smp = '--smp' if self.ui.smp_button.isChecked() else ''
        build_cmd = os.path.join(
            mfix_home, 'configure_mfix %s %s && make pymfix' % (smp, dmp))
        self.build_thread.start_command(build_cmd, self.get_project_dir())

    def clear_output(self):
        """ build mfix """
        self.run_thread.start_command(
            'echo "Removing:";'
            ' ls *.LOG *.OUT *.RES *.SP? *.pvd *vtp;'
            ' rm -f *.LOG *.OUT *.RES *.SP? *.pvd *vtp',
            self.get_project_dir())

    def run_mfix(self):
        """ build mfix """
        pymfix_exe = os.path.join(self.get_project_dir(), 'pymfix')
        self.run_thread.start_command(pymfix_exe, self.get_project_dir())

    def connect_mfix(self):
        """ connect to running instance of mfix """
        self.ui.mfix_browser.load(QUrl("http://localhost:5000"))

    def save_project(self):
        project_dir = self.settings.value('project_dir')

        self.project._keywordDict['run_name'] = self.ui.run_name.text()
        self.project._keywordDict['description'] = self.ui.description.text()

        self.setWindowTitle('MFIX - %s' % project_dir)
        self.project.save(os.path.join(project_dir, 'mfix.dat'))

    def unsaved(self):
        project_dir = self.settings.value('project_dir')
        self.setWindowTitle('MFIX - %s *' % project_dir)

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
                               'directory.\n'
                               'Please fix and try again.'),
                         buttons=['ok'],
                         default='ok',
                         )
            return

        self.settings.setValue('project_dir', project_dir)

        # self.recentProjects.append(self.projectDir)
        # self.updateRecentProjectMenu()

        self.setWindowTitle('MFIX - %s' % project_dir)

        # Look for mfix.dat file
        mfix_dat = os.path.abspath(os.path.join(project_dir, 'mfix.dat'))
        # mfix_dat = os.path.abspath(os.path.join(project_dir, 'mfix.dat'))

        if os.path.exists(mfix_dat):
            # mylogger.debug('found mfix.dat: {}'.format(mfix_dat))
            # check to see if file is already open
            # if not self.codeEditor.is_file_opened(mfix_dat):
                # self.codeEditor.load(mfix_dat)
            # parse and update widget values
            # self.setWidgetInit(False)
            # self.projectManager.loadMfixDat(mfix_dat)
            # self.setWidgetInit(True)

            src = open(mfix_dat).read()
            self.ui.mfix_dat_source.setPlainText(src)
            self.mode_changed('developer')

            self.project = Project(mfix_dat)
            self.ui.run_name.setText(str(self.project['run_name']))
            self.ui.description.setText(str(self.project['description']))

        else:
            print("mfix.dat doesn't exist")
            # self.newMfixDat()


# --- Threads ---
class MfixThread(QThread):

    line_printed = pyqtSignal(str)

    def __init__(self, parent):
        super(MfixThread, self).__init__(parent)
        self.cmd = None
        self.cwd = None

    def start_command(self, cmd, cwd):
        self.cmd = cmd
        self.cwd = cwd
        self.start()

    def run(self):
        if self.cmd:
            popen = subprocess.Popen(
                self.cmd, stdout=subprocess.PIPE, shell=True, cwd=self.cwd)
            lines_iterator = iter(popen.stdout.readline, b"")
            for line in lines_iterator:
                self.line_printed.emit(line)


class RunThread(MfixThread):
    line_printed = pyqtSignal(str)


class BuildThread(MfixThread):
    line_printed = pyqtSignal(str)


class ClearThread(MfixThread):
    line_printed = pyqtSignal(str)

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
