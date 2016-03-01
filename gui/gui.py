""" MFIX GUI """


# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals
import logging
import os
import signal
import sys
import subprocess
import re

from qtpy import QtCore, QtGui
from qtpy.QtCore import QObject, QThread, pyqtSignal, QUrl, QTimer, QSettings

# TODO: add pyside?
try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic

logging.basicConfig(stream=sys.stdout,
                    filemode='w', level=logging.DEBUG,
                    format='%(name)s - %(levelname)s - %(message)s')

SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), ))
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))
LOG = logging.getLogger(__name__)
LOG.debug(SCRIPT_DIRECTORY)

# local imports
# from pyqtnode import NodeWidget
from widgets.vtkwidget import VtkWidget

def get_image_path(name):
    " get path to images "
    path = os.path.join(SCRIPT_DIRECTORY, 'icons', name)

    if os.name == 'nt':
        path = path.replace('\\', '//')

    return path

def make_callback(func, param):
    '''
    Helper function to make sure lambda functions are cached and not lost.
    '''
    return lambda: func(param)

def get_icon(name, default=None, resample=False):
    """Return image inside a QIcon object
    default: default image name or icon
    resample: if True, manually resample icon pixmaps for usual sizes
    (16, 24, 32, 48, 96, 128, 256). This is recommended for QMainWindow icons
    created from SVG images on non-Windows platforms due to a Qt bug (see
    http://code.google.com/p/spyderlib/issues/detail?id=1314)."""
    if default is None:
        icon = QtGui.QIcon(get_image_path(name))
    elif isinstance(default, QtGui.QIcon):
        icon_path = get_image_path(name)
        icon = default if icon_path is None else QtGui.QIcon(icon_path)
    else:
        icon = QtGui.QIcon(get_image_path(name, default))
    if resample:
        icon0 = QtGui.QIcon()
        for size in (16, 24, 32, 48, 96, 128, 256, 512):
            icon0.addPixmap(icon.pixmap(size, size))
        return icon0
    else:
        return icon

def get_unique_string(base, listofstrings):
    " uniquify a string "
    if base in listofstrings:
        # look for number at end
        nums = re.findall('[\d]+', base)
        if nums:
            number = int(nums[-1]) + 1
            base = base.replace(nums[-1], '')
        else:
            number = 1

        base = get_unique_string(''.join([base, str(number)]), listofstrings)

    return base

def widget_iter(widget):
    for child in widget.children():
        if child.children():
            for child2 in widget_iter(child):
                yield child2
        yield child

# --- Main Gui ---
class MfixGui(QtGui.QMainWindow):
    '''
    Main window class handling all gui interactions
    '''
    def __init__(self, app, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.app = app

        self.ui = uic.loadUi('uifiles/gui.ui', self)

        self.settings = QSettings('MFIX', 'MFIX')

        self.setWindowTitle('MFIX')
        self.setWindowIcon(get_icon('mfix.png'))

        # --- data ---
        self.modebuttondict = {'modeler':self.ui.pushButtonModeler,
                               'workflow':self.ui.pushButtonWorkflow,
                               'developer':  self.ui.pushButtonDeveloper,
                              }

        self.booleanbtndict = {
            'union':      self.ui.toolButtonGeometryUnion,
            'intersection':  self.ui.toolButtonGeometryIntersect,
            'difference': self.ui.toolButtonGeometryDifference,
            }

        # --- icons ---
        self.ui.toolButtonNew.setIcon(get_icon('newfolder.png'))
        self.ui.toolButtonOpen.setIcon(get_icon('openfolder.png'))
        self.ui.toolButtonSave.setIcon(get_icon('save.png'))
        
        self.ui.toolButtonAddGeometry.setIcon(get_icon('add.png'))
        self.ui.toolButtonRemoveGeometry.setIcon(get_icon('remove.png'))
        self.ui.toolButtonGeometryUnion.setIcon(get_icon('union.png'))
        self.ui.toolButtonGeometryIntersect.setIcon(get_icon('intersect.png'))
        self.ui.toolButtonGeometryDifference.setIcon(get_icon('difference.png'))

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
        self.ui.toolButtonTFMSolidsSpeciesDelete.setIcon(get_icon('remove.png'))
        self.ui.toolButtonTFMSolidsSpeciesCopy.setIcon(get_icon('copy.png'))
        
        self.ui.toolButtonOpen.pressed.connect(self.open_project)


        self.ui.build_mfix_button.pressed.connect(self.build_mfix)
        self.ui.run_mfix_button.pressed.connect(self.run_mfix)
        self.ui.connect_mfix_button.pressed.connect(self.connect_mfix)
        self.ui.clear_output_button.pressed.connect(self.clear_output)

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

        self.build_thread.line_printed.connect(make_handler(self.ui.build_output))
        self.run_thread.line_printed.connect(make_handler(self.ui.run_output))
        self.clear_thread.line_printed.connect(make_handler(self.ui.run_output))

        
        # --- vtk setup ---
        self.__setup_vtk_widget()

        # --- workflow setup ---
        # self.nodeChart = NodeWidget(showtoolbar=False)
        # Build defualt node library
        # self.nodeChart.nodeLibrary.buildDefualtLibrary()
        # self.ui.horizontalLayoutPyqtnode.addWidget(self.nodeChart)

        # --- signals ---
        for mode, btn in self.modebuttondict.items():
            btn.pressed.connect(make_callback(self.mode_changed, mode))

        self.ui.treeWidgetModelNavigation.itemSelectionChanged.connect(self.navigationChanged)

        # --- default ---
        self.ui.pushButtonModeler.setChecked(True)
        self.mode_changed('modeler')
        top = self.ui.treeWidgetModelNavigation.topLevelItem(0)
        self.ui.treeWidgetModelNavigation.setCurrentItem(top)

        self.open_project(self.get_project_dir())

    def __setup_vtk_widget(self):
        " setup the vtk widget "

        self.vtkwidget = VtkWidget()

        # --- geometry buttons ---
        self.add_geometry_menu = QtGui.QMenu(self)
        self.ui.toolButtonAddGeometry.setMenu(self.add_geometry_menu)

        action = QtGui.QAction('STL File',  self.add_geometry_menu)
        action.triggered.connect(self.vtkwidget.add_stl)
        self.add_geometry_menu.addAction(action)

        self.add_geometry_menu.addSeparator()

        for geo in self.vtkwidget.primitivedict.keys():
            action = QtGui.QAction(geo,  self.add_geometry_menu)
            action.triggered.connect(
                make_callback(self.vtkwidget.add_primitive, geo))
            self.add_geometry_menu.addAction(action)

        self.ui.treeWidgetGeometry.setStyleSheet(
            "QTreeView::indicator:unchecked {image: url(%s);}"
            "QTreeView::indicator:checked {image: url(%s);}"
            % (get_image_path('visibilityofftransparent.png'),
               get_image_path('visibility.png'))
            )

        # setup signals
        self.ui.toolButtonRemoveGeometry.pressed.connect(
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
        for i in range(self.ui.stackedWidgetMode.count()):
            widget = self.ui.stackedWidgetMode.widget(i)
            if mode == str(widget.objectName()):
                current_index = i
                break

        self.ui.stackedWidgetMode.setCurrentIndex(current_index)

        for key, btn in self.modebuttondict.items():
            btn.setChecked(mode == key)

    def navigation_changed(self):
        current_selection = self.ui.treeWidgetModelNavigation.selectedItems()

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
        mfix_home = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        dmp = '--dmp' if self.ui.dmp_button.isChecked() else ''
        smp = '--smp' if self.ui.smp_button.isChecked() else ''
        build_cmd = os.path.join(mfix_home, 'configure_mfix %s %s && make pymfix' % (smp, dmp))
        self.build_thread.start_command(build_cmd, self.get_project_dir())

    def clear_output(self):
        """ build mfix """
        self.run_thread.start_command('echo "Removing:";'
                                      ' ls *.LOG *.OUT *.RES *.SP? *.pvd *vtp;'
                                      ' rm -f *.LOG *.OUT *.RES *.SP? *.pvd *vtp',
                                      self.get_project_dir())

    def run_mfix(self):
        """ build mfix """
        pymfix_exe = os.path.join(self.get_project_dir(), 'pymfix')
        self.run_thread.start_command(pymfix_exe, self.get_project_dir())

    def connect_mfix(self):
        """ build mfix """
        self.ui.mfix_browser.load(QUrl("http://localhost:5000"))

    def open_project(self, project_dir=None):
        """
        Open MFiX Project
        """

        if not project_dir:
            project_dir = str(
                 QtGui.QFileDialog.getExistingDirectory(self, 'Open Project Directory',
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
            #check to see if file is already open
            # if not self.codeEditor.is_file_opened(mfix_dat):
                # self.codeEditor.load(mfix_dat)
            # parse and update widget values
            # self.setWidgetInit(False)
            # self.projectManager.loadMfixDat(mfix_dat)
            # self.setWidgetInit(True)

            src = open(mfix_dat).read()
            self.ui.mfix_dat_source.setPlainText(src)
            self.mode_changed('developer')
            # self.ui.stackedWidgetMode.setCurrentIndex(2)
        else:
            print("mfix.dat doesn't exist")
            # self.newMfixDat()

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
            popen = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, shell=True, cwd=self.cwd)
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
    # mfix.vtkiren.Initialize()

    # exit with Ctrl-C at the terminal
    timer = QTimer()
    timer.start(500)  # You may change this if you wish.
    timer.timeout.connect(lambda: None)  # Let the interpreter run each 500 ms.
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    qapp.exec_()

    qapp.deleteLater()
    sys.exit()
