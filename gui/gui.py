""" MFIX GUI """


# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals
import logging
import os
import shutil
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

        self.ui.toolbutton_add_geometry.setIcon(get_icon('add.png'))
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

<<<<<<< HEAD
        self.ui.toolButtonTFMSolidsSpeciesAdd.setIcon(get_icon('add.png'))
        self.ui.toolButtonTFMSolidsSpeciesDelete.setIcon(
            get_icon('remove.png'))
        self.ui.toolButtonTFMSolidsSpeciesCopy.setIcon(get_icon('copy.png'))

        # --- Connect Signals to Slots---
        # open/save/new project
        self.ui.toolbutton_open.pressed.connect(self.open_project)
        self.ui.toolbutton_save.pressed.connect(self.save_project)
=======
        # --- icons ---
        self.ui.toolButtonNew.setIcon(get_icon('newfolder.png'))
        self.ui.toolButtonOpen.setIcon(get_icon('openfolder.png'))
        self.ui.toolButtonSave.setIcon(get_icon('save.png'))
        self.ui.toolButtonNew.pressed.connect(self.new_project)
        self.ui.toolButtonOpen.pressed.connect(self.open_project)
        self.ui.toolButtonSave.pressed.connect(self.save_project)
>>>>>>> b9e087e375f994b7cd206bb15ef419f089d9bf06

        # mode (modeler, workflow, developer)
        for mode, btn in self.modebuttondict.items():
            btn.released.connect(make_callback(self.mode_changed, mode))

        # navigation tree
        self.ui.treeWidgetModelNavigation.itemSelectionChanged.connect(
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
        top = self.ui.treeWidgetModelNavigation.topLevelItem(0)
        self.ui.treeWidgetModelNavigation.setCurrentItem(top)

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

<<<<<<< HEAD
    def __setup_workflow_widget(self):

        self.nodeChart = NodeWidget(showtoolbar=False)
        # Build defualt node library
        self.nodeChart.nodeLibrary.buildDefualtLibrary()
        self.ui.horizontalLayoutPyqtnode.addWidget(self.nodeChart)
=======
        last_project = self.get_project_dir()
        if os.path.exists(os.path.join(last_project,'mfix.dat')):
            self.open_project(last_project)
>>>>>>> b9e087e375f994b7cd206bb15ef419f089d9bf06

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

<<<<<<< HEAD
=======
    def treeWidgetGeometryChanged(self):

        current_selection = self.ui.treeWidgetGeometry.selectedItems()


        enableboolbtn = False
        if len(current_selection) == 2 and all(
                [self.ui.treeWidgetGeometry.indexOfTopLevelItem(select) > -1
                 for select in current_selection]
        ):

            enableboolbtn = True

        for btn in self.booleanbtndict.values():
            btn.setEnabled(enableboolbtn)

        if current_selection:
            text = str(current_selection[-1].text(0)).lower()
            self.ui.toolButtonRemoveGeometry.setEnabled(True)

            current_index = 0
            for i in range(self.ui.stackedWidgetGeometryDetails.count()):
                widget = self.ui.stackedWidgetGeometryDetails.widget(i)
                if str(widget.objectName()) in text:
                    current_index = i
                    break

            for child in widget_iter(widget):
                if isinstance(child, QtGui.QLineEdit):
                    for key, value in self.geometrydict[text].items():
                        if key in str(child.objectName()).lower():
                            child.setText(str(value))

            self.ui.groupBoxGeometryParameters.setTitle(text)

            self.ui.stackedWidgetGeometryDetails.setCurrentIndex(current_index)
        else:
            self.ui.stackedWidgetGeometryDetails.setCurrentIndex(
                self.ui.stackedWidgetGeometryDetails.count()-1
                )

            self.ui.groupBoxGeometryParameters.setTitle('Parameters')
            self.ui.toolButtonRemoveGeometry.setEnabled(False)

    def geometryClicked(self, item):

        name = str(item.text(0)).lower()
        if item.checkState(0) == QtCore.Qt.Unchecked:
            self.geometrydict[name]['actor'].VisibilityOff()
        else:
            self.geometrydict[name]['actor'].VisibilityOn()

        self.vtkRenderWindow.Render()


    def addStl(self):

        # filename = str(QtGui.QFileDialog.getOpenFileName(
        #     self, 'Select an STL File',
        #     self.currentDirectory,
        #     'STL File (*.stl)',))

        if filename:
            # TODO: check to see if name is unique in self.geometrydict
            name = os.path.basename(filename).lower()
            name = get_unique_string(name, list(self.geometrydict.keys()))

            reader = vtk.vtkSTLReader()
            reader.SetFileName(filename)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(reader.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetRepresentationToWireframe()

            self.vtkrenderer.AddActor(actor)

            self.vtkrenderer.ResetCamera()
            self.vtkRenderWindow.Render()

            # Add to dict
            self.geometrydict[name] = {
                'reader': reader,
                'mapper': mapper,
                'actor':  actor,
                'type':   'stl',
                'centerx': 0.0,
                'centery': 0.0,
                'centerz': 0.0,
                }

            # Add to tree
            item = QtGui.QTreeWidgetItem([name])
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(0, QtCore.Qt.Checked)
            self.ui.treeWidgetGeometry.addTopLevelItem(item)
            self.ui.treeWidgetGeometry.setCurrentItem(item)

    def primitive_edited(self, widget):

        current_selection = self.ui.treeWidgetGeometry.selectedItems()
        parameter = str(widget.objectName()).lower()

        if current_selection:
            name = str(current_selection[-1].text(0)).lower()

            for key in self.geometrydict[name].keys():
                if key in parameter:
                    string = str(widget.text())
                    if 'resolution' in parameter:
                        self.geometrydict[name][key] = int(string)
                    else:
                        self.geometrydict[name][key] = float(string)

            self.update_primitive(name)
            self.update_transform(name)
            self.vtkRenderWindow.Render()

    def update_primitive(self, name):
        primtype = self.geometrydict[name]['type']
        source = self.geometrydict[name]['source']

        if primtype == 'sphere':
            # Create source
            source.SetRadius(self.geometrydict[name]['radius'])
            source.SetThetaResolution(self.geometrydict[name]['thetaresolution'])
            source.SetPhiResolution(self.geometrydict[name]['phiresolution'])

        elif primtype == 'box':
            source.SetXLength(self.geometrydict[name]['lengthx'])
            source.SetYLength(self.geometrydict[name]['lengthy'])
            source.SetZLength(self.geometrydict[name]['lengthz'])

        elif primtype == 'cone':
            source.SetRadius(self.geometrydict[name]['radius'])
            source.SetHeight(self.geometrydict[name]['height'])
            source.SetDirection(self.geometrydict[name]['directionx'],
                                self.geometrydict[name]['directiony'],
                                self.geometrydict[name]['directionz'])
            source.SetResolution(self.geometrydict[name]['resolution'])
            source.CappingOn()

        elif primtype == 'cylinder':
            source.SetRadius(self.geometrydict[name]['radius'])
            source.SetHeight(self.geometrydict[name]['height'])
            source.SetResolution(self.geometrydict[name]['resolution'])

        else:
            return

        # common props
        source.SetCenter(self.geometrydict[name]['centerx'],
                         self.geometrydict[name]['centery'],
                         self.geometrydict[name]['centerz'])

        source.Update()

        return source

    def update_transform(self, name):
        transform = self.geometrydict[name]['transform']
        transform_filter = self.geometrydict[name]['transformFilter']

        # rotation
        transform.RotateWXYZ(self.geometrydict[name]['rotationx'], 1, 0, 0)
        transform.RotateWXYZ(self.geometrydict[name]['rotationy'], 0, 1, 0)
        transform.RotateWXYZ(self.geometrydict[name]['rotationz'], 0, 0, 1)
#        transformFilter.SetTransform(transform)
        transform_filter.Update()

        return transform_filter

    def add_primitive(self, primtype='sphere'):

        name = get_unique_string(primtype, list(self.geometrydict.keys()))

        # create primative
        if primtype in self.primitivedict:
            source = self.primitivedict[primtype]()
        else:
            return

        self.geometrydict[name] = {
            'centerx':  0.0,
            'centery':  0.0,
            'centerz':  0.0,
            'rotationx': 0.0,
            'rotationy': 0.0,
            'rotationz': 0.0,
            'radius':  1.0,
            'directionx': 1.0,
            'directiony': 0.0,
            'directionz': 0.0,
            'lengthx': 1.0,
            'lengthy': 1.0,
            'lengthz': 1.0,
            'height':  1.0,
            'resolution': 10,
            'thetaresolution':10,
            'phiresolution':10,
            'type':primtype,
            'source': source,
        }

        source = self.update_primitive(name)

        # convert to triangles
        trianglefilter = vtk.vtkTriangleFilter()
        trianglefilter.SetInputData(source.GetOutput())

        # Create transformer
        transform = vtk.vtkTransform()
        transformFilter = vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        transformFilter.SetInputConnection(source.GetOutputPort())

        self.geometrydict[name]['transform'] = transform
        self.geometrydict[name]['transformFilter'] = transformFilter

        self.update_transform(name)


        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(transformFilter.GetOutputPort())

        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
#        actor.GetProperty().SetOpacity(0.25)
        actor.GetProperty().SetRepresentationToWireframe()

        self.vtkrenderer.AddActor(actor)

        self.vtkrenderer.ResetCamera()
        self.vtkRenderWindow.Render()

        # add to dict
        self.geometrydict[name]['mapper'] = mapper
        self.geometrydict[name]['actor'] = actor
        self.geometrydict[name]['trianglefiler'] = trianglefilter
        self.geometrydict[name]['source'] = source


        # Add to tree
        item = QtGui.QTreeWidgetItem([name])
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(0, QtCore.Qt.Checked)
        self.ui.treeWidgetGeometry.addTopLevelItem(item)
        self.ui.treeWidgetGeometry.setCurrentItem(item)


    def boolean_operation(self, booltype):

        current_selection = self.ui.treeWidgetGeometry.selectedItems()

        if len(current_selection) == 2:

            # Save references
            boolname = get_unique_string(booltype, list(self.geometrydict.keys()))
            self.geometrydict[boolname] = {
                'type':     booltype,
                'children': [],
            }

            boolean_operation = vtk.vtkBooleanOperationPolyDataFilter()

            if booltype == 'union':
                boolean_operation.SetOperationToUnion()
            elif booltype == 'intersection':
                boolean_operation.SetOperationToIntersection()
            else:
                boolean_operation.SetOperationToDifference()

            for i, selection in enumerate(current_selection):
                name = str(selection.text(0)).lower()
                self.geometrydict[boolname]['children'].append(name)

                if 'trianglefiler' in self.geometrydict[name]:
                    geometry = self.geometrydict[name]['trianglefiler']
                elif 'booleanoperation' in self.geometrydict[name]:
                    geometry = self.geometrydict[name]['booleanoperation']
                elif 'reader' in self.geometrydict[name]:
                    geometry = self.geometrydict[name]['reader']

                # hide the sources
                self.geometrydict[name]['actor'].VisibilityOff()

                if vtk.VTK_MAJOR_VERSION <= 5:
                    boolean_operation.SetInputConnection(i, geometry.GetOutputPort())
                else:
                    boolean_operation.SetInputData(i, geometry.GetOutput())

            boolean_operation.Update()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(boolean_operation.GetOutputPort())
            mapper.ScalarVisibilityOff()

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetRepresentationToWireframe()

            self.vtkrenderer.AddActor(actor)

            self.vtkRenderWindow.Render()

            # save references
            self.geometrydict[boolname]['booleanoperation'] = boolean_operation
            self.geometrydict[boolname]['mapper'] = mapper
            self.geometrydict[boolname]['actor'] = actor

            # Add to tree
            toplevel = QtGui.QTreeWidgetItem([boolname])
            toplevel.setFlags(toplevel.flags() | QtCore.Qt.ItemIsUserCheckable)
            toplevel.setCheckState(0, QtCore.Qt.Checked)

            # remove children from tree
            for select in current_selection:
                toplevelindex = self.ui.treeWidgetGeometry.indexOfTopLevelItem(select)
                item = self.ui.treeWidgetGeometry.takeTopLevelItem(toplevelindex)
                item.setCheckState(0, QtCore.Qt.Unchecked)
                toplevel.addChild(item)

            self.ui.treeWidgetGeometry.addTopLevelItem(toplevel)
            self.ui.treeWidgetGeometry.setCurrentItem(toplevel)

    def remove_geometry(self):
        " unfinished "
        pass

        # current_selection = self.ui.treeWidgetGeometry.selectedItems()
        # if current_selection:
        #     text = str(current_selection[-1].text(0)).lower()

    def message(self,
                title = 'Warning',
                icon = 'warning',
                text = 'This is a warning.',
                buttons = ['ok'],
                default = 'ok',
                infoText = None,
                detailedtext = None,
                ):
        '''
        Create a message box:
        title = 'title'
        icon = 'warning' or 'info'
        text = 'test to show'
        buttons = ['ok',...] where value is 'ok', 'yes', 'no', 'cancel', 'discard'
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

>>>>>>> b9e087e375f994b7cd206bb15ef419f089d9bf06
    def build_mfix(self):
        """ build mfix """
        mfix_home = os.path.dirname(
            os.path.dirname(os.path.realpath(__file__)))
        dmp = '--dmp' if self.ui.dmp_button.isChecked() else ''
        smp = '--smp' if self.ui.smp_button.isChecked() else ''
<<<<<<< HEAD
        build_cmd = os.path.join(
            mfix_home, 'configure_mfix %s %s && make pymfix' % (smp, dmp))
=======
        build_cmd = os.path.join(mfix_home, 'configure_mfix --python %s %s && make pymfix' % (smp, dmp))
>>>>>>> b9e087e375f994b7cd206bb15ef419f089d9bf06
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
        url = self.ui.mfix_url.text()
        qurl = QUrl(url)
        self.ui.mfix_browser.load(qurl)

        current_index = 0
        for i in range(self.ui.stackedWidgetTaskPane.count()):
            widget = self.ui.stackedWidgetTaskPane.widget(i)
            if 'job' == str(widget.objectName()):
                current_index = i
                break

        self.ui.stackedWidgetTaskPane.setCurrentIndex(current_index)
        top = self.ui.treeWidgetModelNavigation.topLevelItem(9).child(0)
        self.ui.treeWidgetModelNavigation.setCurrentItem(top)

    def save_project(self):
        project_dir = self.settings.value('project_dir')

        self.project._keywordDict['run_name'] = self.ui.run_name.text()
        self.project._keywordDict['description'] = self.ui.description.text()

        self.setWindowTitle('MFIX - %s' % project_dir)
        self.project.save(os.path.join(project_dir, 'mfix.dat'))

    def unsaved(self):
        project_dir = self.settings.value('project_dir')
        self.setWindowTitle('MFIX - %s *' % project_dir)

    def new_project(self, project_dir=None):
        if not project_dir:
            project_dir = str(
                QtGui.QFileDialog.getExistingDirectory(self, 'Create Project in Directory',
                                                       "",
                                                       QtGui.QFileDialog.ShowDirsOnly))
        try:
            shutil.copyfile('mfix.dat.template', os.path.join(project_dir, 'mfix.dat'))
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
                         text=('mfix.dat file does not exist in this directory.\n'),
                         buttons=['ok'],
                         default='ok',
            )
            return

        self.settings.setValue('project_dir', project_dir)
        self.setWindowTitle('MFIX - %s' % project_dir)

<<<<<<< HEAD
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
=======
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

        self.project = Project(mfix_dat)
        self.ui.run_name.setText(str(self.project['run_name']))
        self.ui.description.setText(str(self.project['description']))
>>>>>>> b9e087e375f994b7cd206bb15ef419f089d9bf06



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
