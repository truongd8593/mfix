""" MFIX GUI """


# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals
import logging
import os
import shutil
import signal
import sys
import subprocess
import re
import time
import urllib2

from qtpy import QtCore, QtGui
from qtpy.QtCore import QObject, QThread, pyqtSignal, QUrl, QTimer, QSettings

from tools.mfixproject import Project

try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic

# VTK imports
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

logging.basicConfig(stream=sys.stdout,
                    filemode='w', level=logging.DEBUG,
                    format='%(name)s - %(levelname)s - %(message)s')

SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), ))
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))
LOG = logging.getLogger(__name__)
LOG.debug(SCRIPT_DIRECTORY)

# local imports
# from pyqtnode import NodeWidget

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

# --- vtk stuff ---
class CustomInteractorStyleTrackballCamera(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self,parent=None):
        self.AddObserver("LeftButtonPressEvent",self.leftButtonPressEvent)

        self.LastPickedActor = None
        self.LastPickedProperty = vtk.vtkProperty()

    def leftButtonPressEvent(self,obj,event):
        clickPos = self.GetInteractor().GetEventPosition()

        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())

        # get the new
        self.NewPickedActor = picker.GetActor()

        # If something was selected
        if self.NewPickedActor:
            # If we picked something before, reset its property
            if self.LastPickedActor:
                self.LastPickedActor.GetProperty().DeepCopy(self.LastPickedProperty)


            # Save the property of the picked actor so that we can
            # restore it next time
            self.LastPickedProperty.DeepCopy(self.NewPickedActor.GetProperty())
            # Highlight the picked actor by changing its properties
            self.NewPickedActor.GetProperty().SetColor(0.0, 0.0, 1.0)
            self.NewPickedActor.GetProperty().SetDiffuse(1.0)
            self.NewPickedActor.GetProperty().SetSpecular(0.0)

            # save the last picked actor
            self.LastPickedActor = self.NewPickedActor

        # clear selection
        elif self.LastPickedActor:
            self.LastPickedActor.GetProperty().DeepCopy(self.LastPickedProperty)
            self.LastPickedActor = None


        self.OnLeftButtonDown()
        return

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

        self.geometrydict = {}

        self.primitivedict = {
            'sphere': vtk.vtkSphereSource,
            'box':    vtk.vtkCubeSource,
            'cylinder': vtk.vtkCylinderSource,
            'cone': vtk.vtkConeSource,
            }

        self.booleanbtndict = {
            'union':      self.ui.toolButtonGeometryUnion,
            'intersection':  self.ui.toolButtonGeometryIntersect,
            'difference': self.ui.toolButtonGeometryDifference,
            }

        # self.ui.toolButtonNew.setIcon(self.style().standardIcon(QtGui.QStyle.SP_FileIcon))
        # self.ui.toolButtonOpen.setIcon(self.style().standardIcon(QtGui.QStyle.SP_DialogOpenButton))
        # self.ui.toolButtonSave.setIcon(self.style().standardIcon(QtGui.QStyle.SP_DialogSaveButton))

        # --- icons ---

        self.ui.toolButtonNew.setIcon(get_icon('newfolder.png'))
        self.ui.toolButtonOpen.setIcon(get_icon('openfolder.png'))
        self.ui.toolButtonSave.setIcon(get_icon('save.png'))

        self.ui.toolButtonNew.pressed.connect(self.new_project)
        self.ui.toolButtonOpen.pressed.connect(self.open_project)
        self.ui.toolButtonSave.pressed.connect(self.save_project)

        self.ui.build_mfix_button.pressed.connect(self.build_mfix)
        self.ui.run_mfix_button.pressed.connect(self.run_mfix)
        self.ui.connect_mfix_button.pressed.connect(self.connect_mfix)
        self.ui.clear_output_button.pressed.connect(self.clear_output)

        self.ui.nodes_i.valueChanged.connect(self.set_mpi)
        self.ui.nodes_j.valueChanged.connect(self.set_mpi)
        self.ui.nodes_k.valueChanged.connect(self.set_mpi)

        self.ui.dmp_button.pressed.connect(self.set_nodes)

        self.ui.run_name.textChanged.connect(self.unsaved)
        self.ui.description.textChanged.connect(self.unsaved)

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

        self.build_thread.line_printed.connect(make_handler(self.ui.command_output))
        self.run_thread.line_printed.connect(make_handler(self.ui.command_output))
        self.clear_thread.line_printed.connect(make_handler(self.ui.command_output))

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

        # --- geometry setup ---
        self.add_geometry_menu = QtGui.QMenu(self)
        self.ui.toolButtonAddGeometry.setMenu(self.add_geometry_menu)
        self.ui.toolButtonAddGeometry.setPopupMode(QtGui.QToolButton.InstantPopup)

        action = QtGui.QAction('STL File',  self.add_geometry_menu)
        action.triggered.connect(self.addStl)
        self.add_geometry_menu.addAction(action)

        self.add_geometry_menu.addSeparator()

        for geo in self.primitivedict.keys():
            action = QtGui.QAction(geo,  self.add_geometry_menu)
            action.triggered.connect(make_callback(self.add_primitive, geo))
            self.add_geometry_menu.addAction(action)

        self.ui.treeWidgetGeometry.setStyleSheet(
            "QTreeView::indicator:unchecked {image: url(%s);}"
            "QTreeView::indicator:checked {image: url(%s);}"
            % (get_image_path('visibilityofftransparent.png'), get_image_path('visibility.png'))
            )

        # --- vtk setup ---
        self.vtkWidget = QVTKRenderWindowInteractor(self.ui.widgetModelGraphics)
        self.ui.horizontalLayoutModelGraphics.addWidget(self.vtkWidget)

        self.vtkrenderer = vtk.vtkRenderer()
        self.vtkrenderer.GradientBackgroundOn()
        self.vtkrenderer.SetBackground(.4, .4, .4)
        self.vtkrenderer.SetBackground2(1.0, 1.0, 1.0)

        self.vtkRenderWindow = self.vtkWidget.GetRenderWindow()
        self.vtkRenderWindow.AddRenderer(self.vtkrenderer)
        self.vtkiren = self.vtkWidget.GetRenderWindow().GetInteractor()

        self.style = CustomInteractorStyleTrackballCamera()
        self.style.SetDefaultRenderer(self.vtkrenderer)
        self.vtkiren.SetInteractorStyle(self.style)

        # Orientation Marker Widget
        self.axes = vtk.vtkAxesActor()
        self.axes.AxisLabelsOn()
        self.axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1,0,0)
        self.axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
        self.axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,1,0)
        self.axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
        self.axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,0,1)
        self.axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()

        self.orientationWidget = vtk.vtkOrientationMarkerWidget()
        self.orientationWidget.SetOutlineColor( 0.9300, 0.5700, 0.1300 )
        self.orientationWidget.SetOrientationMarker(self.axes)
        self.orientationWidget.SetInteractor(self.vtkiren)
        self.orientationWidget.SetViewport( 0.0, 0.0, 0.2, 0.2 )
        self.orientationWidget.SetEnabled( 1 )
        self.orientationWidget.InteractiveOff()

        self.vtkrenderer.ResetCamera()

        # --- workflow setup ---
        # self.nodeChart = NodeWidget(showtoolbar=False)
        # Build defualt node library
        # self.nodeChart.nodeLibrary.buildDefualtLibrary()
        # self.ui.horizontalLayoutPyqtnode.addWidget(self.nodeChart)

        # --- workflow setup ---
        # self.nodeChart = NodeWidget(showtoolbar=False)
        # Build defualt node library
        # self.nodeChart.nodeLibrary.buildDefualtLibrary()
        # self.ui.horizontalLayoutPyqtnode.addWidget(self.nodeChart)

        # --- signals ---
        for mode, btn in self.modebuttondict.items():
            btn.pressed.connect(make_callback(self.mode_changed, mode))

        self.ui.treeWidgetModelNavigation.itemSelectionChanged.connect(self.navigationChanged)

        self.ui.toolButtonRemoveGeometry.pressed.connect(self.remove_geometry)
        self.ui.treeWidgetGeometry.itemSelectionChanged.connect(self.treeWidgetGeometryChanged)
        self.ui.treeWidgetGeometry.itemClicked.connect(self.geometryClicked)

        # connect boolean
        for key, btn in self.booleanbtndict.items():
            btn.pressed.connect(make_callback(self.boolean_operation, key))

        # connect primitive lineedits
        for widget in widget_iter(self.ui.stackedWidgetGeometryDetails):
            if isinstance(widget, QtGui.QLineEdit):
                widget.editingFinished.connect(make_callback(self.primitive_edited, widget))

        # --- default ---
        self.ui.pushButtonModeler.setChecked(True)
        self.mode_changed('modeler')
        top = self.ui.treeWidgetModelNavigation.topLevelItem(0)
        self.ui.treeWidgetModelNavigation.setCurrentItem(top)

        last_project = self.get_project_dir()
        if last_project and os.path.exists(os.path.join(last_project,'mfix.dat')):
            self.open_project(last_project)

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

    def navigationChanged(self):
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
        mfix_home = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        dmp = '--dmp' if self.ui.dmp_button.isChecked() else ''
        smp = '--smp' if self.ui.smp_button.isChecked() else ''
        return os.path.join(mfix_home, 'configure_mfix --python %s %s && make -j pymfix' % (smp, dmp))

    def build_mfix(self):
        """ build mfix """
        self.build_thread.start_command(self.make_build_cmd(), self.get_project_dir())

    def clear_output(self):
        """ build mfix """
        self.run_thread.start_command('echo "Removing:";'
                                      ' ls *.LOG *.OUT *.RES *.SP? *.pvd *vtp;'
                                      ' rm -f *.LOG *.OUT *.RES *.SP? *.pvd *vtp',
                                      self.get_project_dir())

    def run_mfix(self):
        """ build mfix """
        if self.ui.dmp_button.isChecked():
            pymfix_exe = os.path.join(self.get_project_dir(),'pymfix')
        else:
            nodesi = int(self.ui.nodes_i.text())
            nodesj = int(self.ui.nodes_j.text())
            nodesk = int(self.ui.nodes_k.text())
            total = nodesi*nodesj*nodesk
            pymfix_exe = 'mpirun -np {} pymfix NODESI={} NODESJ={} NODESK={}'.format(total, nodesi, nodesj, nodesk)

        build_and_run_cmd = '{} && {}'.format(self.make_build_cmd(), pymfix_exe)
        self.run_thread.start_command(build_and_run_cmd, self.get_project_dir())

    def update_residuals(self):
        self.ui.residuals.setText(self.updater.residuals)

    def connect_mfix(self):
        """ connect to running instance of mfix """
        url = "http://{}:{}".format(self.ui.mfix_host.text(), self.ui.mfix_port.text())
        log = logging.getLogger(__name__)
        log.debug("trying to connect to {}".format(url))
        qurl = QUrl(url)
        self.ui.mfix_browser.load(qurl)

        self.updater = UpdateResidualsThread()
        self.updater.sig.connect(self.update_residuals)
        self.updater.start()

        current_index = 0
        for i in range(self.ui.stackedWidgetTaskPane.count()):
            widget = self.ui.stackedWidgetTaskPane.widget(i)
            if 'interact' == str(widget.objectName()):
                current_index = i
                break

        self.ui.stackedWidgetTaskPane.setCurrentIndex(current_index)
        top = self.ui.treeWidgetModelNavigation.topLevelItem(10)
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
        if len(project_dir) < 1:
            return
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

class UpdateResidualsThread(QThread):

    sig = QtCore.pyqtSignal(object)

    def run(self):
        while True:
            self.residuals = urllib2.urlopen('http://localhost:5000/residuals').read()
            time.sleep(1)
            self.sig.emit('update')

if __name__ == '__main__':
    qapp = QtGui.QApplication(sys.argv)

    mfix = MfixGui(qapp)

    mfix.show()
    mfix.vtkiren.Initialize()

    # exit with Ctrl-C at the terminal
    timer = QTimer()
    timer.start(500)  # You may change this if you wish.
    timer.timeout.connect(lambda: None)  # Let the interpreter run each 500 ms.
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    qapp.exec_()

    qapp.deleteLater()
    sys.exit()
