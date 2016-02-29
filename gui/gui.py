# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals
import os
import sys
import subprocess
import re

SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname( __file__ ),))
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))
print(SCRIPT_DIRECTORY)
# Qt Imports
from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtCore import QObject, QThread, pyqtSignal, QUrl

# VTK imports
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

# local imports
# from pyqtnode import NodeWidget

def get_image_path(name, default=None):
    path = os.path.join(SCRIPT_DIRECTORY,'icons',name)

    if os.name == 'nt':
        path = path.replace('\\', '//')

    return path

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
        icon_path = get_image_path(name, default=None)
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
    if base in listofstrings:
        # look for number at end
        nums = re.findall('[\d]+',base)
        if nums:
            number = int(nums[-1]) + 1
            base = base.replace(nums[-1], '')
        else:
            number=1

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

        self.setWindowTitle('MFIX')
        self.setWindowIcon(get_icon('mfix.png'))

        # --- data ---
        self.currentDirectory = './'

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

        self.ui.toolButtonOpen.setIcon(get_icon('openfolder.png'))

        # --- icons ---
        self.ui.toolButtonNew.setIcon(get_icon('newfolder.png'))
        self.ui.toolButtonOpen.setIcon(get_icon('openfolder.png'))
        self.ui.toolButtonSave.setIcon(get_icon('save.png'))
        self.ui.toolButtonOpen.pressed.connect(self.open_project)


        self.ui.build_mfix_button.pressed.connect(self.build_mfix)
        self.ui.run_mfix_button.pressed.connect(self.run_mfix)
        self.ui.connect_mfix_button.pressed.connect(self.connect_mfix)

        self.build_thread = BuildThread(self)
        self.run_thread = RunThread(self)

        def make_handler(qtextbrowser):
            def handle_line(line):
                print(qtextbrowser.objectName()," : ",line)
                cursor = qtextbrowser.textCursor()
                cursor.movePosition(cursor.End)
                cursor.insertText(line)
                qtextbrowser.ensureCursorVisible()

            return handle_line

        self.build_thread.line_printed.connect(make_handler(self.ui.build_output))

        self.run_thread.line_printed.connect(make_handler(self.ui.run_output))

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
        self.addGeomatryMenu = QtGui.QMenu(self)
        self.ui.toolButtonAddGeometry.setMenu(self.addGeomatryMenu)
        self.ui.toolButtonAddGeometry.setPopupMode(QtGui.QToolButton.InstantPopup)

        action = QtGui.QAction('STL File',  self.addGeomatryMenu)
        action.triggered.connect(self.addStl)
        self.addGeomatryMenu.addAction(action)

        self.addGeomatryMenu.addSeparator()

        for geo in self.primitivedict.keys():
            action = QtGui.QAction(geo,  self.addGeomatryMenu)
            action.triggered.connect(self._makeCallback(self.addPrimitive, geo))
            self.addGeomatryMenu.addAction(action)

        self.ui.treeWidgetGeometry.setStyleSheet(
            "QTreeView::indicator:unchecked {image: url("+get_image_path('visibilityofftransparent.png')+");}"
            "QTreeView::indicator:checked {image: url("+get_image_path('visibility.png')+");}"
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
            btn.pressed.connect(self._makeCallback(self.modeChanged, mode))

        self.ui.treeWidgetModelNavigation.itemSelectionChanged.connect(self.navigationChanged)

        self.ui.toolButtonRemoveGeometry.pressed.connect(self.removeGeometry)
        self.ui.treeWidgetGeometry.itemSelectionChanged.connect(self.treeWidgetGeometryChanged)
        self.ui.treeWidgetGeometry.itemClicked.connect(self.geometryClicked)

        # connect boolean
        for key, btn in self.booleanbtndict.items():
            btn.pressed.connect(self._makeCallback(self.booleanOperation,key))

        # connect primitive lineedits
        for widget in widget_iter(self.ui.stackedWidgetGeometryDetails):
            if isinstance(widget, QtGui.QLineEdit):
                widget.editingFinished.connect(self._makeCallback(self.primitiveEdited, widget))

        # --- default ---
        self.ui.pushButtonModeler.setChecked(True)
        self.modeChanged('modeler')
        self.ui.treeWidgetModelNavigation.setCurrentItem(self.ui.treeWidgetModelNavigation.topLevelItem(0))

    def _makeCallback(self, func, param):
        '''
        Helper function to make sure lambda functions are cached and not lost.
        '''
        return lambda: func(param)

    def modeChanged(self, mode):

        for i in range(self.ui.stackedWidgetMode.count()):
            widget = self.ui.stackedWidgetMode.widget(i)
            if mode == str(widget.objectName()):
                break

        self.ui.stackedWidgetMode.setCurrentIndex(i)

        for key, btn in self.modebuttondict.items():
            btn.setChecked(mode == key)

    def navigationChanged(self):
        currentSelection = self.ui.treeWidgetModelNavigation.selectedItems()

        if currentSelection:

            text = str(currentSelection[-1].text(0))
            text = '_'.join(text.lower().split(' '))

            for i in range(self.ui.stackedWidgetTaskPane.count()):
                widget = self.ui.stackedWidgetTaskPane.widget(i)
                if text == str(widget.objectName()):
                    break

            self.ui.stackedWidgetTaskPane.setCurrentIndex(i)

    def treeWidgetGeometryChanged(self):

        currentSelection = self.ui.treeWidgetGeometry.selectedItems()


        enableboolbtn = False
        if len(currentSelection)==2 and all(
                [self.ui.treeWidgetGeometry.indexOfTopLevelItem(select)>-1 for select in currentSelection]
                ):

            enableboolbtn = True

        for btn in self.booleanbtndict.values():
            btn.setEnabled(enableboolbtn)

        if currentSelection:
            text = str(currentSelection[-1].text(0)).lower()
            self.ui.toolButtonRemoveGeometry.setEnabled(True)

            for i in range(self.ui.stackedWidgetGeometryDetails.count()):
                widget = self.ui.stackedWidgetGeometryDetails.widget(i)
                if str(widget.objectName()) in text:
                    break

            for child in widget_iter(widget):
                if isinstance(child, QtGui.QLineEdit):
                    for key, value in self.geometrydict[text].items():
                        if key in str(child.objectName()).lower():
                            child.setText(str(value))

            self.ui.groupBoxGeometryParameters.setTitle(text)

            self.ui.stackedWidgetGeometryDetails.setCurrentIndex(i)
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

        filename = str(QtGui.QFileDialog.getOpenFileName(
            self, 'Select an STL File',
            self.currentDirectory,
            'STL File (*.stl)',))

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
            self.geometrydict[name]={
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

    def primitiveEdited(self, widget):

        currentSelection = self.ui.treeWidgetGeometry.selectedItems()
        parameter = str(widget.objectName()).lower()

        if currentSelection:
            name = str(currentSelection[-1].text(0)).lower()

            for key in self.geometrydict[name].keys():
                if key in parameter:
                    string = str(widget.text())
                    if 'resolution' in parameter:
                        self.geometrydict[name][key] = int(string)
                    else:
                        self.geometrydict[name][key] = float(string)

            self.updatePrimitive(name)
            self.updateTransform(name)
            self.vtkRenderWindow.Render()

    def updatePrimitive(self, name):
        primtype = self.geometrydict[name]['type']
        source = self.geometrydict[name]['source']

        if primtype =='sphere':
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

    def updateTransform(self, name):
        transform = self.geometrydict[name]['transform']
        transformFilter = self.geometrydict[name]['transformFilter']

        # rotation
        transform.RotateWXYZ(self.geometrydict[name]['rotationx'],1,0,0)
        transform.RotateWXYZ(self.geometrydict[name]['rotationy'],0,1,0)
        transform.RotateWXYZ(self.geometrydict[name]['rotationz'],0,0,1)
#        transformFilter.SetTransform(transform)
        transformFilter.Update()

        return transformFilter

    def addPrimitive(self, primtype='sphere'):

        name = get_unique_string(primtype, list(self.geometrydict.keys()))

        # create primative
        if primtype in self.primitivedict:
            source = self.primitivedict[primtype]()
        else:
            return

        self.geometrydict[name]={
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

        source = self.updatePrimitive(name)

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

        self.updateTransform(name)


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


    def booleanOperation(self, booltype):

        currentSelection = self.ui.treeWidgetGeometry.selectedItems()

        if len(currentSelection)==2:

            # Save references
            boolname = get_unique_string(booltype, list(self.geometrydict.keys()))
            self.geometrydict[boolname]={
                'type':     booltype,
                'children': [],
            }

            booleanOperation = vtk.vtkBooleanOperationPolyDataFilter()

            if booltype == 'union':
                booleanOperation.SetOperationToUnion()
            elif booltype == 'intersection':
                booleanOperation.SetOperationToIntersection()
            else:
                booleanOperation.SetOperationToDifference()

            for i, selection in enumerate(currentSelection):
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

#                if vtk.VTK_MAJOR_VERSION <= 5:
                booleanOperation.SetInputConnection(i, geometry.GetOutputPort())
#                else:
#                    booleanOperation.SetInputData(i, geometry.GetOutput())

            booleanOperation.Update()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(booleanOperation.GetOutputPort())
            mapper.ScalarVisibilityOff()

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetRepresentationToWireframe()

            self.vtkrenderer.AddActor(actor)

            self.vtkRenderWindow.Render()

            # save references
            self.geometrydict[boolname]['booleanoperation'] = booleanOperation
            self.geometrydict[boolname]['mapper'] = mapper
            self.geometrydict[boolname]['actor'] = actor

            # Add to tree
            toplevel = QtGui.QTreeWidgetItem([boolname])
            toplevel.setFlags(toplevel.flags() | QtCore.Qt.ItemIsUserCheckable)
            toplevel.setCheckState(0, QtCore.Qt.Checked)

            # remove children from tree
            for select in currentSelection:
                toplevelindex = self.ui.treeWidgetGeometry.indexOfTopLevelItem(select)
                item = self.ui.treeWidgetGeometry.takeTopLevelItem(toplevelindex)
                item.setCheckState(0, QtCore.Qt.Unchecked)
                toplevel.addChild(item)

            self.ui.treeWidgetGeometry.addTopLevelItem(toplevel)
            self.ui.treeWidgetGeometry.setCurrentItem(toplevel)

    def removeGeometry(self):

        currentSelection = self.ui.treeWidgetGeometry.selectedItems()
        if currentSelection:
            text = str(currentSelection[-1].text(0)).lower()

    def build_mfix(self):
        """ build mfix """
        curr_dir = os.path.dirname(os.path.realpath(__file__))
        self.build_thread.start_command(os.path.join(curr_dir, '../configure_mfix && make -j pymfix'))

    def run_mfix(self):
        """ build mfix """
        curr_dir = os.path.dirname(os.path.realpath(__file__))
        self.run_thread.start_command(os.path.join(curr_dir, 'pymfix'))

    def connect_mfix(self):
        """ build mfix """
        self.ui.mfix_browser.load(QUrl("http://localhost:5000"))

    def open_project(self, project_dir=None):
        """
        Open MFiX Project
        """

        if not project_dir:
            project_dir = str(QtGui.QFileDialog.getExistingDirectory(self, 'Open Project Directory',
                                                                       "",
                                                                       QtGui.QFileDialog.ShowDirsOnly,
                                                                       ))
        if len(project_dir)<1:return

        writable = True
        try:
            import tempfile
            testfile = tempfile.TemporaryFile(project_dir)
            testfile.close()
        except IOError:
            writable = False

        if not writable:
            self.message(title = 'Warning',
                         icon = 'warning',
                         text = ('You do not have write access to this '
                                 'directory.\n'
                                 'Please fix and try again.'),
                         buttons = ['ok'],
                         default = 'ok',
                         )
            mylogger.debug('User does not have write access: {}'.format(project_dir))
            return

        self.project_dir = project_dir

        # if not self.projectDir == project_dir:
        #     self.codeEditor.closeAll()

        self.projectDir = project_dir

        # self.recentProjects.append(self.projectDir)
        # self.updateRecentProjectMenu()

        self.setWindowTitle('MFIX - %s' % self.projectDir)

        # # look for gui data file
        # guidatafile = glob.glob(os.path.join(self.projectDir, self.guiDataFileName))
        # if guidatafile:
        #     mylogger.debug('Opening gui data file: {}'.format(guidatafile[0]))
        #     with open(guidatafile[0], 'r') as f:
        #         self.projectGuiData.update(json.load(f))

        # Look for mfix.dat file
        mfixDat = os.path.abspath(os.path.join(self.projectDir,'mfix.dat'))

        if os.path.exists(mfixDat):
            # mylogger.debug('found mfix.dat: {}'.format(mfixDat))
            #check to see if file is already open
            # if not self.codeEditor.is_file_opened(mfixDat):
                # self.codeEditor.load(mfixDat)
            # parse and update widget values
            # self.setWidgetInit(False)
            # self.projectManager.loadMfixDat(mfixDat)
            # self.setWidgetInit(True)

            src = open(mfixDat).read()
            self.ui.mfix_dat_source.setPlainText(src)
            self.modeChanged('developer')
            # self.ui.stackedWidgetMode.setCurrentIndex(2)
        else:
            print("mfix.dat doesn't exist")
            # self.newMfixDat()

        return

        # #Look for usr fortran files
        # usrFiles = glob.glob(os.path.join(self.projectDir,'*.f'))
        # if usrFiles:
        #     for usrFile in usrFiles:
        #         if not self.codeEditor.is_file_opened(usrFile):
        #             self.codeEditor.load(usrFile)

        # if mfixDat:
        #     self.codeEditor.go_to_file(mfixDat, 0, None)

        # self.editorWindow.showWindow()
        # self.enableWidgets()

        # Look for stl file
        stlFileList = glob.glob(os.path.join(self.projectDir, 'geometry.stl'))
        if stlFileList:
            self.openStlFile(stlFileList[0])

        # Look for log file
        runlogPath = os.path.join(self.projectDir, self.CONF.get('run','logfilename'))
        if os.path.exists(runlogPath):

            if hasattr(self, 'logwatcher'):
                self.logwatcher.cancel()

            self.logwatcher = WatchRunDir(self.projectDir, logfname=self.CONF.get('run','logfilename'), updatetime = float(self.CONF.get('run','refreshrate')))
            self.logwatcher.update.connect(self.logChanged)
            self.logwatcher.finished.connect(self.runfinished)
            self.logwatcher.start()

            self.running = True
            self.stalled = False

            self.progressBar.setMaximum(100)
            self.progressBar.setMinimum(0)
            self.progressBar.setValue(0)
            self.progressBar.setFormat('Running %p%')
            self.progressBar.cancel.connect(lambda:self.cancelRun())
            self.progressBar.show()

            self.mfixprocesstimer.start(self.CONF.get('run', 'checkrunning'))

        os.chdir(self.projectDir)
        mylogger.debug('CWD: {}'.format(os.getcwd()))
        self.gridInputChange()
        self.updateDocs()

class RunThread(QThread):

    line_printed = pyqtSignal(str)

    def __init__(self, parent):
        super(RunThread, self).__init__(parent)
        self.cmd = None

    def start_command(self, cmd):
        self.cmd = cmd
        self.start()

    def run(self):
        if self.cmd:
            popen = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, shell=True)
            lines_iterator = iter(popen.stdout.readline, b"")
            for line in lines_iterator:
                self.line_printed.emit(line)

class BuildThread(QThread):

    line_printed = pyqtSignal(str)

    def __init__(self, parent):
        super(BuildThread, self).__init__(parent)
        self.cmd = None

    def start_command(self, cmd):
        self.cmd = cmd
        self.start()

    def run(self):
        if self.cmd:
            popen = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, shell=True)
            lines_iterator = iter(popen.stdout.readline, b"")
            for line in lines_iterator:
                self.line_printed.emit(line)


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    mfix = MfixGui(app)

    mfix.show()
    mfix.vtkiren.Initialize()

    app.exec_()

    app.deleteLater()
    sys.exit()
