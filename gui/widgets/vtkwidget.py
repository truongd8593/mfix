# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

import os

# Qt imports
from qtpy import QtCore, QtGui

# VTK imports
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

# local imports
from tools.general import get_unique_string, widget_iter


class CustomInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self, parent=None):
        self.AddObserver("LeftButtonPressEvent", self.left_button_press_event)

        self.LastPickedActor = None
        self.LastPickedProperty = vtk.vtkProperty()

    def left_button_press_event(self, obj, event):
        clickPos = self.GetInteractor().GetEventPosition()

        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())

        # get the new actor
        self.NewPickedActor = picker.GetActor()

        # If something was selected
        if self.NewPickedActor:
            # If we picked something before, reset its property
            if self.LastPickedActor:
                self.LastPickedActor.GetProperty().DeepCopy(
                    self.LastPickedProperty)

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
            self.LastPickedActor.GetProperty().DeepCopy(
                self.LastPickedProperty)
            self.LastPickedActor = None

        self.OnLeftButtonDown()
        return


class VtkWidget(QtGui.QWidget):
    def __init__(self,parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.parent = parent
        self.geometrytree = self.parent.ui.treeWidgetGeometry
        self.booleanbtndict = self.parent.booleanbtndict

        # --- data ---
        self.geometrydict = {}

        self.primitivedict = {}
        self.primitivedict = {
            'sphere': vtk.vtkSphereSource,
            'box':    vtk.vtkCubeSource,
            'cylinder': vtk.vtkCylinderSource,
            'cone': vtk.vtkConeSource,
            }

        # --- layout ---
        self.hlayout = QtGui.QHBoxLayout(self)
        self.setLayout(self.hlayout)

        self.vtkWindowWidget = QVTKRenderWindowInteractor(self)
        self.hlayout.addWidget(self.vtkWindowWidget)

        # --- setup vtk stuff ---
        self.vtkrenderer = vtk.vtkRenderer()
        self.vtkrenderer.GradientBackgroundOn()
        self.vtkrenderer.SetBackground(.4, .4, .4)
        self.vtkrenderer.SetBackground2(1.0, 1.0, 1.0)

        self.vtkRenderWindow = self.vtkWindowWidget.GetRenderWindow()
        self.vtkRenderWindow.AddRenderer(self.vtkrenderer)
        self.vtkiren = self.vtkWindowWidget.GetRenderWindow().GetInteractor()

        self.style = CustomInteractorStyle()
        self.style.SetDefaultRenderer(self.vtkrenderer)
        self.vtkiren.SetInteractorStyle(self.style)

        # Orientation Marker Widget
        self.axes = vtk.vtkAxesActor()
        self.axes.AxisLabelsOn()
        self.axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(
            1, 0, 0)
        self.axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
        self.axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(
            0, 1, 0)
        self.axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
        self.axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(
            0, 0, 1)
        self.axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()

        self.orientationWidget = vtk.vtkOrientationMarkerWidget()
        self.orientationWidget.SetOutlineColor(0.9300, 0.5700, 0.1300)
        self.orientationWidget.SetOrientationMarker(self.axes)
        self.orientationWidget.SetInteractor(self.vtkiren)
        self.orientationWidget.SetViewport(0.0, 0.0, 0.2, 0.2)
        self.orientationWidget.SetEnabled(1)
        self.orientationWidget.InteractiveOff()

        self.vtkrenderer.ResetCamera()

        # connect events
        self.geometrytree.itemSelectionChanged.connect(
            self.tree_widget_geometry_changed)
        self.geometrytree.itemClicked.connect(self.geometry_clicked)

    def tree_widget_geometry_changed(self):

        current_selection = self.geometrytree.selectedItems()

        enableboolbtn = False
        if len(current_selection) == 2 and all(
                [self.geometrytree.indexOfTopLevelItem(select) > -1
                 for select in current_selection]
        ):

            enableboolbtn = True

        for btn in self.booleanbtndict.values():
            btn.setEnabled(enableboolbtn)

        if current_selection:
            text = str(current_selection[-1].text(0)).lower()
            self.parent.ui.toolButtonRemoveGeometry.setEnabled(True)

            current_index = 0
            for i in range(self.parent.ui.stackedWidgetGeometryDetails.count()):
                widget = self.parent.ui.stackedWidgetGeometryDetails.widget(i)
                if str(widget.objectName()) in text:
                    current_index = i
                    break

            for child in widget_iter(widget):
                if isinstance(child, QtGui.QLineEdit):
                    for key, value in self.geometrydict[text].items():
                        if key in str(child.objectName()).lower():
                            child.setText(str(value))

            self.parent.ui.groupBoxGeometryParameters.setTitle(text)

            self.parent.ui.stackedWidgetGeometryDetails.setCurrentIndex(current_index)
        else:
            self.parent.ui.stackedWidgetGeometryDetails.setCurrentIndex(
                self.parent.ui.stackedWidgetGeometryDetails.count()-1
                )

            self.parent.ui.groupBoxGeometryParameters.setTitle('Parameters')
            self.parent.ui.toolButtonRemoveGeometry.setEnabled(False)

    def geometry_clicked(self, item):

        name = str(item.text(0)).lower()
        if item.checkState(0) == QtCore.Qt.Unchecked:
            self.geometrydict[name]['actor'].VisibilityOff()
        else:
            self.geometrydict[name]['actor'].VisibilityOn()

        self.vtkRenderWindow.Render()

    def add_stl(self):

        filename = str(QtGui.QFileDialog.getOpenFileName(
            self, 'Select an STL File',
            self.currentDirectory,
            'STL File (*.stl)',))

        if filename:
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
            self.geometrytree.addTopLevelItem(item)
            self.geometrytree.setCurrentItem(item)

    def primitive_edited(self, widget):

        current_selection = self.geometrytree.selectedItems()
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
            # self.vtkRenderWindow.Render()

    def update_primitive(self, name):
        primtype = self.geometrydict[name]['type']
        source = self.geometrydict[name]['source']

        if primtype == 'sphere':
            # Create source
            source.SetRadius(self.geometrydict[name]['radius'])
            source.SetThetaResolution(
                self.geometrydict[name]['thetaresolution'])
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
            'centerx':         0.0,
            'centery':         0.0,
            'centerz':         0.0,
            'rotationx':       0.0,
            'rotationy':       0.0,
            'rotationz':       0.0,
            'radius':          1.0,
            'directionx':      1.0,
            'directiony':      0.0,
            'directionz':      0.0,
            'lengthx':         1.0,
            'lengthy':         1.0,
            'lengthz':         1.0,
            'height':          1.0,
            'resolution':      10,
            'thetaresolution': 10,
            'phiresolution':   10,
            'type':            primtype,
            'source':          source,
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
        self.geometrytree.addTopLevelItem(item)
        self.geometrytree.setCurrentItem(item)

    def boolean_operation(self, booltype):

        current_selection = self.geometrytree.selectedItems()

        if len(current_selection) == 2:

            # Save references
            boolname = get_unique_string(booltype,
                                         list(self.geometrydict.keys()))
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

#                if vtk.VTK_MAJOR_VERSION <= 5:
                boolean_operation.SetInputConnection(i,
                                                     geometry.GetOutputPort())
#                else:
#                    boolean_operation.SetInputData(i, geometry.GetOutput())

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
                toplevelindex = self.geometrytree.indexOfTopLevelItem(select)
                item = self.geometrytree.takeTopLevelItem(toplevelindex)
                item.setCheckState(0, QtCore.Qt.Unchecked)
                toplevel.addChild(item)

            self.geometrytree.addTopLevelItem(toplevel)
            self.geometrytree.setCurrentItem(toplevel)

    def remove_geometry(self):
        currentSelection = self.geometrytree.selectedItems()
        if currentSelection:
            text = str(currentSelection[-1].text(0)).lower()

            # remove tree item
            toplevelindex = self.geometrytree.indexOfTopLevelItem(
                currentSelection[-1])
            item = self.geometrytree.takeTopLevelItem(toplevelindex)

            # move children to toplevel, make visible
            for child in item.takeChildren():
                self.geometrytree.addTopLevelItem(child)
                self.geometrydict[
                    str(child.text(0)).lower()]['actor'].VisibilityOn()
                child.setCheckState(0, QtCore.Qt.Checked)

            # remove graphics
            geo = self.geometrydict.pop(text)
            self.vtkrenderer.RemoveActor(geo['actor'])

            self.vtkRenderWindow.Render()
