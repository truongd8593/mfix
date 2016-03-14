# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

import os
import copy
from collections import OrderedDict

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
            self.NewPickedActor.GetProperty().SetColor(255/255.0, 140/255.0, 0)
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
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.parent = parent
        self.geometrytree = self.parent.ui.treeWidgetGeometry
        self.booleanbtndict = self.parent.booleanbtndict

        # --- data ---
        self.geometrydict = {}

        self.primitivedict = {}
        self.primitivedict = OrderedDict([
            ('sphere',   vtk.vtkSphereSource),
            ('box',      vtk.vtkCubeSource),
            ('cylinder', vtk.vtkCylinderSource),
            ('cone',     vtk.vtkConeSource),
            ])

        self.parametricdict = OrderedDict([
            ('torus',    vtk.vtkParametricTorus),
            ('boy', vtk.vtkParametricBoy),
            ('conic_spiral', vtk.vtkParametricConicSpiral),
            ('cross_cap', vtk.vtkParametricCrossCap),
            ('dini', vtk.vtkParametricDini),
            ('ellipsoid', vtk.vtkParametricEllipsoid),
            ('enneper', vtk.vtkParametricEnneper),
            ('figure_8_klein', vtk.vtkParametricFigure8Klein),
            ('klein', vtk.vtkParametricKlein),
            ('mobius', vtk.vtkParametricMobius),
            ('random_hills', vtk.vtkParametricRandomHills),
            ('roman', vtk.vtkParametricRoman),
            ('super_ellipsoid', vtk.vtkParametricSuperEllipsoid),
            ('super_toroid', vtk.vtkParametricSuperToroid),
            ('torus', vtk.vtkParametricTorus),
            ])

        self.filterdict = OrderedDict([
            ('clean',              vtk.vtkCleanPolyData),
            ('fill_holes',         vtk.vtkFillHolesFilter),
            ('triangle',           vtk.vtkTriangleFilter),
            ('decimate',           vtk.vtkDecimatePro),
            ('quadric_decimation', vtk.vtkQuadricDecimation),
            ('quadric_clustering', vtk.vtkQuadricClustering)
        ])

        # --- layout ---
        self.hlayout = QtGui.QHBoxLayout(self)
        self.hlayout.setContentsMargins(0, 0, 0, 0)
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

    # --- geometry ---
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

        # enable/disable delete/copy/filter button
        if len(current_selection) == 1 and \
                self.parent.ui.treeWidgetGeometry.indexOfTopLevelItem(
                current_selection[0]) > -1:
            self.parent.ui.toolbutton_remove_geometry.setEnabled(True)
            self.parent.ui.toolbutton_add_filter.setEnabled(True)
            self.parent.ui.toolbutton_copy_geometry.setEnabled(True)
        else:
            self.parent.ui.toolbutton_remove_geometry.setEnabled(False)
            self.parent.ui.toolbutton_add_filter.setEnabled(False)
            self.parent.ui.toolbutton_copy_geometry.setEnabled(False)

        if current_selection:
            text = str(current_selection[-1].text(0)).lower()

            current_index = 0
            for i in range(
                    self.parent.ui.stackedWidgetGeometryDetails.count()):
                widget = self.parent.ui.stackedWidgetGeometryDetails.widget(i)
                if str(widget.objectName()) in text:
                    current_index = i
                    break

            # set the widget parameters
            for child in widget_iter(widget):
                name = str(child.objectName()).lower().replace('_', '')
                for key, value in self.geometrydict[text].items():
                    if key in name:
                        break

                if isinstance(child, QtGui.QLineEdit):
                    child.setText(str(value))
                elif isinstance(child, QtGui.QCheckBox):
                    child.setChecked(value)

            self.parent.ui.groupBoxGeometryParameters.setTitle(text)

        else:
            current_index = 0

            self.parent.ui.groupBoxGeometryParameters.setTitle('Parameters')
            self.parent.ui.toolbutton_remove_geometry.setEnabled(False)

        self.parent.animate_stacked_widget(
            self.parent.ui.stackedWidgetGeometryDetails,
            self.parent.ui.stackedWidgetGeometryDetails.currentIndex(),
            current_index,
            'horizontal',
            )

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
            self.parent.get_project_dir(),
            'STL File (*.stl)',))

        if filename:
            name = os.path.basename(filename).lower()
            name = get_unique_string(name, list(self.geometrydict.keys()))

            # reader
            reader = vtk.vtkSTLReader()
            reader.SetFileName(filename)

            # Create transformer
            transform = vtk.vtkTransform()
            transform_filter = vtk.vtkTransformPolyDataFilter()
            transform_filter.SetTransform(transform)
            transform_filter.SetInputConnection(reader.GetOutputPort())

            # mapper
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(transform_filter.GetOutputPort())

            # actor
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetRepresentationToWireframe()

            self.vtkrenderer.AddActor(actor)

            self.vtkrenderer.ResetCamera()
            self.vtkRenderWindow.Render()

            # find center of mass
            center_filter = vtk.vtkCenterOfMass()
            center_filter.SetInputConnection(transform_filter.GetOutputPort())
            center_filter.SetUseScalarsAsWeights(False)
            center_filter.Update()
            center = center_filter.GetCenter()

            # Add to dict
            self.geometrydict[name] = {
                'reader':          reader,
                'transform':       transform,
                'transformfilter': transform_filter,
                'center_filter':   center_filter,
                'mapper':          mapper,
                'actor':           actor,
                'type':            'stl',
                'centerx':         center[0],
                'centery':         center[1],
                'centerz':         center[2],
                'rotationx':       0.0,
                'rotationy':       0.0,
                'rotationz':       0.0,
                'translationx':    0.0,
                'translationy':    0.0,
                'translationz':    0.0,
                }

            # Add to tree
            item = QtGui.QTreeWidgetItem([name])
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(0, QtCore.Qt.Checked)
            self.geometrytree.addTopLevelItem(item)
            self.geometrytree.setCurrentItem(item)

    def parameter_edited(self, widget):

        current_selection = self.geometrytree.selectedItems()

        if current_selection:
            name = str(current_selection[-1].text(0)).lower()
            value = None

            parameters = str(widget.objectName()).lower().split('_')
            parameter = ''.join(parameters[1:])

            for key in self.geometrydict[name].keys():
                if key in parameter:
                    break

            # if widget is a lineedit
            if isinstance(widget, QtGui.QLineEdit):
                string = str(widget.text())
                if 'resolution' in parameter or 'divisions' in parameter:
                    value = int(string)
                else:
                    value = float(string)

            if value is not None:
                self.geometrydict[name][key] = value

                if self.geometrydict[name]['type'] in \
                        list(self.primitivedict.keys()):
                    self.update_primitive(name)
                elif self.geometrydict[name]['type'] in \
                        list(self.filterdict.keys()):
                    self.update_filter(name)
                elif self.geometrydict[name]['type'] in \
                        list(self.parametricdict.keys()):
                    self.update_parametric(name)

                if 'transform' in self.geometrydict[name]:
                    self.update_transform(name)
                self.vtkRenderWindow.Render()

    def update_primitive(self, name):
        primtype = self.geometrydict[name]['type']

        if 'source' in self.geometrydict[name]:
            source = self.geometrydict[name]['source']
        else:
            source = None

        # update source
        if primtype == 'sphere':
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

        elif primtype == 'stl':
            pass

        else:
            return

        # common props
        if source is not None:
            source.SetCenter(self.geometrydict[name]['centerx'],
                             self.geometrydict[name]['centery'],
                             self.geometrydict[name]['centerz'])

            source.Update()

        return source

    def update_transform(self, name):
        transform = self.geometrydict[name]['transform']
        transform_filter = self.geometrydict[name]['transformfilter']

        # reset to Identity
        transform.Identity()
        transform.PostMultiply()

        # translate to center
        transform.Translate(-self.geometrydict[name]['centerx'],
                            -self.geometrydict[name]['centery'],
                            -self.geometrydict[name]['centerz'])

        # rotation
        transform.RotateWXYZ(self.geometrydict[name]['rotationx'], 1, 0, 0)
        transform.RotateWXYZ(self.geometrydict[name]['rotationy'], 0, 1, 0)
        transform.RotateWXYZ(self.geometrydict[name]['rotationz'], 0, 0, 1)

        # back to position
        transform.Translate(self.geometrydict[name]['centerx'],
                            self.geometrydict[name]['centery'],
                            self.geometrydict[name]['centerz'])

        # translate stl files
        if self.geometrydict[name]['type'] in ['stl'] + \
                list(self.parametricdict.keys()):
            transform.Translate(
                self.geometrydict[name]['translationx'],
                self.geometrydict[name]['translationy'],
                self.geometrydict[name]['translationz'],
                )

        # update
        transform_filter.Update()

        return transform_filter

    def add_primitive(self, primtype='sphere'):

        name = get_unique_string(primtype, list(self.geometrydict.keys()))

        # create primitive
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

        # Create transformer
        transform = vtk.vtkTransform()
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputConnection(source.GetOutputPort())

        self.geometrydict[name]['transform'] = transform
        self.geometrydict[name]['transformfilter'] = transform_filter

        # update transform
        self.update_transform(name)

        # convert to triangles
        trianglefilter = vtk.vtkTriangleFilter()
        trianglefilter.SetInputConnection(transform_filter.GetOutputPort())

        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(trianglefilter.GetOutputPort())

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
        self.geometrydict[name]['trianglefilter'] = trianglefilter
        self.geometrydict[name]['source'] = source

        # Add to tree
        item = QtGui.QTreeWidgetItem([name])
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(0, QtCore.Qt.Checked)
        self.geometrytree.addTopLevelItem(item)
        self.geometrytree.setCurrentItem(item)

    def update_parametric(self, name):
        paratype = self.geometrydict[name]['type']
        para_object = self.geometrydict[name]['parametric_object']
        source = self.geometrydict[name]['source']

        if paratype == 'torus':
            para_object.SetRingRadius(self.geometrydict[name]['ringradius'])
            para_object.SetCrossSectionRadius(
                self.geometrydict[name]['crosssectionradius'])
        elif paratype == 'boy':
            para_object.SetZScale(self.geometrydict[name]['zscale'])
        elif paratype == 'conic_spiral':
            para_object.SetA(self.geometrydict[name]['ascale'])
            para_object.SetB(self.geometrydict[name]['bfunc'])
            para_object.SetC(self.geometrydict[name]['cfunc'])
            para_object.SetN(self.geometrydict[name]['nfunc'])
        elif paratype == 'dini':
            para_object.SetA(self.geometrydict[name]['ascale'])
            para_object.SetB(self.geometrydict[name]['bscale'])
        elif paratype == 'ellipsoid':
            para_object.SetXRadius(self.geometrydict[name]['radiusx'])
            para_object.SetYRadius(self.geometrydict[name]['radiusy'])
            para_object.SetZRadius(self.geometrydict[name]['radiusz'])
        elif paratype == 'figure_8_klein':
            para_object.SetRadius(self.geometrydict[name]['radius'])
        elif paratype == 'mobius':
            para_object.SetRadius(self.geometrydict[name]['radius'])

        source.Update()

        return source

    def add_parametric(self, name):

        parametric_object = self.parametricdict[name]()
        source = vtk.vtkParametricFunctionSource()
        source.SetParametricFunction(parametric_object)

        self.geometrydict[name] = {
            'translationx':       0.0,
            'translationy':       0.0,
            'translationz':       0.0,
            'centerx':            0.0,
            'centery':            0.0,
            'centerz':            0.0,
            'rotationx':          0.0,
            'rotationy':          0.0,
            'rotationz':          0.0,
            'radius':             1.0,
            'radiusx':            1.0,
            'radiusy':            1.0,
            'radiusz':            1.0,
            'ringradius':         1.0,
            'crosssectionradius': 0.5,
            'zscale':             0.125,
            'ascale':             0.8,
            'bscale':             0.2,
            'bfunc':              1.0,
            'cfunc':              0.1,
            'nfunc':              2.0,
            'type':               name,
            'parametric_object':  parametric_object,
            'source':             source,
        }

        source = self.update_parametric(name)

        # Create transformer
        transform = vtk.vtkTransform()
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputConnection(source.GetOutputPort())

        self.geometrydict[name]['transform'] = transform
        self.geometrydict[name]['transformfilter'] = transform_filter

        # update transform
        self.update_transform(name)

        # convert to triangles
        trianglefilter = vtk.vtkTriangleFilter()
        trianglefilter.SetInputConnection(transform_filter.GetOutputPort())

        # mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(trianglefilter.GetOutputPort())

        # actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetRepresentationToWireframe()

        self.vtkrenderer.AddActor(actor)

        self.vtkrenderer.ResetCamera()
        self.vtkRenderWindow.Render()

        # add to dict
        self.geometrydict[name]['mapper'] = mapper
        self.geometrydict[name]['actor'] = actor
        self.geometrydict[name]['trianglefilter'] = trianglefilter
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

                if 'trianglefilter' in self.geometrydict[name]:
                    geometry = self.geometrydict[name]['trianglefilter']
                elif 'booleanoperation' in self.geometrydict[name]:
                    geometry = self.geometrydict[name]['booleanoperation']
                elif 'reader' in self.geometrydict[name]:
                    geometry = self.geometrydict[name]['transformfilter']

                # hide the sources
                self.geometrydict[name]['actor'].VisibilityOff()

                boolean_operation.SetInputConnection(i,
                                                     geometry.GetOutputPort())
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

    def copy_geometry(self):
        """ duplicate the selected geometry """
        raise NotImplementedError

#        currentSelection = self.geometrytree.selectedItems()
#        if currentSelection:
#            text = str(currentSelection[-1].text(0)).lower()
#
#            new = get_unique_string(text, self.geometrydict.keys())
#
#            # copy values
#            self.geometrydict[new] = {}
#            for key, value in self.geometrydict[text].items():
#                if not isinstance(value, vtk.vtkObject):
#                    self.geometrydict[new][key] = copy.deepcopy(value)
#
#
#            self.geometrydict[new] = copy.deepcopy(self.geometrydict[text])
#
#            self.vtkrenderer.AddActor(self.geometrydict[new]['actor'])

    def update_filter(self, name):

        filtertype = self.geometrydict[name]['type']
        vtkfilter = self.geometrydict[name]['filter']

        if filtertype == 'clean':
            if self.geometrydict[name]['linestopoints']:
                vtkfilter.ConvertLinesToPointsOn()
            else:
                vtkfilter.ConvertLinesToPointsOff()

            if self.geometrydict[name]['polystolines']:
                vtkfilter.ConvertPolysToLinesOn()
            else:
                vtkfilter.ConvertPolysToLinesOff()

            if self.geometrydict[name]['stripstopolys']:
                vtkfilter.ConvertStripsToPolysOn()
            else:
                vtkfilter.ConvertStripsToPolysOff()

        elif filtertype == 'fill_holes':
            vtkfilter.SetHoleSize(self.geometrydict[name]['maximumholesize'])

        elif filtertype == 'triangle':
            if self.geometrydict[name]['processvertices']:
                vtkfilter.PassVertsOn()
            else:
                vtkfilter.PassVertsOff()

            if self.geometrydict[name]['processlines']:
                vtkfilter.PassLinesOn()
            else:
                vtkfilter.PassLinesOff()

        elif filtertype == 'decimate':
            vtkfilter.SetTargetReduction(
                self.geometrydict[name]['targetreduction'])

        elif filtertype == 'quadric_decimation':
            vtkfilter.SetTargetReduction(
                self.geometrydict[name]['targetreduction'])

        elif filtertype == 'quadric_clustering':
            vtkfilter.SetNumberOfXDivisions(
                self.geometrydict[name]['divisionsx'])
            vtkfilter.SetNumberOfYDivisions(
                self.geometrydict[name]['divisionsy'])
            vtkfilter.SetNumberOfZDivisions(
                self.geometrydict[name]['divisionsz'])

            if self.geometrydict[name]['autoadjustdivisions']:
                vtkfilter.AutoAdjustNumberOfDivisionsOn()
            else:
                vtkfilter.AutoAdjustNumberOfDivisionsOff()

        vtkfilter.Update()

    def add_filter(self, filtertype):

        current_selection = self.geometrytree.selectedItems()
        if current_selection:
            selection_text = str(current_selection[-1].text(0)).lower()

            name = get_unique_string(filtertype,
                                     list(self.geometrydict.keys()))

            vtkfilter = self.filterdict[filtertype]()

            self.geometrydict[name] = {
                'type':   filtertype,
                'filter': vtkfilter,
                'linestopoints': True,
                'polystolines': True,
                'stripstopolys': True,
                'maximumholesize': 1.0,
                'processvertices': True,
                'processlines': True,
                'targetreduction': 0.2,
                'preservetopology': True,
                'splitmesh': True,
                'deletevertices': False,
                'divisionsx': 10,
                'divisionsy': 10,
                'divisionsz': 10,
                'autoadjustdivisions': True,
                }

            if 'trianglefilter' in self.geometrydict[selection_text]:
                inputdata = self.geometrydict[selection_text]['trianglefilter']
            elif 'booleanoperation' in self.geometrydict[selection_text]:
                inputdata = self.geometrydict[selection_text][
                    'booleanoperation']
            elif 'reader' in self.geometrydict[selection_text]:
                inputdata = self.geometrydict[selection_text][
                    'transformfilter']
            elif 'filter' in self.geometrydict[selection_text]:
                inputdata = self.geometrydict[selection_text]['filter']

            # set input data
            vtkfilter.SetInputConnection(inputdata.GetOutputPort())

            # update filter
            self.update_filter(name)

            # hide the source
            self.geometrydict[selection_text]['actor'].VisibilityOff()

            # Create a mapper
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(vtkfilter.GetOutputPort())

            # Create an actor
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetRepresentationToWireframe()

            # add actor to render
            self.vtkrenderer.AddActor(actor)

            # update
            self.vtkrenderer.ResetCamera()
            self.vtkRenderWindow.Render()

            # save references
            self.geometrydict[name]['actor'] = actor
            self.geometrydict[name]['mapper'] = mapper

            # Add to tree
            toplevel = QtGui.QTreeWidgetItem([name])
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

    # --- mesh ---
    def change_mesh_tab(self, tabnum, btn):

        self.parent.animate_stacked_widget(
            self.parent.ui.stackedwidget_mesh,
            self.parent.ui.stackedwidget_mesh.currentIndex(),
            tabnum,
            direction='horizontal',
            line=self.parent.ui.line_mesh,
            to_btn=btn,
            btn_layout=self.parent.ui.gridlayout_mesh_tab_btns,
            )
