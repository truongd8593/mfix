# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

import os
from collections import OrderedDict

# 3rd pary imports
import numpy as np

# Qt imports
from qtpy import QtCore, QtGui

# VTK imports
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

# local imports
from tools.general import get_unique_string, widget_iter, get_icon


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

        self.extent_widgets = [self.parent.ui.lineedit_mesh_min_x,
                               self.parent.ui.lineedit_mesh_max_x,
                               self.parent.ui.lineedit_mesh_min_y,
                               self.parent.ui.lineedit_mesh_max_y,
                               self.parent.ui.lineedit_mesh_min_z,
                               self.parent.ui.lineedit_mesh_max_z,
                               ]
        self.cell_widgets = [self.parent.ui.lineedit_mesh_cells_x,
                             self.parent.ui.lineedit_mesh_cells_y,
                             self.parent.ui.lineedit_mesh_cells_z,
                             ]

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
            ])

        self.filterdict = OrderedDict([
            ('clean',              vtk.vtkCleanPolyData),
            ('fill_holes',         vtk.vtkFillHolesFilter),
            ('triangle',           vtk.vtkTriangleFilter),
            ('decimate',           vtk.vtkDecimatePro),
            ('quadric_decimation', vtk.vtkQuadricDecimation),
            ('quadric_clustering', vtk.vtkQuadricClustering)
        ])

        self.rectilinear_grid = vtk.vtkRectilinearGrid()

        self.grid_viewer_dict = {
            'filters': [],
            'mappers': [],
            'actors':  [],
            }

        # --- layout ---
        self.vlayout = QtGui.QVBoxLayout(self)
        self.vlayout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.vlayout)

        self.button_bar = QtGui.QWidget()
        self.button_bar_layout = QtGui.QHBoxLayout(self.button_bar)
        self.button_bar_layout.setContentsMargins(0, 0, 0, 0)
#        self.button_bar.setLayout(self.button_bar_layout)
        self.vlayout.addWidget(self.button_bar)

        self.vtkWindowWidget = QVTKRenderWindowInteractor(self)
        self.vlayout.addWidget(self.vtkWindowWidget)

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

        # Orientation Arrows Marker Widget
#        self.axes = vtk.vtkAxesActor()
#        self.axes.AxisLabelsOn()
#        self.axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(
#            1, 0, 0)
#        self.axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
#        self.axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(
#            0, 1, 0)
#        self.axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
#        self.axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(
#            0, 0, 1)
#        self.axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()

        # Orientation Cube Marker Widget
        self.axes = vtk.vtkAnnotatedCubeActor();
        self.axes.SetXPlusFaceText('E')
        self.axes.SetXMinusFaceText('W')
        self.axes.SetYMinusFaceText('N')
        self.axes.SetYPlusFaceText('S')
        self.axes.SetZMinusFaceText('T')
        self.axes.SetZPlusFaceText('B')
        self.axes.GetTextEdgesProperty().SetColor(1,1,1)
        self.axes.GetTextEdgesProperty().SetLineWidth(2)
        self.axes.GetCubeProperty().SetColor(.39, .71, .965)

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

        self.__add_buttons()

    # --- setup ---
    def __add_buttons(self):

        self.toolbutton_reset = QtGui.QToolButton()
        self.toolbutton_reset.pressed.connect(self.reset_view)
        self.toolbutton_reset.setIcon(get_icon('overscan.png'))

        self.toolbutton_perspective = QtGui.QToolButton()
        self.toolbutton_perspective.pressed.connect(self.perspective)
        self.toolbutton_perspective.setIcon(get_icon('perspective.png'))

        self.toolbutton_view_xy = QtGui.QToolButton()
        self.toolbutton_view_xy.pressed.connect(lambda: self.set_view('xy'))
        self.toolbutton_view_xy.setIcon(get_icon('-z.png'))

        self.toolbutton_view_yz = QtGui.QToolButton()
        self.toolbutton_view_yz.pressed.connect(lambda: self.set_view('yz'))
        self.toolbutton_view_yz.setIcon(get_icon('-x.png'))

        self.toolbutton_view_xz = QtGui.QToolButton()
        self.toolbutton_view_xz.pressed.connect(lambda: self.set_view('xz'))
        self.toolbutton_view_xz.setIcon(get_icon('-y.png'))

        for btn in [self.toolbutton_reset,
                    self.toolbutton_view_xy,
                    self.toolbutton_view_yz,
                    self.toolbutton_view_xz,
                    self.toolbutton_perspective]:
            self.button_bar_layout.addWidget(btn)
            btn.setAutoRaise(True)

        self.button_bar_layout.addStretch()

    # --- geometry ---
    def tree_widget_geometry_changed(self):
        """
        The selected geometry changed, update UI
        """
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
                if str(widget.objectName()) == self.geometrydict[text]['type']:
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
        """
        Hide/Show the clicked geometry
        """

        name = str(item.text(0)).lower()
        if item.checkState(0) == QtCore.Qt.Unchecked:
            self.geometrydict[name]['actor'].VisibilityOff()
        else:
            self.geometrydict[name]['actor'].VisibilityOn()

        self.vtkRenderWindow.Render()

    def add_stl(self):
        """
        Open browse dialog and load selected stl file
        """

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
        """
        Update the value of edited paraemter in the geometrydict
        """
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
                if any([s in parameter for s in ['divisions', 'resolution',
                                                 'nhills']]):
                    value = int(string)
                else:
                    value = float(string)
            elif isinstance(widget, QtGui.QCheckBox):
                value = widget.isChecked()

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
        """
        Update the specified primitive
        """
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
        """
        Update the specified object's transform filter.
        """
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
        """
        Add the specified primitive
        """

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
        """
        Update the specified parameteric object.
        """
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
        elif paratype == 'random_hills':
            para_object.SetHillXVariance(self.geometrydict[name]['variancex'])
            para_object.SetXVarianceScaleFactor(
                self.geometrydict[name]['scalex'])
            para_object.SetHillYVariance(self.geometrydict[name]['variancey'])
            para_object.SetYVarianceScaleFactor(
                self.geometrydict[name]['scaley'])
            para_object.SetHillAmplitude(self.geometrydict[name]['amplitude'])
            para_object.SetAmplitudeScaleFactor(
                self.geometrydict[name]['scaleamplitude'])
            para_object.SetNumberOfHills(self.geometrydict[name]['nhills'])
            if self.geometrydict[name]['allowrandom']:
                para_object.AllowRandomGenerationOn()
            else:
                para_object.AllowRandomGenerationOff()
        elif paratype == 'roman':
            para_object.SetRadius(self.geometrydict[name]['radius'])
        elif paratype == 'super_ellipsoid':
            para_object.SetXRadius(self.geometrydict[name]['radiusx'])
            para_object.SetYRadius(self.geometrydict[name]['radiusy'])
            para_object.SetZRadius(self.geometrydict[name]['radiusz'])
            para_object.SetN1(self.geometrydict[name]['n1'])
            para_object.SetN2(self.geometrydict[name]['n2'])
        elif paratype == 'super_toroid':
            para_object.SetXRadius(self.geometrydict[name]['radiusx'])
            para_object.SetYRadius(self.geometrydict[name]['radiusy'])
            para_object.SetZRadius(self.geometrydict[name]['radiusz'])
            para_object.SetRingRadius(self.geometrydict[name]['ringradius'])
            para_object.SetCrossSectionRadius(
                self.geometrydict[name]['crosssectionradius'])
            para_object.SetN1(self.geometrydict[name]['n1'])
            para_object.SetN2(self.geometrydict[name]['n2'])

        source.Update()

        return source

    def add_parametric(self, name):
        """
        Add the specified parametric object
        """

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
            'nhills':             30,
            'variancex':          2.5,
            'scalex':             0.3,
            'variancey':          2.5,
            'scaley':             0.3,
            'amplitude':          2.0,
            'scaleamplitude':     0.3,
            'allowrandom':        True,
            'n1':                 1.0,
            'n2':                 1.0,
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
        """
        Apply a boolean operation with the currently selected toplevel items.
        """

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
        """
        Remove the currently selected geometry, filter, or boolean operation
        """
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
        """
        Update the currently selected filter
        """

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
        """
        add the selected filter with the input being the currently selected
        toplevel item
        """

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

            inputdata = self.get_input_data(selection_text)

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

    def get_input_data(self, name):
        """ based on the type of geomtry, return the data """

        if 'trianglefilter' in self.geometrydict[name]:
            inputdata = self.geometrydict[name]['trianglefilter']
        elif 'booleanoperation' in self.geometrydict[name]:
            inputdata = self.geometrydict[name][
                'booleanoperation']
        elif 'reader' in self.geometrydict[name]:
            inputdata = self.geometrydict[name][
                'transformfilter']
        elif 'filter' in self.geometrydict[name]:
            inputdata = self.geometrydict[name]['filter']

        return inputdata

    def collect_toplevel_geometry(self):
        """ collect and append visible toplevel polydata """

        append_filter = vtk.vtkAppendPolyData()
        for top_num in range(self.geometrytree.topLevelItemCount()):
            item = self.geometrytree.topLevelItem(top_num)
            if item.checkState(0) == QtCore.Qt.Checked:
                append_filter.AddInputData(
                    self.get_input_data(str(item.text(0))).GetOutput())
        append_filter.Update()

        return append_filter

    def get_geometry_extents(self):
        """ determine the extents of the visible geometry """

        geometry = self.collect_toplevel_geometry()

        return geometry.GetOutput().GetBounds()

    def export_stl(self, file_name):
        """ expoort visivle toplevel geometry """

        geometry = self.collect_toplevel_geometry()

        # clean
        clean_filter = vtk.vtkCleanPolyData()
        clean_filter.SetInputConnection(geometry.GetOutputPort())
        clean_filter.Update()

        # write file
        stl_writer = vtk.vtkSTLWriter()
        stl_writer.SetFileName(file_name)
        stl_writer.SetInputConnection(clean_filter.GetOutputPort())
        stl_writer.Write()

    # --- mesh ---
    def change_mesh_tab(self, tabnum, btn):
        """ switch mesh stacked widget based on selected """
        self.parent.animate_stacked_widget(
            self.parent.ui.stackedwidget_mesh,
            self.parent.ui.stackedwidget_mesh.currentIndex(),
            tabnum,
            direction='horizontal',
            line=self.parent.ui.line_mesh,
            to_btn=btn,
            btn_layout=self.parent.ui.gridlayout_mesh_tab_btns,
            )

    def auto_size_mesh_extents(self):
        """ collect and set the extents of the visible geometry """
        extents = self.get_geometry_extents()

        for widget, extent in zip(self.extent_widgets, extents):
            widget.setText(str(extent))

        self.update_mesh()

    def update_mesh(self):

        extents = []
        for widget in self.extent_widgets:
            extents.append(float(widget.text()))

        cells = []
        for widget in self.cell_widgets:
            cells.append(int(widget.text()))

        self.rectilinear_grid.SetDimensions(*cells)

        # determine cell spacing
        x_coords = vtk.vtkFloatArray()
        for i in np.linspace(extents[0], extents[1], cells[0]):
            x_coords.InsertNextValue(i)

        y_coords = vtk.vtkFloatArray()
        for i in np.linspace(extents[2], extents[3], cells[1]):
            y_coords.InsertNextValue(i)

        z_coords = vtk.vtkFloatArray()
        for i in np.linspace(extents[4], extents[5], cells[2]):
            z_coords.InsertNextValue(i)

        self.rectilinear_grid.SetXCoordinates(x_coords)
        self.rectilinear_grid.SetYCoordinates(y_coords)
        self.rectilinear_grid.SetZCoordinates(z_coords)

        # update actors
        # remove exsisting
        for actor in self.grid_viewer_dict['actors']:
            self.vtkrenderer.RemoveActor(actor)

        self.grid_viewer_dict['filters'] = []
        self.grid_viewer_dict['actors'] = []
        self.grid_viewer_dict['mappers'] = []

        # add new actors
        for i in range(3):
            filter_ = vtk.vtkRectilinearGridGeometryFilter()
            filter_.SetInputData(self.rectilinear_grid)

            if i == 0:
                filter_.SetExtent(0, 0, 0, cells[1], 0, cells[2])
            elif i == 1:
                filter_.SetExtent(0, cells[0], 0, 0, 0, cells[2])
            else:
                filter_.SetExtent(0, cells[0], 0, cells[1], 0, 0)
            filter_.Update()

            self.grid_viewer_dict['filters'].append(filter_)
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(filter_.GetOutputPort())
            self.grid_viewer_dict['mappers'].append(mapper)

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetRepresentationToWireframe()
            actor.GetProperty().SetColor(.39, .71, .965)
            self.grid_viewer_dict['actors'].append(actor)

            self.vtkrenderer.AddActor(actor)

        self.vtkrenderer.ResetCamera()
        self.vtkRenderWindow.Render()

    def vtk_mesher(self):

        signedDistances = vtk.vtkFloatArray()
        signedDistances.SetNumberOfComponents(1)
        signedDistances.SetName("SignedDistances")

        implicitPolyDataDistance = vtk.vtkImplicitPolyDataDistance()
        source = self.collect_toplevel_geometry()
        implicitPolyDataDistance.SetInput(source.GetOutput())

        # Evaluate the signed distance function at all of the grid points
        for point_id in range(self.rectilinear_grid.GetNumberOfPoints()):
            p = self.rectilinear_grid.GetPoint(point_id)
            signedDistance = implicitPolyDataDistance.EvaluateFunction(p)
            signedDistances.InsertNextValue(signedDistance)

        self.rectilinear_grid.GetPointData().SetScalars(signedDistances)

        clipper = vtk.vtkClipDataSet()
        clipper.SetInputData(self.rectilinear_grid)
        clipper.InsideOutOn()
        clipper.SetValue(0.0)
        clipper.Update()

        clipperMapper = vtk.vtkDataSetMapper()
        clipperMapper.SetInputConnection(clipper.GetOutputPort())

        clipperActor = vtk.vtkActor()
        clipperActor.SetMapper(clipperMapper)
        clipperActor.GetProperty().SetRepresentationToWireframe()
        clipperActor.GetProperty().SetColor(.1, .1, .1)

        self.vtkrenderer.AddActor(clipperActor)
        self.vtkRenderWindow.Render()

    def mesh(self):

        mesher = str(self.parent.ui.combobox_mesher.currentText())

        if mesher == 'vtkClipDataSet':
            self.vtk_mesher()

    # --- view ---
    def perspective(self):
        camera = self.vtkrenderer.GetActiveCamera()

        if camera.GetParallelProjection():
            camera.ParallelProjectionOff()
            self.toolbutton_perspective.setIcon(get_icon('perspective.png'))
        else:
            camera.ParallelProjectionOn()
            self.toolbutton_perspective.setIcon(get_icon('parallel.png'))

        self.vtkRenderWindow.Render()

    def set_view(self, view='xy'):

        camera = self.vtkrenderer.GetActiveCamera()

        if view == 'xy':
            camera.SetPosition(0, 0, 10000000)
            camera.SetViewUp(0, 1, 0)
        elif view == 'yz':
            camera.SetPosition(0, 10000000, 0)
            camera.SetViewUp(1, 0, 0)
        elif view == 'xz':
            camera.SetPosition(10000000, 0, 0)
            camera.SetViewUp(0, 0, 1)
        self.vtkrenderer.ResetCamera()
        self.vtkRenderWindow.Render()

    def reset_view(self):
        self.vtkrenderer.ResetCamera()
        self.vtkRenderWindow.Render()
