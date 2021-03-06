# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

import os
import sys
from collections import OrderedDict

# 3rd pary imports
import numpy as np

# Qt imports
from qtpy import QtCore, QtGui

# VTK imports
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

# local imports
from tools.general import (get_unique_string, widget_iter, get_icon,
                           get_image_path, make_callback,)

CELL_TYPE_ENUM = {
    0:  'empty_cell',
    1:  'vertex',
    2:  'poly_vertex',
    3:  'line',
    4:  'poly_line',
    5:  'triangle',
    6:  'triangle_strip',
    7:  'polygon',
    8:  'pixel',
    9:  'quad',
    10: 'tetra',
    11: 'voxel',
    12: 'hexahedron',
    13: 'wedge',
    14: 'pyramid',
    15: 'pentagonal_prism',
    16: 'hexagonal_prism',
    21: 'quadratic_edge',
    22: 'quadratic_triangle',
    23: 'quadratic_quad',
    24: 'quadratic_tetra',
    25: 'quadratic_hexahedron',
    26: 'quadratic_wedge',
    27: 'quadratic_pyramid',
    28: 'biquadratic_quad',
    29: 'triquadratic_hexahedron',
    30: 'quadratic_linear_quad',
    31: 'quadratic_linear_wedge',
    32: 'biquadratic_quadratic_wedge',
    33: 'biquadratic_quadratic_hexahedron',
    34: 'biquadratic_traingle',
    35: 'cubic_line',
    36: 'quadratic_polygon',
    }


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


class CustomOrientationMarkerWidget(vtk.vtkOrientationMarkerWidget):

    def __init__(self, parent=None):
        self.RemoveAllObservers()
        self.AddObserver("OnLeftButtonPressEvent",
                         self.left_button_press_event)

    def left_button_press_event(self, obj, event):
        print('pressed')


class VtkWidget(QtGui.QWidget):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, project, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.project = project
        self.parent = parent
        self.geometrytree = self.parent.ui.treeWidgetGeometry
        self.booleanbtndict = self.parent.booleanbtndict

        # --- data ---
        self.geometrydict = {}
        self.geometry_visible = True

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

        # set default colors
        self.color_dict = {
            'mesh':                QtGui.QColor(244,  67,  54),
            'background_mesh':     QtGui.QColor(100, 182, 247),
            'geometry':            QtGui.QColor(224, 224, 224),
            'regions':             QtGui.QColor(224, 224, 224),
            'initial_conditions':  QtGui.QColor(224, 224, 224),
            'boundary_conditions': QtGui.QColor(224, 224, 224),
            'point_sources':       QtGui.QColor(224, 224, 224),
            }

        # add edge color
        for color in list(self.color_dict.keys()):
            self.color_dict['_'.join([color, 'edge'])] = \
                self.color_dict[color].darker()

        # --- layout ---
        self.grid_layout = QtGui.QGridLayout(self)
        self.grid_layout.setContentsMargins(0, 0, 0, 0)

        self.button_bar = QtGui.QWidget(self)
        self.button_bar_layout = QtGui.QHBoxLayout(self.button_bar)
        self.button_bar_layout.setContentsMargins(0, 0, 0, 0)
        self.button_bar.setLayout(self.button_bar_layout)
        self.button_bar.setGeometry(QtCore.QRect(0, 0, 300, 300))
        self.grid_layout.addWidget(self.button_bar, 0, 0)

        self.vtkWindowWidget = QVTKRenderWindowInteractor(self)
        self.grid_layout.addWidget(self.vtkWindowWidget, 1, 0)

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
        self.axes = vtk.vtkAnnotatedCubeActor()
        self.axes.SetXPlusFaceText('E')
        self.axes.SetXMinusFaceText('W')
        self.axes.SetYMinusFaceText('N')
        self.axes.SetYPlusFaceText('S')
        self.axes.SetZMinusFaceText('T')
        self.axes.SetZPlusFaceText('B')
        self.axes.GetTextEdgesProperty().SetColor(1, 1, 1)
        self.axes.GetTextEdgesProperty().SetLineWidth(2)
        self.axes.GetCubeProperty().SetColor(.39, .71, .965)
        self.axes.PickableOn()

        self.orientation_widget = CustomOrientationMarkerWidget()
        self.orientation_widget.SetOutlineColor(0.9300, 0.5700, 0.1300)
        self.orientation_widget.SetOrientationMarker(self.axes)
        self.orientation_widget.SetInteractor(self.vtkiren)
        self.orientation_widget.SetViewport(0.0, 0.0, 0.2, 0.2)
        self.orientation_widget.SetEnabled(1)
        self.orientation_widget.InteractiveOff()
        self.orientation_widget.PickingManagedOn()

        self.vtkrenderer.ResetCamera()

        # --- setup vtk mappers/actors ---
        self.mesh = None
        self.mesh_mapper = vtk.vtkDataSetMapper()
        self.mesh_mapper.ScalarVisibilityOff()

        self.mesh_actor = vtk.vtkActor()
        self.mesh_actor.SetMapper(self.mesh_mapper)
        self.mesh_actor.GetProperty().SetRepresentationToWireframe()
        self.mesh_actor.GetProperty().SetColor(
            self.color_dict['mesh'].getRgbF()[:3])
        self.mesh_actor.GetProperty().SetEdgeColor(
            self.color_dict['mesh_edge'].getRgbF()[:3])

        self.vtkrenderer.AddActor(self.mesh_actor)

        self.__connect_events()
        self.__add_tool_buttons()

    # --- setup ---
    def __connect_events(self):

        # -- Geometry Tree ---
        self.geometrytree.setStyleSheet(
            "QTreeView::indicator:unchecked {image: url(%s);}"
            "QTreeView::indicator:checked {image: url(%s);}"
            % (get_image_path('visibilityofftransparent.png'),
               get_image_path('visibility.png'))
            )
        self.geometrytree.itemSelectionChanged.connect(
            self.tree_widget_geometry_changed)
        self.geometrytree.itemClicked.connect(self.geometry_clicked)

        # --- geometry button ---
        self.add_geometry_menu = QtGui.QMenu(self)
        self.parent.ui.toolbutton_add_geometry.setMenu(self.add_geometry_menu)

        action = QtGui.QAction('STL File',  self.add_geometry_menu)
        action.triggered.connect(self.add_stl)
        self.add_geometry_menu.addAction(action)

        self.add_geometry_menu.addSeparator()

        for geo in self.primitivedict.keys():
            action = QtGui.QAction(geo, self.add_geometry_menu)
            action.triggered.connect(
                make_callback(self.add_primitive, geo))
            self.add_geometry_menu.addAction(action)

        self.add_geometry_menu.addSeparator()

        for geo in self.parametricdict.keys():
            action = QtGui.QAction(geo.replace('_', ' '),
                                   self.add_geometry_menu)
            action.triggered.connect(
                make_callback(self.add_parametric, geo))
            self.add_geometry_menu.addAction(action)

        # --- filter button ---
        self.add_filter_menu = QtGui.QMenu(self)
        self.parent.ui.toolbutton_add_filter.setMenu(self.add_filter_menu)

        for geo in self.filterdict.keys():
            action = QtGui.QAction(geo.replace('_', ' '),
                                   self.add_filter_menu)
            action.triggered.connect(
                make_callback(self.add_filter, geo))
            self.add_filter_menu.addAction(action)

        # setup signals
        self.parent.ui.toolbutton_remove_geometry.pressed.connect(
            self.remove_geometry)
        self.parent.ui.toolbutton_copy_geometry.pressed.connect(
            self.copy_geometry)

        # connect boolean
        for key, btn in self.booleanbtndict.items():
            btn.pressed.connect(
                make_callback(self.boolean_operation, key))

        # connect parameter widgets
        for widget in widget_iter(self.parent.ui.stackedWidgetGeometryDetails):
            if isinstance(widget, QtGui.QLineEdit):
                widget.editingFinished.connect(
                    make_callback(self.parameter_edited, widget))
            elif isinstance(widget, QtGui.QCheckBox):
                widget.stateChanged.connect(
                    make_callback(self.parameter_edited, widget))

        # --- mesh ---
        # connect mesh tab btns
        for i, btn in enumerate([self.parent.ui.pushbutton_mesh_uniform,
                                 self.parent.ui.pushbutton_mesh_controlpoints,
                                 self.parent.ui.pushbutton_mesh_mesher]):
            btn.pressed.connect(
                make_callback(self.change_mesh_tab, i, btn))

        self.parent.ui.pushbutton_mesh_autosize.pressed.connect(
            self.auto_size_mesh_extents)

        # connect mesher
        self.parent.ui.combobox_mesher.currentIndexChanged.connect(
                    self.change_mesher_options)

        self.parent.ui.pushbutton_generate_mesh.pressed.connect(self.mesher)

    def __add_tool_buttons(self):

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

        self.toolbutton_screenshot = QtGui.QToolButton()
        self.toolbutton_screenshot.pressed.connect(self.screenshot)
        self.toolbutton_screenshot.setIcon(get_icon('camera.png'))

        self.toolbutton_visible = QtGui.QToolButton()
        self.toolbutton_visible.setIcon(get_icon('visibility.png'))

        self.visible_menu = QtGui.QMenu(self)
        self.toolbutton_visible.setMenu(self.visible_menu)
        self.toolbutton_visible.setPopupMode(QtGui.QToolButton.InstantPopup)

        # --- visual representation menu ---
        layout = QtGui.QGridLayout(self.visible_menu)
        layout.setContentsMargins(5, 5, 5, 5)
        for i, geo in enumerate(['Background Mesh', 'Mesh', 'Geometry',
                                 'Regions', 'Initial Conditions',
                                 'Boundary Conditions', 'Point Sources']):

            # tool button
            toolbutton = QtGui.QToolButton()
            toolbutton.pressed.connect(make_callback(self.change_visibility,
                                                     geo.lower(), toolbutton))
            toolbutton.setCheckable(True)
            toolbutton.setChecked(True)
            toolbutton.setAutoRaise(True)
            toolbutton.setIcon(get_icon('visibility.png'))
            layout.addWidget(toolbutton, i, 0)

            # style
            combobox = QtGui.QComboBox()
            combobox.addItems(['wire', 'solid', 'solid with edges', 'points'])
            combobox.currentIndexChanged.connect(
                make_callback(self.change_representation,
                              geo.lower(), combobox))
            layout.addWidget(combobox, i, 1)

            # color
            toolbutton = QtGui.QToolButton()
            toolbutton.pressed.connect(make_callback(self.change_color,
                                                     geo.lower(), toolbutton))
            toolbutton.setAutoRaise(True)
            toolbutton.setStyleSheet("QToolButton{{ background: {};}}".format(
                self.color_dict[
                    '_'.join(geo.lower().split(' '))].name())
                )
            layout.addWidget(toolbutton, i, 2)

            # label
            label = QtGui.QLabel(geo)
            layout.addWidget(label, i, 3)

        for btn in [self.toolbutton_reset,
                    self.toolbutton_view_xy,
                    self.toolbutton_view_yz,
                    self.toolbutton_view_xz,
                    self.toolbutton_perspective,
                    self.toolbutton_screenshot,
                    self.toolbutton_visible]:
            self.button_bar_layout.addWidget(btn)
            btn.setAutoRaise(True)

        self.button_bar_layout.addStretch()

    def emitUpdatedValue(self, key, value, args=None):
        self.value_updated.emit(self, {key: value}, args)

    def updateValue(self, key, newValue, args=None):

        if key in ['xmin', 'xlength', 'ymin', 'ylength', 'zmin', 'zlength',
                   'imax', 'jmax', 'kmax']:
            self.update_mesh()

    def objectName(self):
        return 'VTK Widget'

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
        if item.checkState(0) == QtCore.Qt.Checked:
            if self.geometry_visible:
                self.geometrydict[name]['actor'].VisibilityOn()
            self.geometrydict[name]['visible'] = True
        else:
            self.geometrydict[name]['actor'].VisibilityOff()
            self.geometrydict[name]['visible'] = False

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
            actor.GetProperty().SetColor(
                self.color_dict['geometry'].getRgbF()[:3])
            actor.GetProperty().SetEdgeColor(
                self.color_dict['geometry_edge'].getRgbF()[:3])

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
                'visible':         True,
                }

            # Add to tree
            item = QtGui.QTreeWidgetItem([name])
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(0, QtCore.Qt.Checked)
            self.geometrytree.addTopLevelItem(item)
            self.geometrytree.setCurrentItem(item)

    def parameter_edited(self, widget):
        """
        Update the value of edited parameter in the geometrydict
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
            'visible':         True,
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

        # copy properties from an exsiting actor
        if len(self.geometrydict) > 1:
            other_actor = list(self.geometrydict.keys())
            other_actor.remove(name)
            other_actor = self.geometrydict[other_actor[0]]['actor']
            actor.GetProperty().DeepCopy(other_actor.GetProperty())
        else:
            actor.GetProperty().SetRepresentationToWireframe()

        # make sure the colors are correct
        actor.GetProperty().SetColor(
            self.color_dict['geometry'].getRgbF()[:3])
        actor.GetProperty().SetEdgeColor(
            self.color_dict['geometry_edge'].getRgbF()[:3])

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
            'visible':            True
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
        actor.GetProperty().SetColor(
            self.color_dict['geometry'].getRgbF()[:3])
        actor.GetProperty().SetEdgeColor(
                self.color_dict['geometry_edge'].getRgbF()[:3])

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
            actor.GetProperty().SetColor(
                self.color_dict['geometry'].getRgbF()[:3])
            actor.GetProperty().SetEdgeColor(
                self.color_dict['geometry_edge'].getRgbF()[:3])

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
                'type':                filtertype,
                'filter':              vtkfilter,
                'linestopoints':       True,
                'polystolines':        True,
                'stripstopolys':       True,
                'maximumholesize':     1.0,
                'processvertices':     True,
                'processlines':        True,
                'targetreduction':     0.2,
                'preservetopology':    True,
                'splitmesh':           True,
                'deletevertices':      False,
                'divisionsx':          10,
                'divisionsy':          10,
                'divisionsz':          10,
                'autoadjustdivisions': True,
                'visible':             True,
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
            actor.GetProperty().SetColor(
                self.color_dict['geometry'].getRgbF()[:3])
            actor.GetProperty().SetEdgeColor(
                self.color_dict['geometry_edge'].getRgbF()[:3])

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
        item_count = self.geometrytree.topLevelItemCount()

        for top_num in range(item_count):
            item = self.geometrytree.topLevelItem(top_num)
            if item.checkState(0) == QtCore.Qt.Checked:
                append_filter.AddInputData(
                    self.get_input_data(str(item.text(0))).GetOutput())

        # check to make sure there is geometry
        if append_filter.GetTotalNumberOfInputConnections() <= 0:
            self.parent.message(title='Warning',
                                icon='warning',
                                text='There is no visible geometry. Aborting',
                                buttons=['ok'],
                                default='ok',
                                infoText=None,
                                detailedtext=None,
                                )
            raise ValueError('There is no visible geometry. Aborting')

        append_filter.Update()

        # clean
        clean_filter = vtk.vtkCleanPolyData()
        clean_filter.SetInputConnection(append_filter.GetOutputPort())
        clean_filter.Update()

        return clean_filter

    def get_geometry_extents(self):
        """ determine the extents of the visible geometry """

        geometry = self.collect_toplevel_geometry()

        return geometry.GetOutput().GetBounds()

    # --- output files ---
    def export_stl(self, file_name):
        """ expoort visivle toplevel geometry """

        geometry = self.collect_toplevel_geometry()

        # write file
        stl_writer = vtk.vtkSTLWriter()
        stl_writer.SetFileName(file_name)
        stl_writer.SetInputConnection(geometry.GetOutputPort())
        stl_writer.Write()

    def export_unstructured(self, fname, grid):
        """ export an unstructured grid """
        gw = vtk.vtkXMLUnstructuredGridWriter()
        gw.SetFileName(fname)
        gw.SetInputConnection(grid)
        gw.Write()

    def screenshot(self, fname=None):

        self.toolbutton_screenshot.setDown(False)

        if fname is None:
            fname = str(QtGui.QFileDialog.getSaveFileName(
                self.parent,
                "Save screenshot",
                self.parent.settings.value('project_dir'),
                ';;'.join(["PNG (*.png)",
                           "JPEG (*.jpg)",
                           "PostScript (*.ps)",
                           "All Files (*.*)",
                           ]),
                ))

        # screenshot code:
        window_image = vtk.vtkWindowToImageFilter()
        window_image.SetInput(self.vtkRenderWindow)
        window_image.Update()

        if fname.endswith('.png'):
            writer = vtk.vtkPNGWriter()
        elif fname.endswith('.jpg'):
            writer = vtk.vtkJPEGWriter()
        elif fname.endswith('.ps'):
            writer = vtk.vtkPostScriptWriter()
        else:
            raise TypeError('No available writer')

        writer.SetFileName(fname)
        writer.SetInputConnection(window_image.GetOutputPort())
        writer.Write()

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

    def change_mesher_options(self):
        """ switch the mesh options stacked widget """

        mesher = str(self.parent.ui.combobox_mesher.currentText()).lower()

        current_index = 0
        for i in range(self.parent.ui.stackedwidget_mesher_options.count()):
            widget = self.parent.ui.stackedwidget_mesher_options.widget(i)
            if mesher == str(widget.objectName()).lower():
                current_index = i
                break

        self.parent.animate_stacked_widget(
            self.parent.ui.stackedwidget_mesher_options,
            self.parent.ui.stackedwidget_mesher_options.currentIndex(),
            current_index,
            direction='horizontal',
            )

    def auto_size_mesh_extents(self):
        """ collect and set the extents of the visible geometry """
        extents = self.get_geometry_extents()

        for key, extent in zip(['xmin', 'xlength', 'ymin', 'ylength',
                                'zmin', 'zlength'],
                               extents):
            self.emitUpdatedValue(key, extent)

        self.update_mesh()

    def update_mesh(self):

        extents = []
        for key in ['xmin', 'xlength', 'ymin', 'ylength', 'zmin', 'zlength']:
            if key in self.project:
                extents.append(float(self.project[key]))
            else:
                extents.append(0.0)

        cells = []
        for key in ['imax', 'jmax', 'kmax']:
            if key in self.project:
                cells.append(int(self.project[key])+1)
            else:
                cells.append(1)

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

        # copy properties of one of the actors
        mesh_is_visible = 1
        actor_property = None
        if self.grid_viewer_dict['actors']:
            actor_property = vtk.vtkProperty()
            actor_property.DeepCopy(
                self.grid_viewer_dict['actors'][0].GetProperty())

            # reset the colors
            actor_property.SetColor(
                    self.color_dict['background_mesh'].getRgbF()[:3])
            actor_property.SetEdgeColor(
                 self.color_dict['background_mesh_edge'].getRgbF()[:3])
                 
            mesh_is_visible = self.grid_viewer_dict['actors'][0].GetVisibility()

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
            mapper.ScalarVisibilityOff()
            self.grid_viewer_dict['mappers'].append(mapper)

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.SetVisibility(mesh_is_visible)
            if actor_property is not None:
                actor.GetProperty().DeepCopy(actor_property)
            else:
                actor.GetProperty().SetRepresentationToWireframe()
                actor.GetProperty().SetColor(
                    self.color_dict['background_mesh'].getRgbF()[:3])
                actor.GetProperty().SetEdgeColor(
                 self.color_dict['background_mesh_edge'].getRgbF()[:3])
            self.grid_viewer_dict['actors'].append(actor)

            self.vtkrenderer.AddActor(actor)

        self.vtkrenderer.ResetCamera()
        self.vtkRenderWindow.Render()

    def vtk_calc_distance_from_geometry(self):

        signed_distances = vtk.vtkFloatArray()
        signed_distances.SetNumberOfComponents(1)
        signed_distances.SetName("SignedDistances")

        implicit_poly_data_distance = vtk.vtkImplicitPolyDataDistance()
        source = self.collect_toplevel_geometry()
        implicit_poly_data_distance.SetInput(source.GetOutput())

        # Evaluate the signed distance function at all of the grid points
        for point_id in range(self.rectilinear_grid.GetNumberOfPoints()):
            p = self.rectilinear_grid.GetPoint(point_id)
            signed_distance = implicit_poly_data_distance.EvaluateFunction(p)
            signed_distances.InsertNextValue(signed_distance)

        self.rectilinear_grid.GetPointData().SetScalars(signed_distances)

    def vtk_mesher_table_based(self):

        self.vtk_calc_distance_from_geometry()

        clipper = vtk.vtkTableBasedClipDataSet()
        clipper.SetInputData(self.rectilinear_grid)
        if self.parent.ui.checkbox_mesh_inside.isChecked():
            clipper.InsideOutOn()
        else:
            clipper.InsideOutOff()

        clipper.SetMergeTolerance(
            float(self.parent.ui.lineedit_vtk_mesh_merge.text())
            )
        clipper.SetValue(0.0)
        clipper.Update()
        self.mesh = clipper.GetOutput()

        self.mesh_mapper.SetInputData(self.mesh)

        self.vtkRenderWindow.Render()

        self.mesh_stats()

        # export geometry
        project_dir = self.parent.settings.value('project_dir')
        self.export_unstructured(os.path.join(project_dir, 'mesh.vtu'),
                                 clipper.GetOutputPort())

    def vtk_mesher_threshold(self):

        self.vtk_calc_distance_from_geometry()

        thresh = vtk.vtkThreshold()
        thresh.SetInputData(self.rectilinear_grid)

        if self.parent.ui.checkbox_threshold_inside.isChecked():
            thresh.ThresholdByLower(0)
        else:
            thresh.ThresholdByUpper(0)

        if self.parent.ui.checkbox_threshold_interface.isChecked():
            thresh.AllScalarsOff()
        else:
            thresh.AllScalarsOn()

        thresh.Update()

        self.mesh = thresh.GetOutput()

        self.mesh_mapper.SetInputData(self.mesh)

        self.vtkRenderWindow.Render()

        self.mesh_stats()

        # export geometry
        project_dir = self.parent.settings.value('project_dir')
        self.export_unstructured(os.path.join(project_dir, 'mesh.vtu'),
                                 thresh.GetOutputPort())

    def mesh_stats(self):

        cell_count = self.mesh.GetNumberOfCells()
        cell_types = self.mesh.GetCellTypesArray()

        cell_list = []
        for i in range(cell_types.GetNumberOfTuples()):
            cell_list.append(cell_types.GetTuple(i))

        cell_type_counts = []
        for cell in set(cell_list):
            cell_type_counts.append([
                CELL_TYPE_ENUM[int(cell[0])],
                cell_list.count(cell)
                ])

        self.parent.ui.lineedit_mesh_cells.setText(str(cell_count))

        self.parent.ui.plaintextedit_mesh_cell_types.setPlainText(
            '\n'.join([':\t'.join([str(cell_type), str(cells)])
                      for cell_type, cells in cell_type_counts])
            )

    def mesher(self):

        mesher = str(self.parent.ui.combobox_mesher.currentText())

        if mesher == 'vtkTableBasedClipDataSet':
            self.vtk_mesher_table_based()
        elif mesher == 'vtkThreshold':
            self.vtk_mesher_threshold()

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

    def change_visibility(self, name, toolbutton):

        actors = None
        if name == 'mesh':
            actors = [self.mesh_actor]
        elif name == 'background mesh':
            actors = self.grid_viewer_dict['actors']
        elif name == 'geometry':
            actors = [geo['actor'] for geo in self.geometrydict.values()
                      if geo['visible']]

            if toolbutton.isChecked():
                self.geometry_visible = False
            else:
                self.geometry_visible = True

        if actors is not None:
            for actor in actors:
                if toolbutton.isChecked():
                    actor.VisibilityOff()
                    toolbutton.setIcon(
                        get_icon('visibilityofftransparent.png'))
                else:
                    actor.VisibilityOn()
                    toolbutton.setIcon(get_icon('visibility.png'))

            self.vtkRenderWindow.Render()

    def change_representation(self, name, combobox):

        representation = str(combobox.currentText())
        actors = None
        if name == 'mesh':
            actors = [self.mesh_actor]
        elif name == 'background mesh':
            actors = self.grid_viewer_dict['actors']
        elif name == 'geometry':
            actors = [geo['actor'] for geo in self.geometrydict.values()]

        if actors is not None:
            for actor in actors:
                if representation == 'wire':
                    actor.GetProperty().SetRepresentationToWireframe()
                elif representation == 'solid':
                    actor.GetProperty().SetRepresentationToSurface()
                    actor.GetProperty().EdgeVisibilityOff()
                elif representation == 'solid with edges':
                    actor.GetProperty().SetRepresentationToSurface()
                    actor.GetProperty().EdgeVisibilityOn()
                elif representation == 'points':
                    actor.GetProperty().SetRepresentationToPoints()

            self.vtkRenderWindow.Render()

    def change_color(self, name, button):

        col = QtGui.QColorDialog.getColor()

        if col.isValid():

            button.setStyleSheet("QToolButton{{ background: {};}}".format(
                col.name()))

            name = '_'.join(name.lower().split(' '))
            name_edge = '_'.join([name, 'edge'])

            self.color_dict[name] = col
            self.color_dict[name_edge] = col.darker()

            actors = None
            if name == 'mesh':
                actors = [self.mesh_actor]
            elif name == 'background_mesh':
                actors = self.grid_viewer_dict['actors']
            elif name == 'geometry':
                actors = [geo['actor'] for geo in self.geometrydict.values()]

            if actors is not None:
                for actor in actors:
                    actor.GetProperty().SetColor(
                        self.color_dict[name].getRgbF()[:3])
                    actor.GetProperty().SetEdgeColor(
                     self.color_dict[name_edge].getRgbF()[:3])

                self.vtkRenderWindow.Render()
