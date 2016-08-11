# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division
import os
import copy
import logging
from functools import partial
import numpy as np
from qtpy import QtCore, QtGui, QtWidgets
import math

# VTK imports
import vtk
try:
    # Try Qt 5.x
    from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
except ImportError:
    # Fall back to Qt 4.x
    from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

# local imports
from tools.general import (get_unique_string, widget_iter, get_icon,
                           get_image_path, topological_sort)
from widgets.base import LineEdit
from project import Equation, ExtendedJSON
from widgets.vtk_constants import *

log = logging.getLogger(__name__)

gui = None


def safe_float(x):
    try:
        return float(x)
    except ValueError as e:
        if gui:
            gui.error(str(e))
        else:
            log.error(str(e))
        return 0.0


def clean_geo_dict(d):
    clean_dict = {}

    for geo, geo_dict in d.items():
        clean_dict[geo] = {}
        geo_type = None
        if 'geo_type' in geo_dict:
            geo_type = geo_dict['geo_type']
        clean_dict[geo]['geo_type'] = geo_type
        for key, value in geo_dict.items():
            # filter out vtk objects/default values
            if key not in ['mapper', 'actor', 'reader', 'transform',
                           'transformfilter', 'center_filter', 'source',
                           'trianglefilter', 'parametric_object', 'filter',
                           'booleanoperation']:
                if isinstance(value, Equation) or value != DEFAULT_PARAMS[geo_type][key]:
                    clean_dict[geo][key] = value
    return clean_dict

def clean_visual_dict(d):
    clean_dict = {}
    for geo, geo_dict in d.items():
        clean_dict[geo] = {}
        for key, value in geo_dict.items():
            if key in ['color', 'edge']:
                clean_dict[geo][key] = d[geo][key].getRgb()
            else:
                clean_dict[geo][key] = d[geo][key]

    return clean_dict


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


class VtkWidget(QtWidgets.QWidget):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, project, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        global gui
        gui = parent

        self.project = project
        self.parent = parent

        self.ui = parent.ui
        self.geometrytree = self.ui.geometry.treeWidgetGeometry

        # --- data ---
        self.geometrydict = {}
        self.defer_render = False
        self.region_dict = {}
        self.parameter_key_map = {}
        self.view_flip = [False]*3
        self.booleanbtndict = {
            'union':        self.ui.geometry.toolbutton_geometry_union,
            'intersection': self.ui.geometry.toolbutton_geometry_intersect,
            'difference':   self.ui.geometry.toolbutton_geometry_difference,
            }
        self.visual_props = copy.deepcopy(DEFAULT_VISUAL_PROPS)

        self.rectilinear_grid = vtk.vtkRectilinearGrid()
        self.grid_viewer_dict = {
            'filters': [],
            'mappers': [],
            'actors':  [],
            }

        self.cell_spacing_widgets = [
            self.parent.ui.mesh.lineedit_mesh_cells_size_x,
            self.parent.ui.mesh.lineedit_mesh_cells_size_y,
            self.parent.ui.mesh.lineedit_mesh_cells_size_z
            ]

        # --- layout ---
        self.grid_layout = QtWidgets.QGridLayout(self)
        self.grid_layout.setContentsMargins(0, 0, 0, 0)

        self.button_bar = QtWidgets.QWidget(self)
        self.button_bar_layout = QtWidgets.QHBoxLayout(self.button_bar)
        self.button_bar_layout.setContentsMargins(0, 0, 0, 0)
        self.button_bar.setLayout(self.button_bar_layout)
        self.button_bar.setGeometry(QtCore.QRect(0, 0, 300, 300))
        self.grid_layout.addWidget(self.button_bar, 0, 0)

        self.vtkWindowWidget = QVTKRenderWindowInteractor(self)
        self.grid_layout.addWidget(self.vtkWindowWidget, 1, 0)

        # enable parameters
        for widget in widget_iter(self.ui.geometry.groupBoxGeometryParameters):
            if isinstance(widget, LineEdit):
                widget.allow_parameters = True

        # --- setup vtk stuff ---
        self.vtkrenderer = vtk.vtkRenderer()
        self.vtkrenderer.GradientBackgroundOn()
        self.vtkrenderer.SetBackground(0.4, 0.4, 0.4)
        self.vtkrenderer.SetBackground2(1.0, 1.0, 1.0)

        self.vtkRenderWindow = self.vtkWindowWidget.GetRenderWindow()
        self.vtkRenderWindow.AddRenderer(self.vtkrenderer)
        self.vtkiren = self.vtkWindowWidget.GetRenderWindow().GetInteractor()

        self.style = CustomInteractorStyle()
        self.style.SetDefaultRenderer(self.vtkrenderer)
        self.vtkiren.SetInteractorStyle(self.style)

        # Orientation Arrows Marker Widget
        self.axes = vtk.vtkAxesActor()
        self.axes.AxisLabelsOn()
        self.axes.SetXAxisLabelText("X")
        self.axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1, 0, 0)
        self.axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
        self.axes.SetYAxisLabelText("Y")
        self.axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0, 1, 0)
        self.axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
        self.axes.SetZAxisLabelText("Z")
        self.axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0, 0, 1)
        self.axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()

        # Orientation Marker Widget
        self.orientation_widget = vtk.vtkOrientationMarkerWidget()
        self.orientation_widget.SetOutlineColor(0.9300, 0.5700, 0.1300)
        self.orientation_widget.SetOrientationMarker(self.axes)
        self.orientation_widget.SetInteractor(self.vtkiren)
        self.orientation_widget.SetViewport(0.0, 0.0, 0.2, 0.2)
        self.orientation_widget.SetEnabled(1)
        self.orientation_widget.InteractiveOff()

        # --- balloon widget ---
        # there seems to be issues with this widget, text doesn't show and the
        # interactor behaves differently...
#        self.balloon_rep = vtk.vtkBalloonRepresentation()
#        self.balloon_rep.SetBalloonLayoutToImageRight()
#        self.balloon_rep.GetTextProperty().SetColor(1,1,1)
#
#        self.balloon_widget = vtk.vtkBalloonWidget()
#        self.balloon_widget.SetInteractor(self.vtkiren)
#        self.balloon_widget.SetRepresentation(self.balloon_rep)
#        self.balloon_widget.EnabledOn()

        # --- setup vtk mappers/actors ---
        self.mesh = None
        self.mesh_actor = None
        self.mesh_mapper = None

        # set types on lineedits
        for widget in widget_iter(self.ui.geometry.stackedWidgetGeometryDetails):
            if isinstance(widget, LineEdit):
                parameter = str(widget.objectName()).lower().split('_')

                if any([s in parameter for s in ['divisions', 'resolution',
                                                 'nhills', 'iterations']]):
                    widget.dtype = int
                else:
                    widget.dtype = float

        self.__connect_events()
        self.__add_tool_buttons()

    # --- setup ---
    def __connect_events(self):

        # --- Geometry Tree ---
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
        self.add_geometry_menu = QtWidgets.QMenu(self)
        self.ui.geometry.toolbutton_add_geometry.setMenu(
            self.add_geometry_menu)

        action = QtWidgets.QAction('STL File',  self.add_geometry_menu)
        action.triggered.connect(self.add_stl)
        self.add_geometry_menu.addAction(action)

        self.add_geometry_menu.addSeparator()

        for geo in PRIMITIVE_DICT.keys():
            action = QtWidgets.QAction(geo, self.add_geometry_menu)
            action.triggered.connect(partial(self.add_primitive, primtype=geo))
            self.add_geometry_menu.addAction(action)

        self.add_geometry_menu.addSeparator()

        for geo in PARAMETRIC_DICT.keys():
            action = QtWidgets.QAction(geo.replace('_', ' '), self.add_geometry_menu)
            action.triggered.connect(partial(self.add_parametric, paramtype=geo))
            self.add_geometry_menu.addAction(action)

        # --- filter button ---
        self.add_filter_menu = QtWidgets.QMenu(self)
        self.ui.geometry.toolbutton_add_filter.setMenu(self.add_filter_menu)

        for geo in FILTER_DICT.keys():
            action = QtWidgets.QAction(geo.replace('_', ' '),
                                       self.add_filter_menu)
            action.triggered.connect(partial(self.add_filter, filtertype=geo))
            self.add_filter_menu.addAction(action)

        # setup signals
        self.ui.geometry.toolbutton_remove_geometry.pressed.connect(
            self.remove_geometry)
        self.ui.geometry.toolbutton_copy_geometry.pressed.connect(
            self.copy_geometry)

        # connect boolean
        for key, btn in self.booleanbtndict.items():
            btn.pressed.connect(partial(self.boolean_operation, booltype=key))

        # connect parameter widgets
        for widget in widget_iter(
                self.ui.geometry.stackedWidgetGeometryDetails):
            if isinstance(widget, QtWidgets.QLineEdit):
                widget.editingFinished.connect(partial(self.parameter_edited, widget))
            elif isinstance(widget, QtWidgets.QCheckBox):
                widget.stateChanged.connect(partial(self.parameter_edited, widget))

        # --- mesh ---
        # connect mesh tab btns
        for i, btn in enumerate([self.ui.mesh.pushbutton_mesh_uniform,
                                 self.ui.mesh.pushbutton_mesh_mesher]):
            btn.pressed.connect(partial(self.change_mesh_tab, i, btn))

        self.ui.geometry.pushbutton_mesh_autosize.pressed.connect(
            self.auto_size_mesh_extents)

        # connect mesher
        self.ui.mesh.combobox_mesher.currentIndexChanged.connect(
                    self.change_mesher_options)

        self.ui.mesh.pushbutton_generate_mesh.pressed.connect(self.mesher)
        self.ui.mesh.pushbutton_remove_mesh.pressed.connect(self.remove_mesh)

    def __add_tool_buttons(self):

        self.toolbutton_reset = QtWidgets.QToolButton()
        self.toolbutton_reset.pressed.connect(self.reset_view)
        self.toolbutton_reset.setIcon(get_icon('overscan.png'))

        self.toolbutton_perspective = QtWidgets.QToolButton()
        self.toolbutton_perspective.pressed.connect(self.perspective)
        self.toolbutton_perspective.setIcon(get_icon('perspective.png'))

        self.toolbutton_view_xy = QtWidgets.QToolButton()
        self.toolbutton_view_xy.pressed.connect(lambda: self.set_view('xy'))
        self.toolbutton_view_xy.setIcon(get_icon('-z.png'))

        self.toolbutton_view_yz = QtWidgets.QToolButton()
        self.toolbutton_view_yz.pressed.connect(lambda: self.set_view('yz'))
        self.toolbutton_view_yz.setIcon(get_icon('-x.png'))

        self.toolbutton_view_xz = QtWidgets.QToolButton()
        self.toolbutton_view_xz.pressed.connect(lambda: self.set_view('xz'))
        self.toolbutton_view_xz.setIcon(get_icon('-y.png'))

        self.toolbutton_screenshot = QtWidgets.QToolButton()
        self.toolbutton_screenshot.pressed.connect(self.screenshot)
        self.toolbutton_screenshot.setIcon(get_icon('camera.png'))

        self.toolbutton_visible = QtWidgets.QToolButton()
        self.toolbutton_visible.setIcon(get_icon('visibility.png'))

        self.visible_menu = QtWidgets.QMenu(self)
        self.toolbutton_visible.setMenu(self.visible_menu)
        self.toolbutton_visible.setPopupMode(
            QtWidgets.QToolButton.InstantPopup)

        # --- visual representation menu ---
        layout = QtWidgets.QGridLayout(self.visible_menu)
        layout.setContentsMargins(5, 5, 5, 5)
        self.visual_btns = {}
        for i, geo in enumerate(['Background Mesh', 'Mesh', 'Geometry',
                                 'Regions', ]):
            geo_name = geo
            geo = geo.lower().replace(' ', '_')
            btns=self.visual_btns[geo]={}
            # tool button
            toolbutton = QtWidgets.QToolButton()
            toolbutton.pressed.connect(partial(self.change_visibility, geo, toolbutton))
            toolbutton.setCheckable(True)
            toolbutton.setChecked(True)
            toolbutton.setAutoRaise(True)
            toolbutton.setIcon(get_icon('visibility.png'))
            layout.addWidget(toolbutton, i, 0)
            btns['visible']=toolbutton

            # style
            combobox = QtWidgets.QComboBox()
            combobox.addItems(['wire', 'solid', 'edges', 'points'])
            combobox.activated.connect(partial(self.change_representation, geo, combobox))
            layout.addWidget(combobox, i, 1)
            btns['rep'] = combobox

            # color
            if not geo=='regions':
                toolbutton = QtWidgets.QToolButton()
                toolbutton.pressed.connect(partial(self.change_color, geo, toolbutton))
                toolbutton.setAutoRaise(True)
                layout.addWidget(toolbutton, i, 2)
                btns['color']=toolbutton

            # opacity
            slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
            slider.setRange(0, 100)
            slider.setFixedWidth(40)
            slider.sliderReleased.connect(partial(self.change_opacity, geo, slider))
            layout.addWidget(slider, i, 3)
            btns['opacity']=slider

            # label
            label = QtWidgets.QLabel(geo_name)
            layout.addWidget(label, i, 4)
        self.set_visual_btn_values()

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

    def set_visual_btn_values(self):
        for geo, info in self.visual_props.items():
            for key, value in info.items():
                if geo in self.visual_btns and key in self.visual_btns[geo]:
                    wid = self.visual_btns[geo][key]
                    if key=='rep':
                        wid.setCurrentIndex(wid.findText(value))
                    elif key=='visible':
                        wid.setChecked(value)
                        self.set_visible_btn_image(wid, value)
                    elif key=='color':
                        wid.setStyleSheet("QToolButton{{ background: {};}}".format(value.name()))
                    elif key=='opacity':
                        wid.setValue(value*100)

    def emitUpdatedValue(self, key, value, args=None):
        self.value_updated.emit(self, {key: value}, args)

    def updateValue(self, key, newValue, args=None):

        if key in ['xmin', 'xlength', 'ymin', 'ylength', 'zmin', 'zlength',
                   'imax', 'jmax', 'kmax']:
            self.update_mesh()
        elif key == 'no_k':
            self.ui.geometry.lineedit_keyword_zlength.setEnabled(not newValue)
            self.ui.mesh.lineedit_keyword_kmax.setEnabled(not newValue)

    def objectName(self):
        return 'VTK Widget'

    def default(self):
        """reset to defaults"""
        self.ui.geometry.lineedit_keyword_zlength.setEnabled(True)
        self.ui.mesh.lineedit_keyword_kmax.setEnabled(True)
        self.vtkrenderer.RemoveAllViewProps()
        self.clear_all_geometry()
        self.render()

    # --- save/load ---
    def geometry_to_str(self):
        """convert geometry to string"""

        # save tree
        tree = {}
        itr = QtWidgets.QTreeWidgetItemIterator(self.geometrytree)
        while itr.value():
            item = itr.value()
            text = item.text(0)
            if text not in tree:
                tree[text] = []
            for i in range(item.childCount()):
                tree[text].append(item.child(i).text(0))
            itr += 1

        data = {'geometry_dict': clean_geo_dict(self.geometrydict),
                'tree': tree}

        return ExtendedJSON.dumps(data)

    def geometry_from_str(self, string):
        """convert string to geometry"""
        try:
            data = ExtendedJSON.loads(string)
            tree = data['tree']
            geo_dict = data['geometry_dict']
        except Exception as e:
            self.parent.message(text='Error loading geometry: %s' % e)
            return

        if not tree or not data:
            return

        # convert lists to sets
        for key, value in tree.items():
            tree[key] = set(tree[key])

        # build geometry
        for nodes in topological_sort(tree):
            for node in nodes:
                if node in geo_dict and 'geo_type' in geo_dict[node]:
                    geo_data = copy.deepcopy(DEFAULT_PARAMS[geo_dict[node]['geo_type']])
                    geo_data.update(geo_dict[node])
                    if geo_dict[node]['geo_type'] == 'primitive':
                        self.add_primitive(name=node, data=geo_data)
                    elif geo_dict[node]['geo_type'] == 'parametric':
                        self.add_parametric(name=node, data=geo_data)
                    elif geo_dict[node]['geo_type'] == 'filter':
                        self.add_filter(name=node, data=geo_data,
                                        child=tree[node].pop())
                    elif geo_dict[node]['geo_type'] =='boolean':
                        self.boolean_operation(boolname=node, data=geo_data,
                                               children=tree[node])
                    elif geo_dict[node]['geo_type'] =='stl':
                        self.add_stl(None, filename=geo_dict[node]['filename'],
                                     name=node, data=geo_data)
                else:
                    self.parent.message(text='Error loading geometry: Geometry does not have parameters.')
                    return

    def visual_props_to_str(self):
        """convert visual props to str"""
        return ExtendedJSON.dumps(clean_visual_dict(self.visual_props))

    def visual_props_from_str(self, string):
        """convert string to visual props"""
        data = ExtendedJSON.loads(string)

        for geo, geo_dict in data.items():
            for key, value in geo_dict.items():
                if key in ['color', 'edge']:
                    data[geo][key] = QtGui.QColor(*value)

        self.visual_props = copy.deepcopy(DEFAULT_VISUAL_PROPS)
        self.visual_props.update(data)
        if not self.mesh_actor is None:
            self.set_mesh_actor_props()

        self.set_visual_btn_values()
        for actor in self.grid_viewer_dict['actors']:
            self.set_background_mesh_actor_props(actor)

    # --- render ---
    def render(self, force_render=False, defer_render=None):
        if defer_render is not None:
            self.defer_render = defer_render
        if not self.defer_render or force_render:
            self.vtkRenderWindow.Render()

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
                self.ui.geometry.treeWidgetGeometry.indexOfTopLevelItem(
                current_selection[0]) > -1:

            self.ui.geometry.toolbutton_remove_geometry.setEnabled(True)
            self.ui.geometry.toolbutton_add_filter.setEnabled(True)
            self.ui.geometry.toolbutton_copy_geometry.setEnabled(True)
        else:
            self.ui.geometry.toolbutton_remove_geometry.setEnabled(False)
            self.ui.geometry.toolbutton_add_filter.setEnabled(False)
            self.ui.geometry.toolbutton_copy_geometry.setEnabled(False)

        if current_selection:
            text = str(current_selection[-1].text(0)).lower()

            current_index = 0
            for i in range(
                    self.ui.geometry.stackedWidgetGeometryDetails.count()):
                widget = self.ui.geometry.stackedWidgetGeometryDetails.widget(
                                                                             i)
                if str(widget.objectName()) == self.geometrydict[text]['type']:
                    current_index = i
                    break

            # set the widget parameters
            for child in widget_iter(widget):
                name = str(child.objectName()).lower().replace('_', '')
                for key, value in self.geometrydict[text].items():
                    if key in name:
                        break

                if isinstance(child, LineEdit):
                    child.updateValue(None, value)
                elif isinstance(child, QtWidgets.QCheckBox):
                    child.setChecked(value)

            self.ui.geometry.groupBoxGeometryParameters.setTitle(text)

        else:
            current_index = 0

            self.ui.geometry.groupBoxGeometryParameters.setTitle('Parameters')
            self.ui.geometry.toolbutton_remove_geometry.setEnabled(False)

        self.parent.animate_stacked_widget(
            self.ui.geometry.stackedWidgetGeometryDetails,
            self.ui.geometry.stackedWidgetGeometryDetails.currentIndex(),
            current_index,
            'horizontal',
            )

    def get_tree_item(self, name):
        """return the tree item with name"""
        itr = QtWidgets.QTreeWidgetItemIterator(self.geometrytree)
        while itr.value():
            item = itr.value()
            if name == item.text(0):
                return item
            itr += 1

    def geometry_clicked(self, item):
        """
        Hide/Show the clicked geometry
        """

        name = str(item.text(0)).lower()
        if item.checkState(0) == QtCore.Qt.Checked:
            if self.visual_props['geometry']['visible']:
                self.geometrydict[name]['actor'].VisibilityOn()
            self.geometrydict[name]['visible'] = True
        else:
            self.geometrydict[name]['actor'].VisibilityOff()
            self.geometrydict[name]['visible'] = False

        self.render()

    def add_stl(self, widget, filename=None, name=None, data=None):
        """
        Open browse dialog and load selected stl file
        """

        if filename is None:
            filename = QtWidgets.QFileDialog.getOpenFileName(
                self, 'Select an STL File',
                self.parent.get_project_dir(),
                'STL File (*.stl)',)

            if isinstance(filename, (tuple, list)):
                filename = filename[0]

            filename = str(filename)

        if filename:
            if name is None:
                name = os.path.basename(filename).lower()
                name = get_unique_string(name, list(self.geometrydict.keys()))

            if data is not None:
                self.geometrydict[name] = data
            else:
                self.geometrydict[name] = copy.deepcopy(DEFAULT_STL_PARAMS)
            self.geometrydict[name]['filename'] = filename

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

            self.set_geometry_actor_props(actor, name)

            self.vtkrenderer.AddActor(actor)

            self.render()

            # find center of mass
            center_filter = vtk.vtkCenterOfMass()
            center_filter.SetInputConnection(transform_filter.GetOutputPort())
            center_filter.SetUseScalarsAsWeights(False)
            center_filter.Update()
            center = center_filter.GetCenter()

            # Add to dict
            self.geometrydict[name]['reader'] = reader
            self.geometrydict[name]['transform'] = transform
            self.geometrydict[name]['transformfilter'] = transform_filter
            self.geometrydict[name]['center_filter'] = center_filter
            self.geometrydict[name]['mapper'] = mapper
            self.geometrydict[name]['actor'] = actor
            self.geometrydict[name]['centerx'] = center[0]
            self.geometrydict[name]['centery'] = center[1]
            self.geometrydict[name]['centerz'] = center[2]

            # Add to tree
            item = QtWidgets.QTreeWidgetItem([name])
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(0, QtCore.Qt.Checked)
            self.geometrytree.addTopLevelItem(item)
            self.geometrytree.setCurrentItem(item)

            self.parent.set_unsaved_flag()

    def parameter_edited(self, widget, name=None, value=None, key=None):
        """
        Update the value of edited parameter in the geometrydict
        """
        if name is None:
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
            if isinstance(widget, LineEdit):
                value = widget.value
            elif isinstance(widget, QtWidgets.QCheckBox):
                value = widget.isChecked()

        if value is not None:
            self.update_parameter_map(value, name, key)
            self.geometrydict[name][key] = value

            if self.geometrydict[name]['type'] in \
                    list(PRIMITIVE_DICT.keys()):
                self.update_primitive(name)
            elif self.geometrydict[name]['type'] in \
                    list(FILTER_DICT.keys()):
                self.update_filter(name)
            elif self.geometrydict[name]['type'] in \
                    list(PARAMETRIC_DICT.keys()):
                self.update_parametric(name)

            if 'transform' in self.geometrydict[name]:
                self.update_transform(name)
            self.render()

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
            source.SetRadius(safe_float(self.geometrydict[name]['radius']))
            source.SetThetaResolution(int(
                self.geometrydict[name]['thetaresolution']))
            source.SetPhiResolution(int(
                self.geometrydict[name]['phiresolution']))

        elif primtype == 'box':
            source.SetXLength(safe_float(self.geometrydict[name]['lengthx']))
            source.SetYLength(safe_float(self.geometrydict[name]['lengthy']))
            source.SetZLength(safe_float(self.geometrydict[name]['lengthz']))

        elif primtype == 'cone':
            source.SetRadius(safe_float(self.geometrydict[name]['radius']))
            source.SetHeight(safe_float(self.geometrydict[name]['height']))
            source.SetDirection(safe_float(self.geometrydict[name]['directionx']),
                                safe_float(self.geometrydict[name]['directiony']),
                                safe_float(self.geometrydict[name]['directionz']))
            source.SetResolution(int(self.geometrydict[name]['resolution']))
            source.CappingOn()

        elif primtype == 'cylinder':
            source.SetRadius(safe_float(self.geometrydict[name]['radius']))
            source.SetHeight(safe_float(self.geometrydict[name]['height']))
            source.SetResolution(int(self.geometrydict[name]['resolution']))

        elif primtype == 'stl':
            pass

        else:
            return

        # common props
        if source is not None:
            source.SetCenter(safe_float(self.geometrydict[name]['centerx']),
                             safe_float(self.geometrydict[name]['centery']),
                             safe_float(self.geometrydict[name]['centerz']))

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
        transform.Translate(-safe_float(self.geometrydict[name]['centerx']),
                            -safe_float(self.geometrydict[name]['centery']),
                            -safe_float(self.geometrydict[name]['centerz']))

        # rotation
        transform.RotateWXYZ(safe_float(self.geometrydict[name]['rotationx']), 1, 0, 0)
        transform.RotateWXYZ(safe_float(self.geometrydict[name]['rotationy']), 0, 1, 0)
        transform.RotateWXYZ(safe_float(self.geometrydict[name]['rotationz']), 0, 0, 1)

        # back to position
        transform.Translate(safe_float(self.geometrydict[name]['centerx']),
                            safe_float(self.geometrydict[name]['centery']),
                            safe_float(self.geometrydict[name]['centerz']))

        # translate stl files
        if self.geometrydict[name]['type'] in ['stl'] + \
                list(PARAMETRIC_DICT.keys()):
            transform.Translate(
                safe_float(self.geometrydict[name]['translationx']),
                safe_float(self.geometrydict[name]['translationy']),
                safe_float(self.geometrydict[name]['translationz']),
                )

        # update
        transform_filter.Update()

        return transform_filter

    def add_primitive(self, primtype='sphere', name=None, data=None):
        """
        Add the specified primitive
        """
        if name is None:
            name = get_unique_string(primtype, list(self.geometrydict.keys()))

        # create primitive
        if data is not None:
            primtype = data['type']
        if primtype in PRIMITIVE_DICT:
            source = PRIMITIVE_DICT[primtype]()
        else:
            return

        if data is None:
            self.geometrydict[name] = copy.deepcopy(DEFAULT_PRIMITIVE_PARAMS)
            self.geometrydict[name]['type'] = primtype
        else:
            self.geometrydict[name] = data
        self.geometrydict[name]['source'] = source

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

        self.set_geometry_actor_props(actor, name)

        self.vtkrenderer.AddActor(actor)

        self.render()

        # add to dict
        self.geometrydict[name]['mapper'] = mapper
        self.geometrydict[name]['actor'] = actor
        self.geometrydict[name]['trianglefilter'] = trianglefilter
        self.geometrydict[name]['source'] = source

        # Add to tree
        item = QtWidgets.QTreeWidgetItem([name])
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(0, QtCore.Qt.Checked)
        self.geometrytree.addTopLevelItem(item)
        self.geometrytree.setCurrentItem(item)

        self.parent.set_unsaved_flag()

    def update_parametric(self, name):
        """
        Update the specified parameteric object.
        """
        paratype = self.geometrydict[name]['type']
        para_object = self.geometrydict[name]['parametric_object']
        source = self.geometrydict[name]['source']

        if paratype == 'torus':
            para_object.SetRingRadius(safe_float(
                self.geometrydict[name]['ringradius']))
            para_object.SetCrossSectionRadius(safe_float(
                self.geometrydict[name]['crosssectionradius']))
        elif paratype == 'boy':
            para_object.SetZScale(safe_float(self.geometrydict[name]['zscale']))
        elif paratype == 'conic_spiral':
            para_object.SetA(safe_float(self.geometrydict[name]['ascale']))
            para_object.SetB(safe_float(self.geometrydict[name]['bfunc']))
            para_object.SetC(safe_float(self.geometrydict[name]['cfunc']))
            para_object.SetN(safe_float(self.geometrydict[name]['nfunc']))
        elif paratype == 'dini':
            para_object.SetA(safe_float(self.geometrydict[name]['ascale']))
            para_object.SetB(safe_float(self.geometrydict[name]['bscale']))
        elif paratype == 'ellipsoid':
            para_object.SetXRadius(safe_float(self.geometrydict[name]['radiusx']))
            para_object.SetYRadius(safe_float(self.geometrydict[name]['radiusy']))
            para_object.SetZRadius(safe_float(self.geometrydict[name]['radiusz']))
        elif paratype == 'figure_8_klein':
            para_object.SetRadius(safe_float(self.geometrydict[name]['radius']))
        elif paratype == 'mobius':
            para_object.SetRadius(safe_float(self.geometrydict[name]['radius']))
        elif paratype == 'random_hills':
            para_object.SetHillXVariance(safe_float(
                self.geometrydict[name]['variancex']))
            para_object.SetXVarianceScaleFactor(safe_float(
                self.geometrydict[name]['scalex']))
            para_object.SetHillYVariance(safe_float(
                self.geometrydict[name]['variancey']))
            para_object.SetYVarianceScaleFactor(safe_float(
                self.geometrydict[name]['scaley']))
            para_object.SetHillAmplitude(safe_float(
                self.geometrydict[name]['amplitude']))
            para_object.SetAmplitudeScaleFactor(safe_float(
                self.geometrydict[name]['scaleamplitude']))
            para_object.SetNumberOfHills(int(
                self.geometrydict[name]['nhills']))
            if self.geometrydict[name]['allowrandom']:
                para_object.AllowRandomGenerationOn()
            else:
                para_object.AllowRandomGenerationOff()
        elif paratype == 'roman':
            para_object.SetRadius(safe_float(self.geometrydict[name]['radius']))
        elif paratype == 'super_ellipsoid':
            para_object.SetXRadius(safe_float(self.geometrydict[name]['radiusx']))
            para_object.SetYRadius(safe_float(self.geometrydict[name]['radiusy']))
            para_object.SetZRadius(safe_float(self.geometrydict[name]['radiusz']))
            para_object.SetN1(safe_float(self.geometrydict[name]['n1']))
            para_object.SetN2(safe_float(self.geometrydict[name]['n2']))
        elif paratype == 'super_toroid':
            para_object.SetXRadius(safe_float(self.geometrydict[name]['radiusx']))
            para_object.SetYRadius(safe_float(self.geometrydict[name]['radiusy']))
            para_object.SetZRadius(safe_float(self.geometrydict[name]['radiusz']))
            para_object.SetRingRadius(safe_float(
                self.geometrydict[name]['ringradius']))
            para_object.SetCrossSectionRadius(safe_float(
                self.geometrydict[name]['crosssectionradius']))
            para_object.SetN1(safe_float(self.geometrydict[name]['n1']))
            para_object.SetN2(safe_float(self.geometrydict[name]['n2']))

        source.Update()

        self.parent.set_unsaved_flag()

        return source

    def add_parametric(self, paramtype=None, name=None, data=None):
        """
        Add the specified parametric object
        """

        if paramtype is None:
            paramtype = data['type']

        if name is None:
            name = get_unique_string(paramtype, list(self.geometrydict.keys()))

        parametric_object = PARAMETRIC_DICT[paramtype]()
        source = vtk.vtkParametricFunctionSource()
        source.SetParametricFunction(parametric_object)

        if data is None:
            self.geometrydict[name] = copy.deepcopy(DEFAULT_PARAMETRIC_PARAMS)
            self.geometrydict[name]['type'] = paramtype
        else:
            self.geometrydict[name] = data
        self.geometrydict[name]['parametric_object'] = parametric_object
        self.geometrydict[name]['source'] = source

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

        self.set_geometry_actor_props(actor, name)

        self.vtkrenderer.AddActor(actor)

        self.render()

        # add to dict
        self.geometrydict[name]['mapper'] = mapper
        self.geometrydict[name]['actor'] = actor
        self.geometrydict[name]['trianglefilter'] = trianglefilter
        self.geometrydict[name]['source'] = source

        # Add to tree
        item = QtWidgets.QTreeWidgetItem([name])
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(0, QtCore.Qt.Checked)
        self.geometrytree.addTopLevelItem(item)
        self.geometrytree.setCurrentItem(item)

    def boolean_operation(self, booltype=None, boolname=None, data=None, children=None):
        """
        Apply a boolean operation with the currently selected toplevel items.
        """

        if children is not None:
            current_selection = []
            for child in children:
                current_selection.append(self.get_tree_item(child))
        else:
            current_selection = self.geometrytree.selectedItems()
            if len(current_selection) != 2:
                return

        if boolname is None:
            boolname = get_unique_string(booltype, list(self.geometrydict.keys()))

        # Save references
        if data is not None:
            self.geometrydict[boolname] = data
            booltype = data['type']
        else:
            self.geometrydict[boolname] = copy.deepcopy(DEFAULT_BOOLEAN_PARAMS)
            self.geometrydict[boolname]['type'] = booltype

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

            geometry = self.get_input_data(name)
            boolean_operation.SetInputConnection(
                i, geometry.GetOutputPort())

            # hide the sources
            self.geometrydict[name]['actor'].VisibilityOff()
            self.geometrydict[name]['visible'] = False

        boolean_operation.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(boolean_operation.GetOutputPort())
        mapper.ScalarVisibilityOff()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        self.set_geometry_actor_props(actor, boolname)

        self.vtkrenderer.AddActor(actor)

        self.render()

        # save references
        self.geometrydict[boolname]['booleanoperation'] = boolean_operation
        self.geometrydict[boolname]['mapper'] = mapper
        self.geometrydict[boolname]['actor'] = actor

        # Add to tree
        toplevel = QtWidgets.QTreeWidgetItem([boolname])
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

        self.parent.set_unsaved_flag()

    def clear_all_geometry(self):
        """ remove all geometry """
        self.geometrytree.clear()
        self.geometrydict = {}

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
                if self.visual_props['geometry']['visible']:
                    name = str(child.text(0)).lower()
                    self.geometrydict[name]['actor'].VisibilityOn()
                    self.geometrydict[name]['visible'] = True
                child.setCheckState(0, QtCore.Qt.Checked)

            # remove graphics
            geo = self.geometrydict.pop(text)
            self.vtkrenderer.RemoveActor(geo['actor'])

            self.render()

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
        elif filtertype == 'smooth':
            vtkfilter.SetEdgeAngle(self.geometrydict[name]['edgeangle'])
            vtkfilter.SetFeatureAngle(self.geometrydict[name]['featureangle'])
            vtkfilter.SetNumberOfIterations(
                self.geometrydict[name]['iterations'])
            vtkfilter.SetRelaxationFactor(
                self.geometrydict[name]['relaxation'])

            if not self.geometrydict[name]['featureedgesmoothing']:
                vtkfilter.FeatureEdgeSmoothingOn()
            else:
                vtkfilter.FeatureEdgeSmoothingOff()

            if self.geometrydict[name]['boundarysmoothing']:
                vtkfilter.BoundarySmoothingOn()
            else:
                vtkfilter.BoundarySmoothingOff()

        elif filtertype == 'windowed_sinc':
            vtkfilter.SetEdgeAngle(self.geometrydict[name]['edgeangle'])
            vtkfilter.SetFeatureAngle(self.geometrydict[name]['featureangle'])
            vtkfilter.SetNumberOfIterations(
                self.geometrydict[name]['iterations'])
            vtkfilter.SetPassBand(self.geometrydict[name]['passband'])

            if not self.geometrydict[name]['featureedgesmoothing']:
                vtkfilter.FeatureEdgeSmoothingOn()
            else:
                vtkfilter.FeatureEdgeSmoothingOff()

            if self.geometrydict[name]['boundarysmoothing']:
                vtkfilter.BoundarySmoothingOn()
            else:
                vtkfilter.BoundarySmoothingOff()

            if self.geometrydict[name]['manifoldsmoothing']:
                vtkfilter.NonManifoldSmoothingOn()
            else:
                vtkfilter.NonManifoldSmoothingOff()

            if self.geometrydict[name]['normalize']:
                vtkfilter.NormalizeCoordinatesOn()
            else:
                vtkfilter.NormalizeCoordinatesOff()

        vtkfilter.Update()

    def add_filter(self, filtertype=None, name=None, data=None, child=None):
        """
        add the selected filter with the input being the currently selected
        toplevel item
        """

        if child is not None:
            selection_text = child
            current_selection = [self.get_tree_item(child)]
        else:
            current_selection = self.geometrytree.selectedItems()
            if current_selection:
                selection_text = str(current_selection[-1].text(0)).lower()

                name = get_unique_string(
                    filtertype, list(self.geometrydict.keys()))
            else:
                return

        if data is None:
            self.geometrydict[name] = copy.deepcopy(DEFAULT_FILTER_PARAMS)
            self.geometrydict[name]['type'] = filtertype
        else:
            self.geometrydict[name] = data

        # init filter
        if data is not None:
            filtertype = data['type']
        vtkfilter = FILTER_DICT[filtertype]()
        self.geometrydict[name]['filter'] = vtkfilter

        # set input data
        inputdata = self.get_input_data(selection_text)
        vtkfilter.SetInputConnection(inputdata.GetOutputPort())

        # update filter
        self.update_filter(name)

        # hide the source
        self.geometrydict[selection_text]['actor'].VisibilityOff()
        self.geometrydict[selection_text]['visible'] = False

        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(vtkfilter.GetOutputPort())
        mapper.ScalarVisibilityOff()

        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        self.set_geometry_actor_props(actor, name)

        # add actor to render
        self.vtkrenderer.AddActor(actor)

        # update
        self.render()

        # save references
        self.geometrydict[name]['actor'] = actor
        self.geometrydict[name]['mapper'] = mapper

        # Add to tree
        toplevel = QtWidgets.QTreeWidgetItem([name])
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

        self.parent.set_unsaved_flag()

    def get_input_data(self, name):
        """ based on the type of geometry, return the data """

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

        geo = None
        if append_filter.GetTotalNumberOfInputConnections() > 0:
            append_filter.Update()

            # clean
            clean_filter = vtk.vtkCleanPolyData()
            clean_filter.SetInputConnection(append_filter.GetOutputPort())
            clean_filter.Update()

            geo = clean_filter

        return geo

    def get_geometry_extents(self):
        """ determine the extents of the visible geometry """

        geometry = self.collect_toplevel_geometry()

        bounds = None
        if geometry:
            bounds = geometry.GetOutput().GetBounds()

        return bounds

    def set_geometry_actor_props(self, actor, name):
        """ set the geometry proprerties to the others in the scene """

        props = self.visual_props['geometry']

        self.set_representation(actor, props['rep'])
        actor.GetProperty().SetColor(props['color'].getRgbF()[:3])
        actor.GetProperty().SetEdgeColor(props['edge'].getRgbF()[:3])
        actor.GetProperty().SetOpacity(props['opacity'])

        # check visibility
        if not props['visible']:
            actor.VisibilityOff()

    # --- parameters ---
    def update_parameter_map(self, new_value, name, key):
        """update the mapping of parameters and keywords"""

        data = self.geometrydict
        name_key = ','.join([name, key])

        # new params
        new_params = []
        if isinstance(new_value, Equation):
            new_params = new_value.get_used_parameters()

        # old params
        old_value = data[name][key]

        old_params = []
        if isinstance(old_value, Equation):
            old_params = old_value.get_used_parameters()

        add = set(new_params)-set(old_params)
        for param in add:
            if param not in self.parameter_key_map:
                self.parameter_key_map[param] = set()
            self.parameter_key_map[param].add(name_key)

        remove = set(old_params)-set(new_params)
        for param in remove:
            self.parameter_key_map[param].remove(name_key)
            if len(self.parameter_key_map[param]) == 0:
                self.parameter_key_map.pop(param)

    def update_parameters(self, params):
        """parameters have changed, update regions"""
        data = self.geometrydict
        self.defer_render = True
        for param in params:
            if param in self.parameter_key_map:
                for var in self.parameter_key_map[param]:
                    name, key = var.split(',')
                    value = data[name][key]
                    self.parameter_edited(None, name, value, key)
        self.render(defer_render=False)

    # --- regions ---
    def update_region_source(self, name):
        """
        Update the specified primitive
        """
        primtype = self.region_dict[name]['type']
        props = self.region_dict[name]

        if 'source' in self.region_dict[name]:
            source = self.region_dict[name]['source']
        else:
            source = None


        lengths = [abs(safe_float(to) - safe_float(from_)) for
                   from_, to in zip(props['from'], props['to'])]
        center = [min(map(safe_float, ft)) + l / 2.0 for ft, l in
                  zip(zip(props['from'], props['to']),
                      lengths)]


        # update source
        if primtype == 'sphere':
            source.SetRadius(min(lengths)/2.0)
        elif primtype == 'point':
            source.SetRadius(.01)
            center = props['from']
        elif primtype == 'box':
            source.SetXLength(lengths[0])
            source.SetYLength(lengths[1])
            source.SetZLength(lengths[2])
        elif primtype == 'XY-plane':
            source.SetXLength(lengths[0])
            source.SetYLength(lengths[1])
            source.SetZLength(0)
            center[2] = props['from'][2]
        elif primtype == 'XZ-plane':
            source.SetXLength(lengths[0])
            source.SetYLength(0)
            source.SetZLength(lengths[2])
            center[1] = props['from'][1]
        elif primtype == 'YZ-plane':
            source.SetXLength(0)
            source.SetYLength(lengths[1])
            source.SetZLength(lengths[2])
            center[0] = props['from'][0]
        elif primtype == 'STL':
            source.SetXLength(lengths[0])
            source.SetYLength(lengths[1])
            source.SetZLength(lengths[2])
        else:
            return

        # common props
        if source is not None:
            source.SetCenter(*center)
            source.Update()

        return source

    def new_region(self, name, region):
        self.region_dict[name] = copy.deepcopy(region)

        if region['type'] == 'point':
            shape = 'sphere'
        else:
            shape = 'box'

        source = PRIMITIVE_DICT[shape]()
        self.region_dict[name]['source'] = source
        self.update_region_source(name)

        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())

        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        self.set_region_actor_props(actor, name, region['color'].color_float)

        self.vtkrenderer.AddActor(actor)
#        self.balloon_widget.AddBalloon(actor, name, None)

        self.region_dict[name]['actor'] = actor
        self.region_dict[name]['mapper'] = mapper
        self.select_facets(name)

        self.change_region_visibility(
            name, self.region_dict[name]['visibility'],)

    def delete_region(self, name):
        region = self.region_dict.pop(name)
        self.vtkrenderer.RemoveActor(region['actor'])
        if 'clip_actor' in region:
            self.vtkrenderer.RemoveActor(region['clip_actor'])

    def update_region(self, name, region):
        self.region_dict[name].update(copy.deepcopy(region))
        self.update_region_source(name)
        self.select_facets(name)
        self.render()

    def change_region_color(self, name, color):
        """ change the color of a region """
        self.region_dict[name]['color'] = copy.deepcopy(color)

        actor = self.region_dict[name]['actor']
        actor.GetProperty().SetColor(
            *self.region_dict[name]['color'].color_float)

        if 'clip_actor' in self.region_dict[name]:
            actor = self.region_dict[name]['clip_actor']
            actor.GetProperty().SetColor(
                *self.region_dict[name]['color'].color_float)

        self.render()

    def change_region_type(self, name, region):
        """ change the type of a region """

        self.region_dict[name].update(copy.deepcopy(region))

        if region['type'] == 'point':
            shape = 'sphere'
        else:
            shape = 'box'

        source = PRIMITIVE_DICT[shape]()
        self.region_dict[name]['source'] = source
        self.update_region_source(name)
        self.select_facets(name)

        self.region_dict[name]['mapper'].SetInputConnection(
            source.GetOutputPort())

        self.render()

    def change_region_name(self, old_name, new_name):
        """ change the name of a region """

        region = self.region_dict.pop(old_name)
        self.region_dict[new_name] = region

    def change_region_visibility(self, name, visible):
        """ change the visibility of a region """

        if visible and self.visual_props['regions']['visible']:
            self.region_dict[name]['actor'].VisibilityOn()
            if 'clip_actor' in self.region_dict[name]:
                self.region_dict[name]['clip_actor'].VisibilityOn()
        else:
            self.region_dict[name]['actor'].VisibilityOff()
            if 'clip_actor' in self.region_dict[name]:
                self.region_dict[name]['clip_actor'].VisibilityOff()
        self.region_dict[name]['visible'] = visible

        self.render()

    def set_region_actor_props(self, actor, name, color=None):
        """ set the geometry properties to the others in the scene """

        props = self.visual_props['regions']
        self.set_representation(actor, props['rep'])
        actor.GetProperty().SetOpacity(props['opacity'])

        if color:
            actor.GetProperty().SetColor(*color)
        else:
            actor.GetProperty().SetColor(props['color'].getRgbF()[:3])

        # check visibility
        if not props['visible']:
            actor.VisibilityOff()

    def select_facets(self, name):
        """ select facets with an implicit function """

        region = self.region_dict[name]
        # remove old objects
        if 'clip_actor' in region:
            self.vtkrenderer.RemoveActor(region['clip_actor'])
            for key in ['clip_actor',  'clip_mapper', 'clipper', 'implicit']:
                region.pop(key)

        if region['type'].lower() != 'stl':
            return

        # bounds
        bounds = [0.0]*6
        bounds[::2] = region['from']
        bounds[1::2] = region['to']
        lengths = [abs(safe_float(to) - safe_float(f)) for
                   f, to in zip(region['from'], region['to'])]
        center = [min(f) + l / 2.0 for f, l in
                  zip(zip(region['from'], region['to']),
                      lengths)]

        implicit = None
        if region['stl_shape'] == 'box':
            implicit = vtk.vtkBox()
            implicit.SetBounds(bounds)
        elif region['stl_shape'] == 'ellipsoid':
            trans = vtk.vtkTransform()
            # for some reason, everything is inverted?
            trans.Scale([1/l if l > 0 else 1 for l in lengths])
            trans.Translate([-c for c in center])
            implicit = vtk.vtkSphere()
            implicit.SetTransform(trans)

        region['implicit'] = implicit

        geo = self.collect_toplevel_geometry()
        if implicit is not None and geo is not None:
            if region['slice']:
                clipper = vtk.vtkClipPolyData()
                clipper.SetClipFunction(implicit)
                clipper.SetInputConnection(geo.GetOutputPort())
                clipper.GenerateClippedOutputOn()
                clipper_output = clipper.GetClippedOutputPort()
            else:
                clipper = vtk.vtkExtractPolyDataGeometry()
                clipper.SetImplicitFunction(implicit)
                clipper.SetInputConnection(geo.GetOutputPort())
                clipper.ExtractInsideOn()
                clipper_output = clipper.GetOutputPort()

            clip_mapper = vtk.vtkPolyDataMapper()
            clip_mapper.SetInputConnection(clipper_output)
            clip_mapper.ScalarVisibilityOff()
            clip_actor = vtk.vtkActor()
            clip_actor.SetMapper(clip_mapper)

            region['clipper'] = clipper
            region['clip_mapper'] = clip_mapper
            region['clip_actor'] = clip_actor
            region['actor'].GetProperty().SetRepresentationToWireframe()

            self.set_region_actor_props(clip_actor, name,
                                        color=region['color'].color_float)

            self.vtkrenderer.AddActor(clip_actor)

    # --- output files ---
    def export_stl(self, file_name):
        """ export visible toplevel geometry """

        geometry = self.collect_toplevel_geometry()

        if geometry:
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
            fname = str(QtWidgets.QFileDialog.getSaveFileName(
                self.parent,
                "Save screenshot",
                self.parent.settings.value('project_dir'),
                ';;'.join(["PNG (*.png)",
                           "JPEG (*.jpg)",
                           "PostScript (*.ps)",
                           "All Files (*.*)",
                           ]),
                ))

        if not fname:
            return

        # screenshot code:
        window_image = vtk.vtkWindowToImageFilter()
        window_image.SetInput(self.vtkRenderWindow)
        window_image.SetMagnification(3)
        window_image.SetInputBufferTypeToRGBA()
#        window_image.ReadFrontBufferOff()
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
            self.ui.mesh.stackedwidget_mesh,
            self.ui.mesh.stackedwidget_mesh.currentIndex(),
            tabnum,
            direction='horizontal',
            line=self.ui.mesh.line_mesh,
            to_btn=btn,
            btn_layout=self.ui.mesh.gridlayout_mesh_tab_btns,
            )

    def change_mesher_options(self):
        """ switch the mesh options stacked widget """

        mesher = str(self.ui.mesh.combobox_mesher.currentText()).lower()

        current_index = 0
        for i in range(self.ui.mesh.stackedwidget_mesher_options.count()):
            widget = self.ui.mesh.stackedwidget_mesher_options.widget(i)
            if mesher == str(widget.objectName()).lower():
                current_index = i
                break

        self.parent.animate_stacked_widget(
            self.ui.mesh.stackedwidget_mesher_options,
            self.ui.mesh.stackedwidget_mesher_options.currentIndex(),
            current_index,
            direction='horizontal',
            )

    def auto_size_mesh_extents(self):
        """ collect and set the extents of the visible geometry """
        extents = self.get_geometry_extents()

        if extents:
            for key, extent in zip(['xmin', 'xlength', 'ymin', 'ylength',
                                    'zmin', 'zlength'],
                                   extents):
                if 'min' not in key:  # mfix doesn't support mins yet
                    self.emitUpdatedValue(key, extent)

            self.update_mesh()

    def update_mesh(self):

        extents = []
        for key in ['xmin', 'xlength', 'ymin', 'ylength', 'zmin', 'zlength']:
            try:
                extents.append(safe_float(self.project[key]))
            except TypeError:
                extents.append(0.0)
            except KeyError:
                extents.append(0.0)

        cells = []
        for key in ['imax', 'jmax', 'kmax']:
            if key in self.project:
                cells.append(int(self.project[key])+1)
            else:
                cells.append(1)

        # average cell width
        for (f, t), c, wid in zip(zip(extents[::2], extents[1::2]), cells, self.cell_spacing_widgets):
            wid.setText('{0:.2e}'.format((t-f)/c))
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
            self.set_background_mesh_actor_props(actor)
            self.grid_viewer_dict['actors'].append(actor)

            self.vtkrenderer.AddActor(actor)

        self.render()

    def set_background_mesh_actor_props(self, actor):
        props = self.visual_props['background_mesh']
        self.set_representation(actor, props['rep'])
        actor.GetProperty().SetColor(props['color'].getRgbF()[:3])
        actor.GetProperty().SetEdgeColor(props['color'].getRgbF()[:3])
        actor.GetProperty().SetOpacity(props['opacity'])
        actor.SetVisibility(int(props['visible']))

    def vtk_calc_distance_from_geometry(self):

        source = self.collect_toplevel_geometry()

        if source:
            signed_distances = vtk.vtkFloatArray()
            signed_distances.SetNumberOfComponents(1)
            signed_distances.SetName("SignedDistances")

            implicit_poly_data_distance = vtk.vtkImplicitPolyDataDistance()

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
        if self.ui.mesh.checkbox_mesh_inside.isChecked():
            clipper.InsideOutOn()
        else:
            clipper.InsideOutOff()

        clipper.SetMergeTolerance(
            safe_float(self.ui.mesh.lineedit_vtk_mesh_merge.text())
            )
        clipper.SetValue(0.0)
        clipper.Update()
        self.export_mesh(clipper)

        self.mesh = clipper.GetOutput()
        self.show_mesh()
        self.mesh_stats()

    def vtk_mesher_threshold(self):

        self.vtk_calc_distance_from_geometry()

        thresh = vtk.vtkThreshold()
        thresh.SetInputData(self.rectilinear_grid)

        if self.ui.mesh.checkbox_threshold_inside.isChecked():
            thresh.ThresholdByLower(0)
        else:
            thresh.ThresholdByUpper(0)

        if self.ui.mesh.checkbox_threshold_interface.isChecked():
            thresh.AllScalarsOff()
        else:
            thresh.AllScalarsOn()

        thresh.Update()
        self.export_mesh(thresh)

        self.mesh = thresh.GetOutput()
        self.show_mesh()
        self.mesh_stats()

    def export_mesh(self, mesh, name=DEFAULT_MESH_NAME):
        # export geometry
        project_dir = self.parent.get_project_dir()
        self.export_unstructured(os.path.join(project_dir, name),
                                 mesh.GetOutputPort())

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

        self.ui.mesh.lineedit_mesh_cells.setText(str(cell_count))

        self.ui.mesh.plaintextedit_mesh_cell_types.setPlainText(
            '\n'.join([':\t'.join([str(cell_type), str(cells)])
                      for cell_type, cells in cell_type_counts])
            )

    def mesher(self):

        mesher = str(self.ui.mesh.combobox_mesher.currentText())

        if mesher == 'vtkTableBasedClipDataSet':
            self.vtk_mesher_table_based()
        elif mesher == 'vtkThreshold':
            self.vtk_mesher_threshold()

    def remove_mesh(self, name=DEFAULT_MESH_NAME):
        project_dir = self.parent.get_project_dir()
        path = os.path.join(project_dir, name)
        if not os.path.exists(path): return
        key = gui.message(text='Remove %s?' % name, buttons=['yes', 'no'], default='no')
        if key == 'no': return

        os.remove(path)

    def set_mesh_actor_props(self):
        props = self.visual_props['mesh']
        self.set_representation(self.mesh_actor, props['rep'])
        self.mesh_actor.GetProperty().SetColor(props['color'].getRgbF()[:3])
        self.mesh_actor.GetProperty().SetEdgeColor(props['edge'].getRgbF()[:3])
        self.mesh_actor.SetVisibility(int(props['visible']))

    def show_mesh(self):
        if self.mesh is None: return

        if not self.mesh_actor is None:
            self.vtkrenderer.RemoveActor(self.mesh_actor)
        self.mesh_mapper = vtk.vtkDataSetMapper()
        self.mesh_mapper.ScalarVisibilityOff()
        self.mesh_mapper.SetInputData(self.mesh)
        self.mesh_actor = vtk.vtkActor()
        self.mesh_actor.SetMapper(self.mesh_mapper)
        self.vtkrenderer.AddActor(self.mesh_actor)
        self.set_mesh_actor_props()
        self.render()

    # --- view ---
    def perspective(self, parallel=None):
        camera = self.vtkrenderer.GetActiveCamera()

        if parallel is None:
            parallel = not camera.GetParallelProjection()

        if parallel:
            camera.ParallelProjectionOn()
            self.toolbutton_perspective.setIcon(get_icon('parallel.png'))
        else:
            camera.ParallelProjectionOff()
            self.toolbutton_perspective.setIcon(get_icon('perspective.png'))

        self.render()

    def set_view(self, view='xy'):
        self.perspective(parallel=True)
        camera = self.vtkrenderer.GetActiveCamera()
        if view == 'xy':
            camera.SetPosition(0, 0, 10000000)
            camera.SetViewUp(0, 1, 0)
            if self.view_flip[1]:
                camera.Azimuth(180)
            self.view_flip[1] = not self.view_flip[1]
        elif view == 'yz':
            camera.SetPosition(0, 10000000, 0)
            camera.SetViewUp(1, 0, 0)
            if self.view_flip[2]:
                camera.Azimuth(180)
            self.view_flip[2] = not self.view_flip[2]
        elif view == 'xz':
            camera.SetPosition(10000000, 0, 0)
            camera.SetViewUp(0, 1, 0)
            if self.view_flip[2]:
                camera.Azimuth(180)
            self.view_flip[2] = not self.view_flip[2]

        self.reset_view()

    def reset_view(self):
        self.vtkrenderer.ResetCamera()
        self.render()

    def change_visibility(self, name, toolbutton):
        actors = None
        if name == 'mesh':
            actors = [self.mesh_actor]
        elif name == 'background_mesh':
            actors = self.grid_viewer_dict['actors']
        elif name == 'geometry':
            actors = [geo['actor'] for geo in self.geometrydict.values()
                      if geo['visible']]
        elif name == 'regions':
            actors = [geo['actor'] for geo in self.region_dict.values()
                      if geo['visible']]

            actors += [geo['clip_actor'] for geo in self.region_dict.values()
                       if 'clip_actor' in geo and geo['visible']]

        visible = self.visual_props[name]['visible'] = not toolbutton.isChecked()
        self.set_visible_btn_image(toolbutton, visible)
        if actors is not None:
            for actor in actors:
                if visible:
                    actor.VisibilityOn()
                else:
                    actor.VisibilityOff()
            self.render()

    def get_actors(self, name):
        actors = []
        if name == 'mesh':
            actors = [self.mesh_actor]
        elif name == 'background_mesh':
            actors = self.grid_viewer_dict['actors']
        elif name == 'geometry':
            actors = [geo['actor'] for geo in self.geometrydict.values()]
        elif name == 'regions':
            actors = [geo['actor'] for geo in self.region_dict.values()]
            actors += [geo['clip_actor'] for geo in self.region_dict.values()
                       if 'clip_actor' in geo]
        for actor in actors:
            yield actor

    def set_representation(self, actor, rep):
        if rep == 'wire':
            actor.GetProperty().SetRepresentationToWireframe()
        elif rep == 'solid':
            actor.GetProperty().SetRepresentationToSurface()
            actor.GetProperty().EdgeVisibilityOff()
        elif rep == 'edges':
            actor.GetProperty().SetRepresentationToSurface()
            actor.GetProperty().EdgeVisibilityOn()
        elif rep == 'points':
            actor.GetProperty().SetRepresentationToPoints()

    def change_representation(self, name, combobox):
        representation = str(combobox.currentText())

        self.visual_props[name]['rep'] = representation

        for actor in self.get_actors(name):
            self.set_representation(actor, representation)

        self.render()

    def change_color(self, name, button):
        col = QtWidgets.QColorDialog.getColor()
        if not col.isValid(): return

        button.setStyleSheet("QToolButton{{ background: {};}}".format(
            col.name()))

        self.visual_props[name]['color'] = col
        dark = self.visual_props[name]['edge'] = col.darker()

        for actor in self.get_actors(name):
            actor.GetProperty().SetColor(col.getRgbF()[:3])
            actor.GetProperty().SetEdgeColor(dark.getRgbF()[:3])

        self.render()

    def change_opacity(self, name, slider):
        value = slider.value()/100.0
        self.visual_props[name]['opacity'] = value

        for actor in self.get_actors(name):
            actor.GetProperty().SetOpacity(value)
        self.render()

    def set_visible_btn_image(self, btn, checked):
        if not checked:
            btn.setIcon(
                        get_icon('visibilityofftransparent.png'))
        else:
            btn.setIcon(get_icon('visibility.png'))
