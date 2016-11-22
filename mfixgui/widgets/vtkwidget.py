"""
This is the vtk widget. It handles displaying the 3D graphics in the GUI as
well as a simple geometry creation tool and region selection.
"""
from __future__ import print_function, absolute_import, unicode_literals, division
import os
import copy
import traceback
import logging
from functools import partial
import numpy as np
from qtpy import QtCore, QtGui, QtWidgets
LOG = logging.getLogger(__name__)

# VTK imports
VTK_AVAILABLE = True
try:
    import vtk
    LOG.info('VTK version: %s' % vtk.vtkVersion.GetVTKVersion())
    VTK_MAJOR_VERSION = vtk.VTK_MAJOR_VERSION
except ImportError:
    VTK_AVAILABLE = False
    e = traceback.format_exc()
    LOG.info("can't import vtk:\n{}".format(e))
try:
    # Try Qt 5.x
    from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
    
    # patch vtk 7
    # known bug with vtk 7
    # fixed with commit: https://gitlab.kitware.com/vtk/vtk/commit/90ee1cca513db11d3401cc25997c9a0b4ee15166
    if vtk.vtkVersion.GetVTKVersion() == '7.0.0':
        def wheelEvent(cls, ev):
            if hasattr(ev, 'delta'):
                cls.__wheelDelta += ev.delta()
            else:
                cls.__wheelDelta += ev.angleDelta().y()
            
            if cls.__wheelDelta >= 120:
                cls._Iren.MouseWheelForwardEvent()
                cls.__wheelDelta = 0
            elif cls.__wheelDelta <= -120:
                cls._Iren.MouseWheelBackwardEvent()
                cls.__wheelDelta = 0
            
        QVTKRenderWindowInteractor.__wheelDelta = 0
        QVTKRenderWindowInteractor.wheelEvent = wheelEvent
                
except ImportError:
    try:
        # Fall back to Qt 4.x
        from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
    except ImportError:
        VTK_AVAILABLE = False
        e = traceback.format_exc()
        LOG.info("Can't import QVTKRenderWindowInteractor:\n{}".format(e))

# local imports
from mfixgui.tools.general import (get_unique_string, widget_iter, get_icon,
                           get_image_path, topological_sort)
from mfixgui.widgets.base import LineEdit
try:
    from mfixgui.widgets.distributed_popup import DistributionPopUp
    from mfixgui.widgets.cyclone_popup import CyclonePopUp
    from mfixgui.widgets.reactor_popup import ReactorPopUp
    from mfixgui.widgets.hopper_popup import HopperPopUp
except ImportError:
    DistributionPopUp = None
    CyclonePopUp = None
    ReactorPopUp = None
    HopperPopUp = None
from mfixgui.project import Equation, ExtendedJSON

try:
    from mfixgui.widgets.vtk_constants import *
except ImportError:
    DEFAULT_MESH_NAME = 'geometry.stl'

GUI = None


def safe_float(value):
    """try to convert the value to a float, if ValueError, send error to gui
    and return 0"""
    try:
        return float(value)
    except ValueError as error:
        if GUI:
            GUI.error(str(error))
        else:
            LOG.error(str(error))
        return 0.0


def safe_int(value):
    """try to convert the value to a int, if ValueError, send error to gui
    and return 0"""
    try:
        return int(value)
    except ValueError as error:
        if GUI:
            GUI.error(str(error))
        else:
            LOG.error(error)
        return 0


def clean_geo_dict(dirty_dict):
    """clean the geometry dictionary so it can be saved"""
    clean_dict = {}
    for geo, geo_dict in dirty_dict.items():
        new_dict = clean_dict[geo] = {}
        geo_type = None
        if 'geo_type' in geo_dict:
            geo_type = geo_dict['geo_type']
        new_dict['geo_type'] = geo_type
        for key, value in geo_dict.items():
            # filter out vtk objects/default values
            if not isinstance(value, vtk.vtkObject) and (isinstance(value, Equation) or value != DEFAULT_PARAMS[geo_type][key]):
                new_dict[key] = value
    return clean_dict


def remove_vtk_objects(dirty_dict):
    """given a dictionary, remove vtkObjects from it"""
    clean_dict = {}
    for key, value in dirty_dict.items():
        if not isinstance(value, vtk.vtkObject):
            clean_dict[key] = value
    return clean_dict


def clean_visual_dict(dirty_dict):
    """remove qcolor objects from visual dict and save the rgb values"""
    clean_dict = {}
    for geo, geo_dict in dirty_dict.items():
        clean_geo = clean_dict[geo] = {}
        for key, value in geo_dict.items():
            if key in ['color', 'edge']:
                clean_geo[key] = geo_dict[key].name()
            else:
                clean_geo[key] = geo_dict[key]

    return clean_dict


def is_stl_ascii(fname):
    """see if the stl file is ASCII by checkign the first line for solid"""
    with open(fname, 'rb') as stlFile:
        try:
            solid = stlFile.readline().strip().lower().startswith(b'solid')
        except:
            solid = False
    return solid


def purge_multi_solids(fname):
    """Remove multiple solids from an stl file, only works with ascii"""
    # if it already ends in .onesolid.stl, assume cleaned
    if fname.endswith('.onesolid.stl'):
        return fname
    name = os.path.splitext(os.path.basename(fname))[0]
    dir = os.path.dirname(fname)
    newfile = os.path.join(dir, name + '.onesolid.stl')
    multi_solid = 0
    with open(fname, "r") as input:
        with open(newfile, "w") as output:
            output.write('solid ascii\n')
            for line in input:
                if 'solid' not in line:
                    output.write(line)
                else:
                    multi_solid += 1
            output.write('endsolid\n')
    if multi_solid > 2:
        LOG.warn('the stl file: %s has multiple solids, removing' % fname)
        return newfile
    else:
        os.remove(newfile)
        return fname


#class CustomInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
#    """custom vtkInteractorStyleTrackballCamera to highlight selected
#    objects"""
#    def __init__(self, parent=None):
#        self.AddObserver("LeftButtonPressEvent", self.left_button_press_event)
#
#        self.last_picked_actor = None
#        self.last_picked_property = vtk.vtkProperty()
#
#    def left_button_press_event(self, obj, event):
#        """on a left mouse press event, see if there is an actor and highlight
#        it"""
#        clickPos = self.GetInteractor().GetEventPosition()
#
#        picker = vtk.vtkPropPicker()
#        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())
#
#        # get the new actor
#        self.new_picked_actor = picker.GetActor()
#
#        # If something was selected
#        if self.new_picked_actor:
#            # If we picked something before, reset its property
#            if self.last_picked_actor:
#                self.last_picked_actor.GetProperty().DeepCopy(
#                    self.last_picked_property)
#
#            # Save the property of the picked actor so that we can
#            # restore it next time
#            self.last_picked_property.DeepCopy(self.new_picked_actor.GetProperty())
#            # Highlight the picked actor by changing its properties
#            self.new_picked_actor.GetProperty().SetColor(255/255.0, 140/255.0, 0)
#            self.new_picked_actor.GetProperty().SetDiffuse(1.0)
#            self.new_picked_actor.GetProperty().SetSpecular(0.0)
#
#            # save the last picked actor
#            self.last_picked_actor = self.new_picked_actor
#
#        # clear selection
#        elif self.last_picked_actor:
#            self.last_picked_actor.GetProperty().DeepCopy(
#                self.last_picked_property)
#            self.last_picked_actor = None
#
#        self.OnLeftButtonDown()
#        return


class VtkWidget(QtWidgets.QWidget):
    """the vtk widget"""
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, project, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        global GUI
        GUI = parent

        self.project = project
        self.parent = parent

        self.ui = parent.ui
        self.geometrytree = self.ui.geometry.treeWidgetGeometry

        # --- data ---
        self.animate = True
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

        #self.style = CustomInteractorStyle()
        self.style_3d = vtk.vtkInteractorStyleTrackballCamera()
        self.style_3d.SetDefaultRenderer(self.vtkrenderer)
        self.style_2d = vtk.vtkInteractorStyleImage()
        self.style_2d.SetDefaultRenderer(self.vtkrenderer)
        self.vtkiren.SetInteractorStyle(self.style_3d)

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
        self.__setup_wizards()

        self.ui.geometry.stackedWidgetGeometryDetails.setCurrentIndex(0)

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
            self.selected_geometry_changed)
        self.geometrytree.itemClicked.connect(self.geometry_clicked)

        # --- geometry button ---
        self.add_geometry_menu = QtWidgets.QMenu(self)
        self.ui.geometry.toolbutton_add_geometry.setMenu(
            self.add_geometry_menu)

        action = QtWidgets.QAction('STL File', self.add_geometry_menu)
        action.triggered.connect(self.add_stl)
        action.setIcon(get_icon('geometry.png'))
        action.setIconVisibleInMenu(True)
        self.add_geometry_menu.addAction(action)

        # --- implicit functions ---
        p_menu = self.add_geometry_implicit = QtWidgets.QMenu(self)
        p_menu.setTitle('Implicits')
        p_menu.setIcon(get_icon('function.png'))
        a = self.add_geometry_menu.addMenu(p_menu)
        a.setIconVisibleInMenu(True)
        for geo in IMPLICIT_DICT.keys():
            action = QtWidgets.QAction(geo.replace('_', ' '), p_menu)
            action.triggered.connect(partial(self.add_implicit, implicittype=geo))
            action.setIcon(get_icon('function.png'))
            action.setIconVisibleInMenu(True)
            p_menu.addAction(action)

        # --- primitives ---
        p_menu = self.add_geometry_primitive = QtWidgets.QMenu(self)
        p_menu.setTitle('Primitives')
        p_menu.setIcon(get_icon('geometry.png'))
        a = self.add_geometry_menu.addMenu(p_menu)
        a.setIconVisibleInMenu(True)
        for geo in PRIMITIVE_DICT.keys():
            action = QtWidgets.QAction(geo, p_menu)
            action.triggered.connect(partial(self.add_primitive, primtype=geo))
            action.setIcon(get_icon('geometry.png'))
            action.setIconVisibleInMenu(True)
            p_menu.addAction(action)

        # --- parametric ---
        p_menu = self.add_geometry_parametric = QtWidgets.QMenu(self)
        p_menu.setTitle('Parametrics')
        p_menu.setIcon(get_icon('geometry.png'))
        a = self.add_geometry_menu.addMenu(p_menu)
        a.setIconVisibleInMenu(True)
        for geo in PARAMETRIC_DICT.keys():
            action = QtWidgets.QAction(geo.replace('_', ' '), p_menu)
            action.triggered.connect(partial(self.add_parametric, paramtype=geo))
            action.setIcon(get_icon('geometry.png'))
            action.setIconVisibleInMenu(True)
            p_menu.addAction(action)


        # --- filter button ---
        self.add_filter_menu = QtWidgets.QMenu(self)
        self.ui.geometry.toolbutton_add_filter.setMenu(self.add_filter_menu)
        self.add_filter_menu.aboutToShow.connect(self.enable_filters)
        self.filter_actions = []
        for geo in FILTER_DICT.keys():
            action = QtWidgets.QAction(geo.replace('_', ' '), self.add_filter_menu)
            action.triggered.connect(partial(self.add_filter, filtertype=geo))
            action.setIcon(get_icon('filter.png'))
            action.setIconVisibleInMenu(True)
            self.filter_actions.append(action)
            self.add_filter_menu.addAction(action)

        # setup signals
        self.ui.geometry.toolbutton_remove_geometry.pressed.connect(
            self.remove_geometry)
        self.ui.geometry.toolbutton_copy_geometry.pressed.connect(
            self.handle_copy_geometry)

        # connect boolean
        for key, btn in self.booleanbtndict.items():
            btn.pressed.connect(partial(self.boolean_operation, booltype=key))

        # connect parameter widgets
        for widget in widget_iter(
                self.ui.geometry.stackedWidgetGeometryDetails):
            if isinstance(widget, QtWidgets.QLineEdit):
                widget.editingFinished.connect(partial(self.parameter_edited, widget))
            elif isinstance(widget, QtWidgets.QCheckBox):
                widget.stateChanged.connect(lambda i, w=widget: self.parameter_edited(w))

        # --- mesh ---
        self.ui.geometry.pushbutton_mesh_autosize.pressed.connect(
            self.auto_size_mesh_extents)

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
        self.toolbutton_view_xy.setIcon(get_icon('xy.png'))

        self.toolbutton_view_yz = QtWidgets.QToolButton()
        self.toolbutton_view_yz.pressed.connect(lambda: self.set_view('yz'))
        self.toolbutton_view_yz.setIcon(get_icon('yz.png'))

        self.toolbutton_view_xz = QtWidgets.QToolButton()
        self.toolbutton_view_xz.pressed.connect(lambda: self.set_view('xz'))
        self.toolbutton_view_xz.setIcon(get_icon('xz.png'))

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
            btns = self.visual_btns[geo] = {}
            # tool button
            toolbutton = QtWidgets.QToolButton()
            toolbutton.pressed.connect(partial(self.change_visibility, geo, toolbutton))
            toolbutton.setCheckable(True)
            toolbutton.setChecked(True)
            toolbutton.setAutoRaise(True)
            toolbutton.setIcon(get_icon('visibility.png'))
            layout.addWidget(toolbutton, i, 0)
            btns['visible'] = toolbutton

            # style
            combobox = QtWidgets.QComboBox()
            combobox.addItems(['wire', 'solid', 'edges', 'points'])
            combobox.activated.connect(partial(self.change_representation, geo, combobox))
            layout.addWidget(combobox, i, 1)
            btns['rep'] = combobox

            # color
            if not geo == 'regions':
                toolbutton = QtWidgets.QToolButton()
                toolbutton.pressed.connect(partial(self.change_color, geo, toolbutton))
                toolbutton.setAutoRaise(True)
                layout.addWidget(toolbutton, i, 2)
                btns['color'] = toolbutton

            # opacity
            opacity = QtWidgets.QDoubleSpinBox()
            opacity.setRange(0, 1)
            opacity.setSingleStep(0.1)
            opacity.valueChanged.connect(partial(self.change_opacity, geo, opacity))
            layout.addWidget(opacity, i, 3)
            btns['opacity'] = opacity

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

    def __setup_wizards(self):

        self.cyclone_popup = None
        self.distribution_popup = None
        self.reactor_popup = None
        self.hopper_popup = None

        # --- wizards ---
        wizard_menu = QtWidgets.QMenu('wizards', self)

        action = QtWidgets.QAction('distributed', wizard_menu)
        action.triggered.connect(self.handle_distributed_wizard)
        wizard_menu.addAction(action)

        action = QtWidgets.QAction('cyclone', wizard_menu)
        action.triggered.connect(self.handle_cyclone_wizard)
        wizard_menu.addAction(action)

        action = QtWidgets.QAction('reactor', wizard_menu)
        action.triggered.connect(self.handle_reactor_wizard)
        wizard_menu.addAction(action)

        action = QtWidgets.QAction('hopper', wizard_menu)
        action.triggered.connect(self.handle_hopper_wizard)
        wizard_menu.addAction(action)

        self.ui.geometry.toolbutton_wizard.setMenu(wizard_menu)

    def set_visual_btn_values(self):
        """change the visual btns to be in sync with self.visual_props"""
        for geo, info in self.visual_props.items():
            for key, value in info.items():
                if geo in self.visual_btns and key in self.visual_btns[geo]:
                    wid = self.visual_btns[geo][key]
                    if key == 'rep':
                        wid.setCurrentIndex(wid.findText(value))
                    elif key == 'visible':
                        wid.setChecked(value)
                        self.set_visible_btn_image(wid, value)
                    elif key == 'color':
                        wid.setStyleSheet("QToolButton{{ background: {};}}".format(value.name()))
                    elif key == 'opacity':
                        wid.setValue(value)

    def emitUpdatedValue(self, key, value, args=None):
        """emit an updates value"""
        self.value_updated.emit(self, {key: value}, args)

    def updateValue(self, key, newValue, args=None):
        """receive keyword changed from project manager"""

        if key == 'no_k':
            self.change_interaction(newValue)
        elif key == 'out_stl_value':
            self.ui.mesh.checkbox_internal_external_flow.setChecked(newValue == 1.0)

    def objectName(self):
        """return the name of this object"""
        return 'VTK Widget'

    def default(self):
        """reset to defaults"""
        self.ui.geometry.lineedit_keyword_zlength.setEnabled(True)
        self.ui.mesh.lineedit_keyword_kmax.setEnabled(True)
        self.vtkrenderer.RemoveAllViewProps()
        self.clear_all_geometry()
        self.change_interaction()
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
        data = ExtendedJSON.loads(string)
        tree = data.get('tree')
        geo_dict = data.get('geometry_dict')

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
                    geo_type = geo_dict[node]['geo_type']
                    if geo_type == 'primitive':
                        self.add_primitive(name=node, data=geo_data, loading=True)
                    elif geo_type == 'implicit':
                        self.add_implicit(name=node, data=geo_data, loading=True)
                    elif geo_type == 'parametric':
                        self.add_parametric(name=node, data=geo_data, loading=True)
                    elif geo_type == 'filter':
                        self.add_filter(name=node, data=geo_data,
                                        child=tree[node].pop(), loading=True)
                    elif geo_type == 'boolean' or geo_type == 'boolean_implicit':
                        self.boolean_operation(boolname=node, data=geo_data,
                                               children=tree[node], loading=True)
                    elif geo_type == 'stl':
                        self.add_stl(None, filename=geo_data['filename'],
                                     name=node, data=geo_data, loading=True)
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
                    if isinstance(value, (list, tuple)):
                        data[geo][key] = QtGui.QColor(*value)
                    else:
                        data[geo][key] = QtGui.QColor(value)

        self.visual_props = copy.deepcopy(DEFAULT_VISUAL_PROPS)
        self.visual_props.update(data)
        if self.mesh_actor is not None:
            self.set_mesh_actor_props()

        self.set_visual_btn_values()
        for actor in self.grid_viewer_dict['actors']:
            self.set_background_mesh_actor_props(actor)

    # --- render ---
    def render(self, force_render=False, defer_render=None):
        """render the vtk scene, checks for defer_render"""
        if defer_render is not None:
            self.defer_render = defer_render
        if not self.defer_render or force_render:
            self.vtkRenderWindow.Render()

    # --- geometry ---
    def selected_geometry_changed(self):
        """The selected geometry changed, update UI"""
        current_selection = self.geometrytree.selectedItems()
        top_level_items = [self.geometrytree.indexOfTopLevelItem(select) > -1
                           for select in current_selection]

        implicits = ['implicit' in self.geometrydict.get(select.text(0)).get('geo_type') for select in current_selection]

        # boolean btns
        enableboolbtn = False
        if len(current_selection) == 2 and all(top_level_items) and not any(implicits):
            enableboolbtn = True
        elif len(current_selection) >= 2 and all(top_level_items) and all(implicits):
            enableboolbtn = True
        for btn in self.booleanbtndict.values():
            btn.setEnabled(enableboolbtn)

        # enable/disable delete/copy/filter button
        enables = [False]*3
        if len(current_selection) == 1:
            if all(top_level_items):
                enables = [True]*3
            else:
                enables[2] = True
        elif any(top_level_items):
            enables[0] = True

        for enable, widget in zip(enables, [self.ui.geometry.toolbutton_remove_geometry,
                                            self.ui.geometry.toolbutton_add_filter,
                                            self.ui.geometry.toolbutton_copy_geometry,]):
            widget.setEnabled(enable)

        if current_selection:
            text = str(current_selection[-1].text(0)).lower()
            data = self.geometrydict.get(text)
            type_ = data['type']
            geo_type = data.get('geo_type', None)
            if geo_type == 'implicit':
                type_ += '_implicit'

            new_index = 0
            for i in range(
                    self.ui.geometry.stackedWidgetGeometryDetails.count()):
                widget = self.ui.geometry.stackedWidgetGeometryDetails.widget(i)
                if str(widget.objectName()) == type_:
                    new_index = i
                    break

            # set the widget values
            for child in widget_iter(widget):
                name = str(child.objectName()).lower().replace('_', '')
                for key, value in data.items():
                    if key in name:
                        break

                if isinstance(child, LineEdit):
                    child.updateValue(None, value)
                elif isinstance(child, QtWidgets.QCheckBox):
                    child.setChecked(value)

            self.ui.geometry.groupBoxGeometryParameters.setTitle(text)

        else:
            new_index = 0
            self.ui.geometry.groupBoxGeometryParameters.setTitle('Parameters')
            self.ui.geometry.toolbutton_remove_geometry.setEnabled(False)


        current_index = self.ui.geometry.stackedWidgetGeometryDetails.currentIndex()

        if self.animate:
            self.parent.animate_stacked_widget(
                self.ui.geometry.stackedWidgetGeometryDetails,
                current_index,
                new_index,
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
        """Hide/Show the clicked geometry"""

        name = str(item.text(0)).lower()
        geo = self.geometrydict.get(name)
        if geo is None:
            return
        if item.checkState(0) == QtCore.Qt.Checked:
            if self.visual_props['geometry']['visible']:
                geo['actor'].VisibilityOn()
            geo['visible'] = True
        else:
            geo['actor'].VisibilityOff()
            geo['visible'] = False

        self.render()

    def get_stl_extents(self, filename):
        '''given an stl file, return the extents of that file'''

        # purge solids
        if is_stl_ascii(filename):
            filename = purge_multi_solids(filename)

        # reader
        reader = vtk.vtkSTLReader()
        reader.SetFileName(filename)
        reader.MergingOn()
        reader.Update()

        return reader.GetOutput().GetBounds()

    def add_stl(self, widget, filename=None, name=None, data=None, loading=False):
        """Open browse dialog and load selected stl file"""

        if filename is None:
            filename = QtWidgets.QFileDialog.getOpenFileName(
                self, 'Select an STL File',
                self.parent.get_project_dir(),
                'STL File (*.stl)',)

            if isinstance(filename, (tuple, list)):
                filename = filename[0]

            filename = str(filename)

        if filename:
            # purge solids
            if is_stl_ascii(filename):
                filename = purge_multi_solids(filename)

            if name is None:
                name = os.path.basename(filename).lower()
                name = get_unique_string(name, list(self.geometrydict.keys()))

            if data is not None:
                self.geometrydict[name] = data
            else:
                self.geometrydict[name] = copy.deepcopy(DEFAULT_STL_PARAMS)

            geo_data = self.geometrydict.get(name)
            geo_data['filename'] = filename

            # reader
            reader = vtk.vtkSTLReader()
            reader.SetFileName(filename)
            reader.MergingOn()

            # Create transformer
            transform = vtk.vtkTransform()
            transform_filter = vtk.vtkTransformPolyDataFilter()
            transform_filter.SetTransform(transform)
            transform_filter.SetInputConnection(reader.GetOutputPort())
            transform_filter.Update()

            bounds = transform_filter.GetOutput().GetBounds()
            for key, bound in zip(['extentxmin', 'extentxmax', 'extentymin', 'extentymax', 'extentzmin', 'extentzmax'],
                                  bounds):
                geo_data[key] = bound

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
            geo_data['reader'] = reader
            geo_data['transform'] = transform
            geo_data['transformfilter'] = transform_filter
            geo_data['center_filter'] = center_filter
            geo_data['mapper'] = mapper
            geo_data['actor'] = actor
            geo_data['centerx'] = center[0]
            geo_data['centery'] = center[1]
            geo_data['centerz'] = center[2]

            # Add to tree
            item = QtWidgets.QTreeWidgetItem([name])
            item.setIcon(0, get_icon('geometry.png'))
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(0, QtCore.Qt.Checked)
            self.geometrytree.addTopLevelItem(item)
            self.geometrytree.setCurrentItem(item)

            if not loading:
                self.parent.set_unsaved_flag()

            return name

    def parameter_edited(self, widget, name=None, value=None, key=None):
        """Update the value of edited parameter in the geometrydict"""
        if name is None:
            current_selection = self.geometrytree.selectedItems()

            if current_selection:
                name = str(current_selection[-1].text(0)).lower()
            value = None

            parameters = str(widget.objectName()).lower().split('_')
            parameter = ''.join(parameters[1:])

            geo_data = self.geometrydict.get(name)

            for key in geo_data.keys():
                if key in parameter:
                    break

            # if widget is a lineedit
            if isinstance(widget, LineEdit):
                value = widget.value
            elif isinstance(widget, QtWidgets.QCheckBox):
                value = widget.isChecked()
        else:
            geo_data = self.geometrydict.get(name)

        if value is not None:
            self.update_parameter_map(value, name, key)
            geo_data[key] = value
            type_ = geo_data['type']
            geo_type = geo_data.get('geo_type')

            if type_ in PRIMITIVE_DICT and geo_type == 'primitive':
                self.update_primitive(name)
            elif type_ in FILTER_DICT:
                self.update_filter(name)
            elif type_ in PARAMETRIC_DICT:
                self.update_parametric(name)
            elif type_ in IMPLICIT_DICT:
                self.update_implicit(name)

            if 'transform' in geo_data and not geo_type in ['filter', 'implicit']:
                self.update_transform(name)
            self.render()

    def update_primitive(self, name):
        """Update the specified primitive"""
        primtype = self.geometrydict[name]['type']
        geo = self.geometrydict.get(name)

        if 'source' in self.geometrydict[name]:
            source = geo['source']
        else:
            source = None

        # update source
        if primtype == 'sphere':
            source.SetRadius(safe_float(geo['radius']))
            source.SetThetaResolution(safe_int(geo['thetaresolution']))
            source.SetPhiResolution(safe_int(geo['phiresolution']))
        elif primtype == 'box':
            source.SetXLength(safe_float(geo['lengthx']))
            source.SetYLength(safe_float(geo['lengthy']))
            source.SetZLength(safe_float(geo['lengthz']))
        elif primtype == 'cone':
            source.SetRadius(safe_float(geo['radius']))
            source.SetHeight(safe_float(geo['height']))
            source.SetDirection(safe_float(geo['directionx']),
                                safe_float(geo['directiony']),
                                safe_float(geo['directionz']))
            source.SetResolution(safe_int(geo['resolution']))
            source.CappingOn()
        elif primtype == 'cylinder':
            source.SetRadius(safe_float(geo['radius']))
            source.SetHeight(safe_float(geo['height']))
            source.SetResolution(safe_int(geo['resolution']))
        elif primtype == 'stl':
            pass
        else:
            return

        # common props
        if source is not None:
            source.SetCenter(safe_float(geo['centerx']),
                             safe_float(geo['centery']),
                             safe_float(geo['centerz']))
            source.Update()

        return source

    def update_transform(self, name):
        """Update the specified object's transform filter."""
        geo = self.geometrydict.get(name)
        transform = geo['transform']
        transform_filter = geo['transformfilter']

        # reset to Identity
        transform.Identity()
        transform.PostMultiply()
        transform.PreMultiply()

        # translate to center
        transform.Translate(-safe_float(geo['centerx']),
                            -safe_float(geo['centery']),
                            -safe_float(geo['centerz']))

        # rotation
        transform.RotateWXYZ(safe_float(geo['rotationx']), 1, 0, 0)
        transform.RotateWXYZ(safe_float(geo['rotationy']), 0, 1, 0)
        transform.RotateWXYZ(safe_float(geo['rotationz']), 0, 0, 1)

        # back to position
        transform.Translate(safe_float(geo['centerx']),
                            safe_float(geo['centery']),
                            safe_float(geo['centerz']))

        # translate stl files
        if geo['type'] in ['stl'] + list(PARAMETRIC_DICT.keys()):
            transform.Translate(
                safe_float(geo['translationx']),
                safe_float(geo['translationy']),
                safe_float(geo['translationz']),
                )

        # update
        transform_filter.Update()

        # update stl extents
        if geo['type'] in ['stl']:
            bounds = transform_filter.GetOutput().GetBounds()
            for key, bound in zip(['extentxmin', 'extentxmax', 'extentymin', 'extentymax', 'extentzmin', 'extentzmax'],
                                  bounds):
                geo[key] = bound

            self.selected_geometry_changed()

        return transform_filter

    def add_primitive(self, primtype='sphere', name=None, data=None, loading=False):
        """Add the specified primitive"""
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
        geo = self.geometrydict.get(name)
        geo['source'] = source

        source = self.update_primitive(name)

        # Create transformer
        transform = vtk.vtkTransform()
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputConnection(source.GetOutputPort())

        geo['transform'] = transform
        geo['transformfilter'] = transform_filter

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
        geo['mapper'] = mapper
        geo['actor'] = actor
        geo['trianglefilter'] = trianglefilter
        geo['source'] = source

        # Add to tree
        item = QtWidgets.QTreeWidgetItem([name])
        item.setIcon(0, get_icon('geometry.png'))
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(0, QtCore.Qt.Checked)
        self.geometrytree.addTopLevelItem(item)
        self.geometrytree.setCurrentItem(item)

        if not loading:
            self.parent.set_unsaved_flag()

        return name

    def update_implicit(self, name):
        """update the implicit function"""
        implicittype = self.geometrydict[name]['type']
        geo = self.geometrydict.get(name)

        source = geo.get('source', None)
        sample = geo.get('sample', None)
        surface = geo.get('surface', None)
        transform = geo.get('transform', None)

        x, y, z = safe_float(geo['centerx']), safe_float(geo['centery']), safe_float(geo['centerz'])
        rotx, roty, rotz = safe_float(geo['rotationx']), safe_float(geo['rotationy']), safe_float(geo['rotationz'])
        r = safe_float(geo['radius'])
        h = safe_float(geo['height'])

        # update transform
        if transform:
            # reset to Identity
            transform.Identity()
            transform.PostMultiply()

            # back to position
            transform.Translate(-x, -y, -z)

            # rotation
            transform.RotateWXYZ(rotx, 1, 0, 0)
            transform.RotateWXYZ(roty, 0, 1, 0)
            transform.RotateWXYZ(rotz, 0, 0, 1)

        # update source
        bounds = [safe_float(geo[k]) for k in ['minx','maxx','miny', 'maxy', 'minz', 'maxz']]
        if implicittype == 'sphere':
            source.SetRadius(r)
            bounds = [-r, r, -r, r, -r, r]
        elif implicittype == 'box':
            dx = safe_float(geo['lengthx'])/2.0
            dy = safe_float(geo['lengthy'])/2.0
            dz = safe_float(geo['lengthz'])/2.0
            bounds = [-dx, dx, -dy, dy, -dz, dz]
            source.SetBounds(bounds)
        elif implicittype == 'cone':
            angle = np.rad2deg(np.arctan(r/h))
            geo['cone_source'].SetAngle(angle)
            geo['plane1'].SetOrigin(h, 0, 0)
            geo['plane1'].SetNormal(-1, 0, 0)
            geo['plane2'].SetOrigin(0, 0, 0)
            geo['plane2'].SetNormal(1, 0, 0)
            transform.Translate(h/2.0, 0, 0)
            bounds = [0, h, -r, r, -r, r]
        elif implicittype == 'cylinder':
            geo['cylinder_source'].SetRadius(r)
            geo['plane1'].SetOrigin(x, y + h/2, z)
            geo['plane1'].SetNormal(0, -1, 0)
            geo['plane2'].SetOrigin(x, y - h/2, z)
            geo['plane2'].SetNormal(0, 1, 0)
            bounds = [-r, r, -h/2, h/2, -r, r]
        elif implicittype == 'quadric':
            source.SetCoefficients(
                safe_float(geo['a0']),
                safe_float(geo['a1']),
                safe_float(geo['a2']),
                safe_float(geo['a3']),
                safe_float(geo['a4']),
                safe_float(geo['a5']),
                safe_float(geo['a6']),
                safe_float(geo['a7']),
                safe_float(geo['a8']),
                safe_float(geo['a9']),
                )
        elif implicittype == 'superquadric':
            source.SetPhiRoundness(safe_float(geo['phi']))
            source.SetThetaRoundness(safe_float(geo['theta']))
            source.SetThickness(safe_float(geo['thickness']))
            source.SetSize(r)
            source.SetToroidal(safe_int(geo['toroidal']))
            bounds = [-r, r, -r, r, -r, r]
        else:
            return

        # common props
        if source is not None and hasattr(source, 'Update'):
            source.Update()

        if sample:
            # transform bounds
            t_mat = transform.GetInverse().GetMatrix()
            bounds_list = []
            for dx in bounds[:2]:
                for dy in bounds[2:4]:
                    for dz in bounds[4:]:
                        bounds_list.append(t_mat.MultiplyFloatPoint([dx, dy, dz, 1]))

            xs = [i[0] for i in bounds_list]
            ys = [i[1] for i in bounds_list]
            zs = [i[2] for i in bounds_list]

            bounds = [min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)]
            geo['bounds'] = bounds
            sample.SetModelBounds(*bounds)
            sample.SetSampleDimensions(IMPLICIT_DEFAULT_RES, IMPLICIT_DEFAULT_RES, IMPLICIT_DEFAULT_RES)
            sample.Update()

        if surface:
            surface.Update()

    def add_implicit(self, implicittype=None, name=None, data=None, loading=False):
        """Add an implicit function"""

        # create implicit
        if data is not None:
            implicittype = data['type']
        if implicittype in IMPLICIT_DICT:
            source = IMPLICIT_DICT[implicittype]()
        else:
            return

        if name is None:
            name = get_unique_string(implicittype, list(self.geometrydict.keys()))

        if data is None:
            self.geometrydict[name] = copy.deepcopy(DEFAULT_IMPLICIT_PARAMS)
            self.geometrydict[name]['type'] = implicittype
        else:
            self.geometrydict[name] = data
        geo = self.geometrydict.get(name)

        # transform
        transform = vtk.vtkTransform()
        source.SetTransform(transform)

        if implicittype == 'cylinder':
            boolean = vtk.vtkImplicitBoolean()
            boolean.SetOperationTypeToDifference()
            boolean.AddFunction(source)
            geo['cylinder_source'] = source
            p1 = geo['plane1'] = vtk.vtkPlane()
            p1.SetTransform(transform)
            boolean.AddFunction(p1)
            p2 = geo['plane2'] = vtk.vtkPlane()
            p2.SetTransform(transform)
            boolean.AddFunction(p2)
            source = boolean
        elif implicittype == 'cone':
            boolean = vtk.vtkImplicitBoolean()
            boolean.SetOperationTypeToDifference()
            boolean.AddFunction(source)
            geo['cone_source'] = source
            p1 = geo['plane1'] = vtk.vtkPlane()
            p1.SetTransform(transform)
            boolean.AddFunction(p1)
            p2 = geo['plane2'] = vtk.vtkPlane()
            p2.SetTransform(transform)
            boolean.AddFunction(p2)
            source = boolean

        sample = vtk.vtkSampleFunction()
        sample.SetSampleDimensions(IMPLICIT_DEFAULT_RES, IMPLICIT_DEFAULT_RES, IMPLICIT_DEFAULT_RES)
        sample.SetImplicitFunction(source)
        sample.ComputeNormalsOff()

        # contour
        surface = vtk.vtkContourFilter()
        surface.SetInputConnection(sample.GetOutputPort())
        surface.SetValue(0, 0.0)

        geo['source'] = source
        geo['sample'] = sample
        geo['surface'] = surface
        geo['transform'] = transform

        self.update_implicit(name)

        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(surface.GetOutputPort())
        mapper.ScalarVisibilityOff()

        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        self.set_geometry_actor_props(actor, name)
        self.vtkrenderer.AddActor(actor)
        self.render()

        geo['mapper'] = mapper
        geo['actor'] = actor

        # Add to tree
        item = QtWidgets.QTreeWidgetItem([name])
        item.setIcon(0, get_icon('function.png'))
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(0, QtCore.Qt.Checked)
        self.geometrytree.addTopLevelItem(item)
        self.geometrytree.setCurrentItem(item)

        if not loading:
            self.parent.set_unsaved_flag()

    def update_parametric(self, name):
        """Update the specified parameteric object."""
        geo = self.geometrydict.get(name)
        paratype = geo['type']
        para_object = geo['parametric_object']
        source = geo['source']

        if paratype == 'torus':
            para_object.SetRingRadius(safe_float(geo['ringradius']))
            para_object.SetCrossSectionRadius(safe_float(geo['crosssectionradius']))
        elif paratype == 'boy':
            para_object.SetZScale(safe_float(geo['zscale']))
        elif paratype == 'conic_spiral':
            para_object.SetA(safe_float(geo['ascale']))
            para_object.SetB(safe_float(geo['bfunc']))
            para_object.SetC(safe_float(geo['cfunc']))
            para_object.SetN(safe_float(geo['nfunc']))
        elif paratype == 'dini':
            para_object.SetA(safe_float(geo['ascale']))
            para_object.SetB(safe_float(geo['bscale']))
        elif paratype == 'ellipsoid':
            para_object.SetXRadius(safe_float(geo['radiusx']))
            para_object.SetYRadius(safe_float(geo['radiusy']))
            para_object.SetZRadius(safe_float(geo['radiusz']))
        elif paratype == 'figure_8_klein':
            para_object.SetRadius(safe_float(geo['radius']))
        elif paratype == 'mobius':
            para_object.SetRadius(safe_float(geo['radius']))
        elif paratype == 'random_hills':
            para_object.SetHillXVariance(safe_float(geo['variancex']))
            para_object.SetXVarianceScaleFactor(safe_float(geo['scalex']))
            para_object.SetHillYVariance(safe_float(geo['variancey']))
            para_object.SetYVarianceScaleFactor(safe_float(geo['scaley']))
            para_object.SetHillAmplitude(safe_float(geo['amplitude']))
            para_object.SetAmplitudeScaleFactor(safe_float(geo['scaleamplitude']))
            para_object.SetNumberOfHills(safe_int(geo['nhills']))
            para_object.SetAllowRandomGeneration(safe_int(geo['allowrandom']))
        elif paratype == 'roman':
            para_object.SetRadius(safe_float(geo['radius']))
        elif paratype == 'super_ellipsoid':
            para_object.SetXRadius(safe_float(geo['radiusx']))
            para_object.SetYRadius(safe_float(geo['radiusy']))
            para_object.SetZRadius(safe_float(geo['radiusz']))
            para_object.SetN1(safe_float(geo['n1']))
            para_object.SetN2(safe_float(geo['n2']))
        elif paratype == 'super_toroid':
            para_object.SetXRadius(safe_float(geo['radiusx']))
            para_object.SetYRadius(safe_float(geo['radiusy']))
            para_object.SetZRadius(safe_float(geo['radiusz']))
            para_object.SetRingRadius(safe_float(geo['ringradius']))
            para_object.SetCrossSectionRadius(safe_float(geo['crosssectionradius']))
            para_object.SetN1(safe_float(geo['n1']))
            para_object.SetN2(safe_float(geo['n2']))

        source.Update()
        self.parent.set_unsaved_flag()
        return source

    def add_parametric(self, paramtype=None, name=None, data=None, loading=False):
        """Add the specified parametric object"""

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
        geo = self.geometrydict.get(name)
        geo['parametric_object'] = parametric_object
        geo['source'] = source

        source = self.update_parametric(name)

        # Create transformer
        transform = vtk.vtkTransform()
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputConnection(source.GetOutputPort())

        geo['transform'] = transform
        geo['transformfilter'] = transform_filter

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
        geo['mapper'] = mapper
        geo['actor'] = actor
        geo['trianglefilter'] = trianglefilter
        geo['source'] = source

        # Add to tree
        item = QtWidgets.QTreeWidgetItem([name])
        item.setIcon(0, get_icon('geometry.png'))
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(0, QtCore.Qt.Checked)
        self.geometrytree.addTopLevelItem(item)
        self.geometrytree.setCurrentItem(item)

        if not loading:
            self.parent.set_unsaved_flag()

        return name

    def boolean_operation(self, booltype=None, boolname=None, data=None, children=None, loading=False):
        """Apply a boolean operation with the currently selected toplevel
        items."""

        if children is not None:
            current_selection = []
            for child in children:
                current_selection.append(self.get_tree_item(child))
        else:
            current_selection = self.geometrytree.selectedItems()

        if boolname is None:
            boolname = get_unique_string(booltype, list(self.geometrydict.keys()))

        # Save references
        if data is not None:
            bool_data = self.geometrydict[boolname] = data
            booltype = data['type']
        else:
            bool_data = self.geometrydict[boolname] = copy.deepcopy(DEFAULT_BOOLEAN_PARAMS)
            bool_data['type'] = booltype

        implicit = all(['implicit' in self.geometrydict.get(select.text(0)).get('geo_type') for select in current_selection])

        if implicit:
            boolean_operation = vtk.vtkImplicitBoolean()
        else:
            boolean_operation = vtk.vtkBooleanOperationPolyDataFilter()

        if booltype == 'union':
            if implicit:
                boolean_operation.SetOperationTypeToUnion()
            else:
                boolean_operation.SetOperationToUnion()
            icon = 'union'
        elif booltype == 'intersection':
            if implicit:
                boolean_operation.SetOperationTypeToIntersection()
            else:
                boolean_operation.SetOperationToIntersection()
            icon = 'intersect'
        else:
            if implicit:
                boolean_operation.SetOperationTypeToDifference()
            else:
                boolean_operation.SetOperationToDifference()
            icon = 'difference'

        union_list = []
        for i, selection in enumerate(current_selection):
            name = str(selection.text(0)).lower()
            bool_data['children'].append(name)
            child = self.geometrydict[name]
            if implicit:
                boolean_operation.AddFunction(child['source'])
                union_list.append(child['bounds'])
            else:
                geometry = self.get_input_data(name)
                boolean_operation.SetInputConnection(
                    i, geometry.GetOutputPort())

            # hide the sources
            child['actor'].VisibilityOff()
            child['visible'] = False

        mapper = vtk.vtkPolyDataMapper()

        if implicit:
            sample = vtk.vtkSampleFunction()
            sample.SetSampleDimensions(IMPLICIT_DEFAULT_RES, IMPLICIT_DEFAULT_RES, IMPLICIT_DEFAULT_RES)
            sample.SetImplicitFunction(boolean_operation)
            sample.ComputeNormalsOff()

            # union the bounds
            array = np.asarray(union_list)
            bounds = []
            for i in range(0, 6, 2):
                bounds += [min(array[:,i]), max(array[:,i+1])]
            sample.SetModelBounds(*bounds)

            # contour
            surface = vtk.vtkContourFilter()
            surface.SetInputConnection(sample.GetOutputPort())
            surface.SetValue(0, 0.0)

            bool_data['sample'] = sample
            bool_data['surface'] = surface
            bool_data['geo_type'] = 'boolean_implicit'
            bool_data['source'] = boolean_operation
            bool_data['bounds'] = bounds

            mapper.SetInputConnection(surface.GetOutputPort())
        else:
            boolean_operation.Update()
            mapper.SetInputConnection(boolean_operation.GetOutputPort())
        mapper.ScalarVisibilityOff()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        self.set_geometry_actor_props(actor, boolname)

        self.vtkrenderer.AddActor(actor)

        self.render()

        # save references
        bool_data['booleanoperation'] = boolean_operation
        bool_data['mapper'] = mapper
        bool_data['actor'] = actor

        # Add to tree
        toplevel = QtWidgets.QTreeWidgetItem([boolname])
        toplevel.setIcon(0, get_icon(icon+'.png'))
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

        if not loading:
            self.parent.set_unsaved_flag()

        return boolname

    def clear_all_geometry(self):
        """remove all geometry"""
        self.geometrytree.clear()
        for name, geo in self.geometrydict.items():
            self.remove_from_parameter_map(name, geo)
        self.geometrydict = {}

    def remove_geometry(self):
        """Remove the currently selected geometry, filter, or boolean operation
        """
        currentSelection = self.geometrytree.selectedItems()
        self.animate = False
        for selection in currentSelection:
            text = str(selection.text(0)).lower()
            # remove tree item
            toplevelindex = self.geometrytree.indexOfTopLevelItem(
                selection)
            item = self.geometrytree.takeTopLevelItem(toplevelindex)

            # move children to toplevel, make visible
            children = item.takeChildren()
            for child in children:
                self.geometrytree.addTopLevelItem(child)
                if self.visual_props['geometry']['visible']:
                    geo = self.geometrydict.get(str(child.text(0)).lower())
                    geo['actor'].VisibilityOn()
                    geo['visible'] = True
                child.setCheckState(0, QtCore.Qt.Checked)

            # remove graphics
            geo = self.geometrydict.pop(text)
            self.remove_from_parameter_map(text, geo)
            self.vtkrenderer.RemoveActor(geo['actor'])

        self.animate = True
        if children:
            self.geometrytree.setCurrentItem(children[0])
        else:
            i = self.geometrytree.topLevelItemCount() - 1
            item = self.geometrytree.topLevelItem(i)
            self.geometrytree.setCurrentItem(item)

        self.render()

    def handle_copy_geometry(self):
        """duplicate the selected geometry"""
        currentSelection = self.geometrytree.selectedItems()
        if currentSelection:
            text = str(currentSelection[-1].text(0)).lower()
            self.copy_geometry(text)

    def copy_geometry(self, name, center=None, rotation=None):
        """given a geometry name, copy it and optionaly change the center or
        rotation"""
        data = copy.deepcopy(remove_vtk_objects(self.geometrydict[name]))

        if center is not None:
            for key, val in zip(['centerx', 'centery', 'centerz'], center):
                data[key] = val

        if rotation is not None:
            for key, val in zip(['rotationx', 'rotationy', 'rotationz'], rotation):
                data[key] = val

        name = None
        geo_type = data['geo_type']
        if geo_type == 'primitive':
            name = self.add_primitive(data=data)
        elif geo_type == 'parametric':
            name = self.add_parametric(data=data)
        elif geo_type == 'filter':
            name = self.add_filter(data=data)
        elif geo_type == 'boolean':
            name = self.boolean_operation(data=data)
        elif geo_type == 'stl':
            name = self.add_stl(None, filename=data['filename'], data=data)
        elif geo_type == 'implicit':
            name = self.add_implicit(data=data)
        return name

    def update_filter(self, name):
        """Update the currently selected filter"""
        geo = self.geometrydict.get(name)
        filtertype = geo['type']
        vtkfilter = geo['filter']

        if filtertype == 'clean':
            vtkfilter.SetConvertLinesToPoints(safe_int(geo['linestopoints']))
            vtkfilter.SetConvertPolysToLines(safe_int(geo['polystolines']))
            vtkfilter.SetConvertStripsToPolys(safe_int(geo['stripstopolys']))
        elif filtertype == 'fill_holes':
            vtkfilter.SetHoleSize(safe_float(geo['maximumholesize']))
        elif filtertype == 'triangle':
            vtkfilter.SetPassVerts(safe_int(geo['processvertices']))
            vtkfilter.SetPassLines(safe_int(geo['processlines']))
        elif filtertype == 'decimate':
            vtkfilter.SetTargetReduction(safe_float(geo['targetreduction']))
        elif filtertype == 'quadric_decimation':
            vtkfilter.SetTargetReduction(safe_float(geo['targetreduction']))
        elif filtertype == 'quadric_clustering':
            vtkfilter.SetNumberOfXDivisions(safe_int(geo['divisionsx']))
            vtkfilter.SetNumberOfYDivisions(safe_int(geo['divisionsy']))
            vtkfilter.SetNumberOfZDivisions(safe_int(geo['divisionsz']))
            vtkfilter.SetAutoAdjustNumberOfDivisions(safe_int(geo['autoadjustdivisions']))
        elif filtertype == 'smooth':
            vtkfilter.SetEdgeAngle(safe_float(geo['edgeangle']))
            vtkfilter.SetFeatureAngle(safe_float(geo['featureangle']))
            vtkfilter.SetNumberOfIterations(safe_int(geo['iterations']))
            vtkfilter.SetRelaxationFactor(safe_float(geo['relaxation']))
            vtkfilter.SetFeatureEdgeSmoothing(safe_int(geo['featureedgesmoothing']))
            vtkfilter.SetBoundarySmoothing(safe_int(geo['boundarysmoothing']))
        elif filtertype == 'windowed_sinc':
            vtkfilter.SetEdgeAngle(safe_float(geo['edgeangle']))
            vtkfilter.SetFeatureAngle(safe_float(geo['featureangle']))
            vtkfilter.SetNumberOfIterations(safe_int(geo['iterations']))
            vtkfilter.SetPassBand(safe_float(geo['passband']))
            vtkfilter.SetFeatureEdgeSmoothing(safe_int(geo['featureedgesmoothing']))
            vtkfilter.SetBoundarySmoothing(safe_int(geo['boundarysmoothing']))
            vtkfilter.SetNonManifoldSmoothing(safe_int(geo['manifoldsmoothing']))
            vtkfilter.SetNormalizeCoordinates(safe_int(geo['normalize']))
        elif filtertype == 'reverse_sense':
            vtkfilter.SetReverseCells(safe_int(geo['reversecells']))
            vtkfilter.SetReverseNormals(safe_int(geo['reversenormals']))
        elif filtertype == 'transform':
            transform = geo['transform']
            # reset to Identity
            transform.Identity()
            transform.PostMultiply()

            # find the center
            polydata = vtkfilter.GetInput()
            com = vtk.vtkCenterOfMass()
            com.SetInputData(polydata)
            com.SetUseScalarsAsWeights(False)
            com.Update()
            x, y, z = com.GetCenter()

            # translate to center
            transform.Translate(-x, -y, -z)

            # scale
            transform.Scale(safe_float(geo['scalex']),
                            safe_float(geo['scaley']),
                            safe_float(geo['scalez']))

            # rotation
            transform.RotateWXYZ(safe_float(geo['rotationx']), 1, 0, 0)
            transform.RotateWXYZ(safe_float(geo['rotationy']), 0, 1, 0)
            transform.RotateWXYZ(safe_float(geo['rotationz']), 0, 0, 1)

            # translate
            transform.Translate(x + safe_float(geo['translatex']),
                                y + safe_float(geo['translatey']),
                                z + safe_float(geo['translatez']))
        elif filtertype == 'sample_implicit':
            geo['samplefunction'].SetSampleDimensions(
                safe_int(geo['samplesx']), safe_int(geo['samplesy']), safe_int(geo['samplesz']))
            bounds = [safe_float(geo[k]) for k in ['minx','maxx','miny', 'maxy', 'minz', 'maxz']]
            geo['samplefunction'].SetModelBounds(bounds)
            geo['samplefunction'].Update()


        vtkfilter.Update()

    def enable_filters(self):
        """filter menu about to show, enable/disable based on geo selected"""

        current_selection = self.geometrytree.selectedItems()
        if current_selection:
            name = str(current_selection[-1].text(0)).lower()
        else:
            return

        geo_data = self.geometrydict.get(name)
        geo_type = geo_data.get('geo_type')

        enable = [True]*len(self.filter_actions)
        if 'implicit' in geo_type:
            enable = [False]*len(self.filter_actions)
            enable[0] = True

        for en, act in zip(enable, self.filter_actions):
            act.setEnabled(en)

    def add_filter(self, filtertype=None, name=None, data=None, child=None, loading=False):
        """add the selected filter with the input being the currently selected
        toplevel item"""

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

        geo = self.geometrydict.get(name)

        # init filter
        if data is not None:
            filtertype = data['type']
        vtkfilter = geo['filter'] = FILTER_DICT[filtertype]()

        if filtertype == 'transform':
            t = geo['transform'] = vtk.vtkTransform()
            vtkfilter.SetTransform(t)

        # set input data
        if 'implicit' in filtertype:
            sample = vtk.vtkSampleFunction()
            sample.SetImplicitFunction(self.geometrydict.get(selection_text).get('source'))
            sample.ComputeNormalsOff()
            geo['samplefunction'] = sample
            vtkfilter.SetInputConnection(sample.GetOutputPort())
            vtkfilter.SetValue(0, 0.0)
        else:
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
        geo['actor'] = actor
        geo['mapper'] = mapper


        self.animate = False
        # Add to tree
        toplevel = QtWidgets.QTreeWidgetItem([name])
        toplevel.setIcon(0, get_icon('filter.png'))
        toplevel.setFlags(toplevel.flags() | QtCore.Qt.ItemIsUserCheckable)
        toplevel.setCheckState(0, QtCore.Qt.Checked)

        # remove children from tree
        for select in current_selection:
            toplevelindex = self.geometrytree.indexOfTopLevelItem(select)
            item = self.geometrytree.takeTopLevelItem(toplevelindex)
            item.setCheckState(0, QtCore.Qt.Unchecked)
            toplevel.addChild(item)

        self.geometrytree.addTopLevelItem(toplevel)
        self.animate = True
        self.geometrytree.setCurrentItem(toplevel)

        if not loading:
            self.parent.set_unsaved_flag()

        return name

    def get_input_data(self, name):
        """based on the type of geometry, return the data"""
        geo = self.geometrydict.get(name)
        if 'trianglefilter' in geo:
            inputdata = geo['trianglefilter']
        elif 'booleanoperation' in geo:
            inputdata = geo['booleanoperation']
        elif 'reader' in geo:
            inputdata = geo['transformfilter']
        elif 'filter' in geo:
            inputdata = geo['filter']
        return inputdata

    def collect_toplevel_geometry(self):
        """collect and append visible toplevel polydata"""

        append_filter = vtk.vtkAppendPolyData()
        item_count = self.geometrytree.topLevelItemCount()

        for top_num in range(item_count):
            item = self.geometrytree.topLevelItem(top_num)
            name = str(item.text(0))
            geo_type = self.geometrydict.get(name).get('geo_type', '')
            if item.checkState(0) == QtCore.Qt.Checked and not 'implicit' in geo_type:
                append_filter.AddInputData(
                    self.get_input_data(name).GetOutput())

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
        """determine the extents of the visible geometry"""

        geometry = self.collect_toplevel_geometry()

        bounds = None
        if geometry:
            bounds = geometry.GetOutput().GetBounds()

        return bounds

    def set_geometry_actor_props(self, actor, name):
        """set the geometry proprerties to the others in the scene"""

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

    def remove_from_parameter_map(self, name, del_geo):
        for key, value in del_geo.items():
            name_key = ','.join([name, key])
            self._remove_key(name_key, value)

    def _remove_key(self, name_key, value):
        if not isinstance(value, Equation): return

        for param in value.get_used_parameters():
            self.parameter_key_map[param].remove(name_key)
            if len(self.parameter_key_map[param]) == 0:
                self.parameter_key_map.pop(param)

    # --- regions ---
    def update_region_source(self, name):
        """Update the specified primitive"""
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
        """create a new region"""
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
        """delete a region"""
        region = self.region_dict.pop(name)
        self.vtkrenderer.RemoveActor(region['actor'])
        if 'clip_actor' in region:
            self.vtkrenderer.RemoveActor(region['clip_actor'])

    def update_region(self, name, region):
        """update a region"""
        self.region_dict[name].update(copy.deepcopy(region))
        self.update_region_source(name)
        self.select_facets(name)
        self.render()

    def change_region_color(self, name, color):
        """change the color of a region"""
        reg = self.region_dict.get(name)
        reg['color'] = copy.deepcopy(color)

        actor = reg['actor']
        actor.GetProperty().SetColor(*reg['color'].color_float)

        if 'clip_actor' in self.region_dict[name]:
            actor = reg['clip_actor']
            actor.GetProperty().SetColor(*reg['color'].color_float)

        self.render()

    def change_region_type(self, name, region):
        """change the type of a region"""

        self.region_dict[name].update(copy.deepcopy(region))

        if region['type'] == 'point':
            shape = 'sphere'
        else:
            shape = 'box'

        source = PRIMITIVE_DICT[shape]()
        region_data = self.region_dict.get(name)
        region_data['source'] = source
        self.update_region_source(name)
        self.select_facets(name)
        region_data['mapper'].SetInputConnection(source.GetOutputPort())

        self.render()

    def change_region_name(self, old_name, new_name):
        """change the name of a region"""
        region = self.region_dict.pop(old_name)
        self.region_dict[new_name] = region

    def change_region_visibility(self, name, visible):
        """change the visibility of a region"""
        reg = self.region_dict.get(name)
        if visible and self.visual_props['regions']['visible']:
            reg['actor'].VisibilityOn()
            if 'clip_actor' in reg:
                reg['clip_actor'].VisibilityOn()
        else:
            reg['actor'].VisibilityOff()
            if 'clip_actor' in reg:
                reg['clip_actor'].VisibilityOff()
        reg['visible'] = visible

        self.render()

    def set_region_actor_props(self, actor, name, color=None):
        """set the geometry properties to the others in the scene"""

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
        """select facets with an implicit function"""

        region = self.region_dict[name]
        # remove old objects
        if 'clip_actor' in region:
            self.vtkrenderer.RemoveActor(region['clip_actor'])
            for key in ['clip_actor', 'clip_mapper', 'clipper', 'implicit']:
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
        """export visible toplevel geometry"""

        basename = os.path.splitext(os.path.basename(file_name))[0]
        basename = os.path.join(os.path.dirname(file_name), basename)

        # Export geometry
        geometry = self.collect_toplevel_geometry()
        if geometry:
            GUI.update_keyword('cartesian_grid', True)
            GUI.update_keyword('use_stl', True)
            # write file
            stl_writer = vtk.vtkSTLWriter()
            stl_writer.SetFileName(file_name)
            stl_writer.SetInputConnection(geometry.GetOutputPort())
            stl_writer.Write()
        else:
            GUI.update_keyword('cartesian_grid', False)
            GUI.update_keyword('use_stl', False)

        # write name_000n.stl for stl regions
        if GUI.bcs:
            for i, bc_data in GUI.bcs.items():
                region = bc_data['region']
                data = {}

                if region in self.region_dict:
                    data = self.region_dict[region]
                else:
                    continue
                if data['type'] == 'STL' and 'clipper' in data:
                    fname = basename + '_{0:04d}.stl'.format(i)
                    stl_writer = vtk.vtkSTLWriter()
                    stl_writer.SetFileName(fname)
                    if data['slice']:
                        stl_writer.SetInputConnection(data['clipper'].GetClippedOutputPort())
                    else:
                        stl_writer.SetInputConnection(data['clipper'].GetOutputPort())
                    stl_writer.Write()

    def export_unstructured(self, fname, grid):
        """export an unstructured grid"""
        gw = vtk.vtkXMLUnstructuredGridWriter()
        gw.SetFileName(fname)
        gw.SetInputConnection(grid)
        gw.Write()

    def screenshot(self, fname=None):
        """take a snapshot of the vtk window"""
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
    def auto_size_mesh_extents(self):
        """collect and set the extents of the visible geometry"""
        extents = self.get_geometry_extents()

        if extents:
            for key, extent in zip(['xmin', 'xlength', 'ymin', 'ylength',
                                    'zmin', 'zlength'],
                                   extents):
                if 'min' not in key:  # mfix doesn't support mins yet
                    self.emitUpdatedValue(key, extent)

#            self.update_background_mesh()

    def update_background_mesh(self, spacing):
        """update the background mesh"""
        cells = [len(c) for c in spacing]
        self.rectilinear_grid.SetDimensions(*cells)

        x_coords = vtk.vtkFloatArray()
        for i in spacing[0]:
            x_coords.InsertNextValue(i)

        y_coords = vtk.vtkFloatArray()
        for i in spacing[1]:
            y_coords.InsertNextValue(i)

        z_coords = vtk.vtkFloatArray()
        for i in spacing[2]:
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
        """set the background mesh properties"""
        props = self.visual_props['background_mesh']
        self.set_representation(actor, props['rep'])
        actor.GetProperty().SetColor(props['color'].getRgbF()[:3])
        actor.GetProperty().SetEdgeColor(props['color'].getRgbF()[:3])
        actor.GetProperty().SetOpacity(props['opacity'])
        actor.SetVisibility(int(props['visible']))

    def vtk_calc_distance_from_geometry(self):
        """for every point in the mesh, calculate the distance from the
        toplevel geometry"""
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
        """use vtkTableBasedClipDataSet to slice the background mesh with the
        toplevel geometry"""
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
        """use vtkThreshold to slice the background mesh with the
        toplevel geometry"""
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
        """export the mesh"""
        project_dir = self.parent.get_project_dir()
        self.export_unstructured(os.path.join(project_dir, name),
                                 mesh.GetOutputPort())

    def mesh_stats(self):
        """calculate mesh statistics"""
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
        """mesh the geometry"""
        mesher = str(self.ui.mesh.combobox_mesher.currentText())

        if mesher == 'vtkTableBasedClipDataSet':
            self.vtk_mesher_table_based()
        elif mesher == 'vtkThreshold':
            self.vtk_mesher_threshold()

    def remove_mesh(self, name=DEFAULT_MESH_NAME):
        """remove the mesh from the vtk scene as well as the file"""
        project_dir = self.parent.get_project_dir()
        path = os.path.join(project_dir, name)
        if not os.path.exists(path):
            return
        key = GUI.message(text='Remove %s?' % name, buttons=['yes', 'no'], default='no')
        if key == 'no':
            return

        os.remove(path)

        if self.mesh_actor is not None:
            self.vtkrenderer.RemoveActor(self.mesh_actor)
            self.mesh_actor = None

    def set_mesh_actor_props(self):
        """change the actor properties of the mesh"""
        props = self.visual_props['mesh']
        self.set_representation(self.mesh_actor, props['rep'])
        self.mesh_actor.GetProperty().SetColor(props['color'].getRgbF()[:3])
        self.mesh_actor.GetProperty().SetEdgeColor(props['edge'].getRgbF()[:3])
        self.mesh_actor.SetVisibility(int(props['visible']))

    def show_mesh(self):
        """show the mesh"""
        if self.mesh is None:
            return

        if self.mesh_actor is not None:
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
    def change_interaction(self, style_2d=False):
        if style_2d:
            self.vtkiren.SetInteractorStyle(self.style_2d)
            self.view_flip[1] = False
            self.set_view()
            enabled = False
        else:
            self.vtkiren.SetInteractorStyle(self.style_3d)
            self.perspective(False)
            enabled = True

        for btn in [self.toolbutton_perspective, self.toolbutton_view_yz,
                    self.toolbutton_view_xz]:
            btn.setEnabled(enabled)

    def perspective(self, parallel=None):
        """change the perspective of the vtk scene"""
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
        """set the 2D view of the scene"""
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
        """reset the camera so that the geometry is visible in the scene"""
        self.vtkrenderer.ResetCamera()
        self.render()

    def change_visibility(self, name, toolbutton):
        """change the visibility of one of the scene objects"""
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
                if actor is not None:
                    if visible:
                        actor.VisibilityOn()
                    else:
                        actor.VisibilityOff()
            self.render()

    def get_actors(self, name):
        """given a scene actor type, return the actors"""
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
        """set the representation of an actor"""
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
        """given a scene actor type, change the representation"""
        representation = str(combobox.currentText())

        self.visual_props[name]['rep'] = representation

        for actor in self.get_actors(name):
            self.set_representation(actor, representation)

        self.render()

    def change_color(self, name, button):
        """given a scene actor type, change the color"""
        col = QtWidgets.QColorDialog.getColor()
        if not col.isValid():
            return

        button.setStyleSheet("QToolButton{{ background: {};}}".format(
            col.name()))

        self.visual_props[name]['color'] = col
        dark = self.visual_props[name]['edge'] = col.darker()

        for actor in self.get_actors(name):
            if actor is not None:
                actor.GetProperty().SetColor(col.getRgbF()[:3])
                actor.GetProperty().SetEdgeColor(dark.getRgbF()[:3])

        self.render()

    def change_opacity(self, name, opacity):
        """given a scene actor type, change the opacity"""
        value = opacity.value()
        self.visual_props[name]['opacity'] = value

        for actor in self.get_actors(name):
            if actor is not None:
                actor.GetProperty().SetOpacity(value)
        self.render()

    def set_visible_btn_image(self, btn, checked):
        """given a button, change the icon"""
        if not checked:
            btn.setIcon(get_icon('visibilityofftransparent.png'))
        else:
            btn.setIcon(get_icon('visibility.png'))

    # --- wizards ---
    def handle_cyclone_wizard(self):
        """show the cyclone wizard"""
        if self.cyclone_popup is None:
            self.cyclone_popup = CyclonePopUp(self)
        self.cyclone_popup.popup()

    def handle_distributed_wizard(self):
        """show the distributed wizard"""
        if self.distribution_popup is None:
            self.distribution_popup = DistributionPopUp(self)
        self.distribution_popup.popup()

    def handle_reactor_wizard(self):
        """show the reactor wizard"""
        if self.reactor_popup is None:
            self.reactor_popup = ReactorPopUp(self)
        self.reactor_popup.popup()

    def handle_hopper_wizard(self):
        """show the hopper wizard"""
        if self.hopper_popup is None:
            self.hopper_popup = HopperPopUp(self)
        self.hopper_popup.popup()
