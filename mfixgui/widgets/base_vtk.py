"""
This is a base class for the vtk widgets.
"""
from __future__ import print_function, absolute_import, unicode_literals, division
import traceback
import logging
LOG = logging.getLogger(__name__)

from qtpy import QtCore, QtGui, QtWidgets

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


class BaseVtkWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

        # --- layout ---
        self.grid_layout = QtWidgets.QGridLayout(self)
        self.grid_layout.setContentsMargins(0, 0, 0, 0)

        self.vtkWindowWidget = QVTKRenderWindowInteractor(self)
        self.grid_layout.addWidget(self.vtkWindowWidget, 1, 0)

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

