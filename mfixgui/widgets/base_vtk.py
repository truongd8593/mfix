"""
This is a base class for the vtk widgets.
"""
from __future__ import print_function, absolute_import, unicode_literals, division
import os
import traceback
import logging
LOG = logging.getLogger(__name__)

from qtpy import QtCore, QtWidgets, PYQT5

# VTK imports
VTK_AVAILABLE = True
try:
    import vtk
    LOG.info('VTK version: %s' % vtk.vtkVersion.GetVTKVersion())
    VTK_MAJOR_VERSION = vtk.VTK_MAJOR_VERSION
except ImportError:
    VTK_AVAILABLE = False
    vtk = None
    VTK_MAJOR_VERSION = None
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

from mfixgui.tools.general import get_icon

SETTINGS = QtCore.QSettings('MFIX', 'MFIX')

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

        self.defer_render = False
        self.view_flip = [False]*3
        self.offscreen_vtkrenderer = None

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

    def init_base_toolbar(self):
        """add base toolbar"""
        self.button_bar = QtWidgets.QWidget(self)
        self.button_bar_layout = QtWidgets.QHBoxLayout(self.button_bar)
        self.button_bar_layout.setContentsMargins(0, 0, 0, 0)
        self.button_bar.setLayout(self.button_bar_layout)
        self.button_bar.setGeometry(QtCore.QRect(0, 0, 300, 300))
        self.grid_layout.addWidget(self.button_bar, 0, 0)

        self.toolbutton_reset = QtWidgets.QToolButton()
        self.toolbutton_reset.clicked.connect(self.reset_view)
        self.toolbutton_reset.setIcon(get_icon('overscan.png'))

        self.toolbutton_perspective = QtWidgets.QToolButton()
        self.toolbutton_perspective.clicked.connect(lambda ignore: self.perspective())
        self.toolbutton_perspective.setIcon(get_icon('perspective.png'))

        self.toolbutton_view_xy = QtWidgets.QToolButton()
        self.toolbutton_view_xy.clicked.connect(lambda: self.set_view('xy'))
        self.toolbutton_view_xy.setIcon(get_icon('xy.png'))

        self.toolbutton_view_yz = QtWidgets.QToolButton()
        self.toolbutton_view_yz.clicked.connect(lambda: self.set_view('yz'))
        self.toolbutton_view_yz.setIcon(get_icon('yz.png'))

        self.toolbutton_view_xz = QtWidgets.QToolButton()
        self.toolbutton_view_xz.clicked.connect(lambda: self.set_view('xz'))
        self.toolbutton_view_xz.setIcon(get_icon('xz.png'))

        self.toolbutton_screenshot = QtWidgets.QToolButton()
        self.toolbutton_screenshot.clicked.connect(self.screenshot)
        self.toolbutton_screenshot.setIcon(get_icon('camera.png'))

        for btn in [self.toolbutton_reset,
                    self.toolbutton_view_xy,
                    self.toolbutton_view_yz,
                    self.toolbutton_view_xz,
                    self.toolbutton_perspective,
                    self.toolbutton_screenshot,]:
            self.button_bar_layout.addWidget(btn)
            btn.setAutoRaise(True)
            btn.setFocusPolicy(QtCore.Qt.ClickFocus)

    # --- render ---
    def render(self, force_render=False, defer_render=None):
        """render the vtk scene, checks for defer_render"""
        if defer_render is not None:
            self.defer_render = defer_render
        if not self.defer_render or force_render:
            self.vtkRenderWindow.Render()

    def screenshot(self, checked, fname=None, size=[1920, 1080]):
        """take a snapshot of the vtk window"""
        self.toolbutton_screenshot.setDown(False)

        if fname is None:
            fname = QtWidgets.QFileDialog.getSaveFileName(
                self,
                "Save screenshot",
                os.path.dirname(SETTINGS.value('project_file')),
                ';;'.join(["PNG (*.png)",
                           "JPEG (*.jpg)",
                           "PostScript (*.ps)",
                           "All Files (*.*)",
                          ]),
                )
            if PYQT5:
                fname = fname[0]
        if not fname:
            return

        notvisible = not self.isVisible()
        if notvisible:
            self.init_offscreen_render(size)

        # screenshot code:
        # TODO: get resolution from user
        window_image = vtk.vtkWindowToImageFilter()
        if notvisible:
            window_image.SetInput(self.offscreen_vtkrenderwindow)
        else:
            window_image.SetInput(self.vtkRenderWindow)
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
            # force to png
            writer = vtk.vtkPNGWriter()
            fname += '.png'

        writer.SetFileName(fname)
        writer.SetInputConnection(window_image.GetOutputPort())
        writer.Write()

    def init_offscreen_render(self, size=[1920, 1080]):
        self.offscreen_vtkrenderer = vtk.vtkRenderer()
        self.offscreen_vtkrenderer.GradientBackgroundOn()
        self.offscreen_vtkrenderer.SetBackground(0.4, 0.4, 0.4)
        self.offscreen_vtkrenderer.SetBackground2(1.0, 1.0, 1.0)
        self.offscreen_vtkrenderwindow = vtk.vtkRenderWindow()
        self.offscreen_vtkrenderwindow.SetAlphaBitPlanes(1)
        self.offscreen_vtkrenderwindow.SetOffScreenRendering(1)
        self.offscreen_vtkrenderwindow.AddRenderer(self.offscreen_vtkrenderer)
        self.offscreen_vtkrenderwindow.SetSize(*size)

        actors = self.vtkrenderer.GetActors()

        for i in range(actors.GetNumberOfItems()):
            actor = actors.GetItemAsObject(i)
            if actor:
                mapper = actor.GetMapper()
                new_actor = vtk.vtkActor()
                new_actor.SetMapper(mapper)
                new_actor.SetProperty(actor.GetProperty())
                self.offscreen_vtkrenderer.AddActor(new_actor)

        self.offscreen_vtkrenderer.ResetCamera()

        self.offscreen_vtkrenderer.Render()

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
