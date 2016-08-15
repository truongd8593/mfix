import logging
import os
import vtk
import numpy as np

from qtpy import QtWidgets, uic

log = logging.getLogger('mfix-gui' if __name__=='__main__' else __name__)

class DistributionPopUp(QtWidgets.QDialog):
    def __init__(self, parent):
        QtWidgets.QDialog.__init__(self, parent)
        
        self.vtk_widget = parent
        
        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        ui = self.ui = uic.loadUi(os.path.join(uidir, 'distributed.ui'), self)
        
        ui.pushbutton_cancel.clicked.connect(self.close)
        ui.pushbutton_apply.clicked.connect(self.apply_)
        ui.combobox_distribution.currentIndexChanged.connect(self.dist_changed)
        
    def dist_changed(self):
        enable = self.ui.combobox_distribution.currentText() == 'random'
        self.ui.lineedit_number.setEnabled(enable)
        if enable:
            self.ui.label_numberspacing.setText('number')
        else:
            self.ui.label_numberspacing.setText('spacing')
        
    def popup(self):
        geo = list(self.vtk_widget.geometrydict.keys()) 
        
        self.ui.combobox_shape.clear()
        self.ui.combobox_shape.addItems(geo)
        
        self.ui.combobox_volume.clear()
        self.ui.combobox_volume.addItems(geo)
        
        self.show()
        self.raise_()
        self.activateWindow()
        
    def apply_(self):
        
        geo = self.ui.combobox_volume.currentText()
        shape = self.ui.combobox_shape.currentText()
        
        if self.ui.combobox_distribution.currentText() =='random':
            self.gen_random(geo, shape)
    
#        self.close()
    
    def gen_random(self, geo, shape):
        
        n_points = 0
        max_itr = 0
        bounds = self.get_bounds(geo)
        last = None
        while n_points < int(self.ui.lineedit_number.text()) and max_itr <1000000000:
            point = np.random.random(3)
            point = [p*(t-f)+f for p, f, t in zip(point, bounds[::2], bounds[1::2])]
            for p in self.points_in_geo(geo, [point]):
                new = self.vtk_widget.copy_geometry(shape, p)
                if self.ui.checkbox_union.isChecked() and last is not None:
                    last = self.vtk_widget.boolean_operation(
                        booltype='union', children=[last, new])
                else:
                    last = new
                n_points+=1
            max_itr +=1
            
    def get_bounds(self, geo):
        return self.vtk_widget.get_input_data(geo).GetOutput().GetBounds()
        
    def points_in_geo(self, geo, points):
        
        vtkpoints = vtk.vtkPoints()

        for point in points:
            vtkpoints.InsertNextPoint(point)
            
        vtkpolydata = vtk.vtkPolyData()
        vtkpolydata.SetPoints(vtkpoints)
        
        enclosed_points = vtk.vtkSelectEnclosedPoints()
        enclosed_points.SetInputData(vtkpolydata)
        enclosed_points.SetSurfaceData(self.vtk_widget.get_input_data(geo).GetOutput())
        enclosed_points.Update()
        
        for i in range(len(points)):
            if enclosed_points.IsInside(i):
                yield points[i]