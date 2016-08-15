from __future__ import print_function, absolute_import, unicode_literals, division
import logging
import os
import vtk
import numpy as np

from qtpy import QtWidgets, uic
from tools.general import get_pixmap

log = logging.getLogger('mfix-gui' if __name__=='__main__' else __name__)

class CyclonePopUp(QtWidgets.QDialog):
    def __init__(self, parent):
        QtWidgets.QDialog.__init__(self, parent)

        self.vtk_widget = parent

        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        ui = self.ui = uic.loadUi(os.path.join(uidir, 'cyclone.ui'), self)

        self.setWindowTitle('Cyclone Wizard')

        ui.pushbutton_close.clicked.connect(self.close)
        pixmap = get_pixmap('cyclone_sketch.png')
        ui.label_image.setPixmap(pixmap)
        ui.label_image.setMask(pixmap.mask())
        ui.lineedit_dc.setFocus()

    def popup(self):

        self.show()
        self.raise_()
        self.activateWindow()
