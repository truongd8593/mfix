from __future__ import print_function, absolute_import, unicode_literals, division
import logging
import os
import vtk
import numpy as np

from qtpy import QtWidgets, uic

log = logging.getLogger('mfix-gui' if __name__=='__main__' else __name__)

class ReactorPopUp(QtWidgets.QDialog):
    def __init__(self, parent):
        QtWidgets.QDialog.__init__(self, parent)

        self.vtk_widget = parent

        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        ui = self.ui = uic.loadUi(os.path.join(uidir, 'reactor.ui'), self)

        self.setWindowTitle('Reactor Wizard')

        ui.pushbutton_close.clicked.connect(self.close)

    def popup(self):

        self.show()
        self.raise_()
        self.activateWindow()