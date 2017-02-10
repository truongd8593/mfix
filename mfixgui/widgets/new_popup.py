# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
import logging
import os
import sys

from qtpy import QtWidgets, QtCore, uic
from mfixgui.tools.general import get_icon

log = logging.getLogger('mfix-gui' if __name__=='__main__' else __name__)
SETTINGS = QtCore.QSettings('MFIX', 'MFIX')

class NewProjectDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        QtWidgets.QDialog.__init__(self, parent)

        self.vtk_widget = parent

        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        ui = self.ui = uic.loadUi(os.path.join(uidir, 'new_project_popup.ui'), self)

        self.setWindowTitle('Create a new project')

        ui.toolbutton_browse.clicked.connect(self.browse)
        ui.toolbutton_browse.setIcon(get_icon('folder.png'))

    def get(self, run_name='new_project'):

        ui = self.ui
        ui.lineedit_project_name.setText(run_name)

        items = SETTINGS.value('project_locations', '')
        if items:
            ui.combobox_location.addItems(items.split(','))

        ok = self.exec_()

        # save items
        items = self.get_items()
        # only save 5 most recent ones
        SETTINGS.setValue('project_locations', ','.join(items[-5:]))

        return ok == QtWidgets.QDialog.Accepted, ui.lineedit_project_name.text(), ui.combobox_location.currentText()

    def get_items(self):
        ui = self.ui
        return [ui.combobox_location.itemText(i) for i in range(ui.combobox_location.count())]

    def browse(self):
        cb = self.ui.combobox_location
        loc = QtWidgets.QFileDialog.getExistingDirectory(self, 'Location', cb.currentText())

        if loc:
            if loc not in self.get_items():
                cb.addItem(loc)
            cb.setCurrentIndex(cb.findText(loc))


if __name__  == '__main__':

     # create the QApplication
    qapp = QtWidgets.QApplication([])

    ok, run_name, project_loc = NewProjectDialog(None).get('project')
    print(ok, run_name, project_loc)
    sys.exit()