# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from qtpy import QtCore, QtWidgets, QtGui, uic
import os
import fnmatch

from mfixgui.constants import RESTART_FILES, SPX_FILES, VTK_FILES, OTHER_FILES

class DeletePopup(QtWidgets.QDialog):
    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self, parent)

        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(thisdir, os.pardir, 'uifiles')
        self.ui = uic.loadUi(os.path.join(uidir, 'delete_popup.ui'), self)

        self.setWindowTitle('Delete Files')

    def exec_(self, files, force_remove=False):

        ui = self.ui
        enable = not force_remove
        ui.groupbox_res.setEnabled(enable)
        ui.groupbox_spx.setEnabled(enable)
        ui.groupbox_vtk.setEnabled(enable)
        ui.groupbox_other.setEnabled(enable)

        self.sort_files(files)
        return QtWidgets.QDialog.exec_(self)

    def sort_files(self, files):
        ui = self.ui
        ui.listwidget_restart.clear()
        ui.listwidget_spx.clear()
        ui.listwidget_vtk.clear()
        ui.listwidget_other.clear()

        ui.groupbox_res.setChecked(True)
        ui.groupbox_spx.setChecked(True)
        ui.groupbox_vtk.setChecked(True)
        ui.groupbox_other.setChecked(True)

        files = sorted(files)

        self.restart = []
        self.spx = []
        self.vtk = []
        self.other = []

        for s in RESTART_FILES:
            ff = fnmatch.filter(files, s)
            self.restart.extend(ff)
            ui.listwidget_restart.addItems(ff)
        for s in SPX_FILES:
            ff = fnmatch.filter(files, s)
            self.spx.extend(ff)
            ui.listwidget_spx.addItems(ff)
        for s in VTK_FILES:
            ff = fnmatch.filter(files, s)
            self.vtk.extend(ff)
            ui.listwidget_vtk.addItems(ff)
        for s in OTHER_FILES:
            ff = fnmatch.filter(files, s)
            self.other.extend(ff)
            ui.listwidget_other.addItems(ff)

    def get_output_files(self):
        ui = self.ui

        files = []
        for gb, fs in [(ui.groupbox_res, self.restart),
                       (ui.groupbox_spx, self.spx),
                       (ui.groupbox_vtk, self.vtk),
                       (ui.groupbox_other, self.other)]:
            if gb.isChecked():
                files.extend(fs)
        return files
