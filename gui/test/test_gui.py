'''
  run with:

  > python -m unittest discover

'''

import os
import time
import unittest

from qtpy.QtTest import QTest
from qtpy import QtCore
from qtpy import QtGui
from qtpy import QtWidgets

import gui

class MfixGuiTests(unittest.TestCase):
    ''' unit tests for the GUI '''

    def dismiss(self):
        ''' dismiss modal QMessageBox '''
        for widget in QtGui.qApp.topLevelWidgets():
            if isinstance(widget, QtWidgets.QMessageBox):
                QTest.keyClick(widget, QtCore.Qt.Key_Enter)

    def setUp(self):
        ''' open FluidBed_DES for testing '''
        qapp = QtWidgets.QApplication([])
        self.mfix = gui.MfixGui(qapp)
        self.mfix.show()
        QTest.qWaitForWindowShown(self.mfix)

        arg = '../tutorials/FluidBed_DES/mfix.dat'

        QtCore.QTimer.singleShot(1000, self.dismiss)
        self.mfix.open_project(arg, True)
        self.assertEqual("DES_FB1", self.mfix.ui.general.lineedit_keyword_run_name.text())
        mfxfile = '../tutorials/FluidBed_DES/DES_FB1.mfx'
        self.assertTrue(os.path.exists(mfxfile))

    def test_save_as(self):
        ''' Test the Save As button on the toolbar '''

        newname = 'DES_FB1_new_name'

        self.mfix.get_save_filename = lambda : newname
        QtCore.QTimer.singleShot(1000, self.dismiss)
        QTest.mouseClick(self.mfix.ui.toolbutton_saveas, QtCore.Qt.LeftButton)

        self.assertEqual(newname, self.mfix.ui.general.lineedit_keyword_run_name.text())

    def test_run(self):
        ''' Test the Run button on the toolbar '''

        QtCore.QTimer.singleShot(1000, self.dismiss)
        QTest.mouseClick(self.mfix.ui.toolbutton_run, QtCore.Qt.LeftButton)
        time.sleep(1)

        logfile = '../tutorials/FluidBed_DES/DES_FB1.LOG'
        self.assertTrue(os.path.exists(logfile))
