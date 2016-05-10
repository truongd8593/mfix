# -*- coding: UTF-8 -*-
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
        for widget in QtWidgets.QApplication.instance().topLevelWidgets():
            if isinstance(widget, QtWidgets.QMessageBox):
                QTest.keyClick(widget, QtCore.Qt.Key_Enter)

    def setUp(self):
        ''' open FluidBed_DES for testing '''
        qapp = QtWidgets.QApplication([])
        self.mfix = gui.MfixGui(qapp)
        self.mfix.show()
        QTest.qWaitForWindowShown(self.mfix)

        project_path = '../tutorials/FluidBed_DES/mfix.dat'

        self.mfix.get_open_filename = lambda : project_path
        QtCore.QTimer.singleShot(100, self.dismiss)
        self.mfix.handle_open_action()

        self.assertEqual("DES_FB1", self.mfix.ui.general.lineedit_keyword_run_name.text())
        mfxfile = '../tutorials/FluidBed_DES/DES_FB1.mfx'
        self.assertTrue(os.path.exists(mfxfile))

    def test_save_as(self):
        ''' Test the Save As button on the toolbar '''

        newname = 'DES_FB1_new_name'

        self.mfix.get_save_filename = lambda : newname
        QtCore.QTimer.singleShot(100, self.dismiss)
        QTest.mouseClick(self.mfix.ui.toolbutton_save_as, QtCore.Qt.LeftButton)

        self.assertEqual(newname, self.mfix.ui.general.lineedit_keyword_run_name.text())

    def test_run(self):
        ''' Test the Run button on the toolbar '''

        QtCore.QTimer.singleShot(100, self.dismiss)
        QTest.mouseClick(self.mfix.ui.toolbutton_run, QtCore.Qt.LeftButton)
        time.sleep(1)

        logfile = '../tutorials/FluidBed_DES/DES_FB1.LOG'
        self.assertTrue(os.path.exists(logfile))

    def test_description_unicode(self):
        ''' Test setting TSTOP on the Run pane '''

        run_treeitem = None
        for item in self.mfix.ui.treewidget_model_navigation.findItems(
                    'general', QtCore.Qt.MatchContains | QtCore.Qt.MatchRecursive, 0):
            run_treeitem = item

        self.mfix.ui.treewidget_model_navigation.setCurrentItem(run_treeitem)

        description = self.mfix.ui.general.combobox_keyword_description.value
        new_description = description + u'Παν語'
        self.mfix.ui.general.combobox_keyword_description.setFocus()
        QTest.keyClicks(self.mfix.ui.general.combobox_keyword_description, 'new description')
        QTest.mouseClick(self.mfix.ui.toolbutton_save, QtCore.Qt.LeftButton)

        found = 0
        with open('../tutorials/FluidBed_DES/DES_FB1.mfx') as ff:
            for line in ff.readlines():
                kv = line.split('=')
                if len(kv) > 1 and kv[0].strip()=='description':
                    self.assertEqual(kv[1].strip()[1:-1], new_description)
                    found += 1

        self.assertEquals(found, 1)

    def test_tstop(self):
        ''' Test setting TSTOP on the Run pane '''

        run_treeitem = None
        for item in self.mfix.ui.treewidget_model_navigation.findItems(
                    'run', QtCore.Qt.MatchContains | QtCore.Qt.MatchRecursive, 0):
            run_treeitem = item

        self.mfix.ui.treewidget_model_navigation.setCurrentItem(run_treeitem)

        dt = self.mfix.ui.run.doublespinbox_keyword_dt.value()
        new_tstop = 5*dt
        self.mfix.ui.run.doublespinbox_keyword_tstop.setValue(new_tstop)
        QTest.mouseClick(self.mfix.ui.toolbutton_save, QtCore.Qt.LeftButton)

        found = 0
        with open('../tutorials/FluidBed_DES/DES_FB1.mfx') as ff:
            for line in ff.readlines():
                kv = line.split('=')
                if len(kv) > 1 and kv[0].strip()=='tstop':
                    self.assertAlmostEqual(float(kv[1]), new_tstop)
                    found += 1

        self.assertEquals(found, 1)
