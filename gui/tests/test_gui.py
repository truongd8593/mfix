# -*- coding: utf-8 -*-
'''
  run with:

  > python -m unittest discover

  or nosetests

'''

import fnmatch
import glob
import os
import time
import unittest
#from xvfbwrapper import Xvfb

from qtpy.QtTest import QTest
from qtpy import QtCore
from qtpy import QtWidgets

from .helper_functions import TestQApplication
import gui

def dismiss():
    ''' dismiss modal QMessageBox '''
    for widget in QtWidgets.QApplication.instance().topLevelWidgets():
        if isinstance(widget, QtWidgets.QMessageBox):
            QTest.keyClick(widget, QtCore.Qt.Key_Enter)


class MfixGuiTests(TestQApplication):
    ''' unit tests for the GUI '''

    def find_exes(self):
        """find all mfix and pymfix executables"""
        matches = []
        for root, dirnames, filenames in os.walk(self.mfix_home):
            for filename in fnmatch.filter(filenames, 'mfix'):
                matches.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, 'pymfix'):
                matches.append(os.path.join(root, filename))
        return matches

    def setUp(self):
        ''' open FluidBed_DES for testing '''

        #self.xvfb = Xvfb(width=1280, height=720)
        #self.addCleanup(self.xvfb.stop)
        #self.xvfb.start()

        self.rundir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.mfix_home = os.path.dirname(self.rundir)
        self.rundir = os.path.join(self.mfix_home, 'tutorials', 'FluidBed_DES')
        mfix_dat = os.path.join(self.rundir, 'mfix.dat')

        patterns = [
            '*.LOG', '*.OUT', '*.RES', '*.SP?', 'DES_FB1*', '*.mfx',
            '*.pvd', '*.vtp', 'VTU_FRAME_INDEX.TXT']
        for paths in [glob.glob(os.path.join(self.rundir, n)) for n in patterns]:
            for path in paths:
                try:
                    os.remove(path)
                except OSError:
                    pass

        TestQApplication.setUp(self)
        self.mfix = gui.MfixGui(self.qapp)
        self.mfix.show()
        QTest.qWaitForWindowShown(self.mfix)

        self.mfix.get_open_filename = lambda : mfix_dat
        QtCore.QTimer.singleShot(100, dismiss)
        QTest.mouseClick(self.mfix.ui.toolbutton_open, QtCore.Qt.LeftButton)

        self.assertEqual("DES_FB1", self.mfix.ui.general.lineedit_keyword_run_name.text())
        mfxfile = os.path.join(self.rundir, 'DES_FB1.mfx')
        self.assertTrue(os.path.exists(mfxfile))

    def tearDown(self):
        patterns = [
            '*.LOG', '*.OUT', '*.RES', '*.SP?', 'DES_FB1*', '*.mfx',
            '*.pvd', '*.vtp', 'VTU_FRAME_INDEX.TXT']
        for paths in [glob.glob(os.path.join(self.rundir, n)) for n in patterns]:
            for path in paths:
                try:
                    os.remove(path)
                except OSError:
                    pass

        TestQApplication.tearDown(self)

    def test_save_as(self):
        self.skipTest("FIXME")
        newname = 'DES_FB1_new_name'
        newpath = os.path.join(self.rundir, newname)

        self.mfix.get_save_filename = lambda : newpath
        QtCore.QTimer.singleShot(100, dismiss)
        QTest.mouseClick(self.mfix.ui.toolbutton_save_as, QtCore.Qt.LeftButton)

        self.assertEqual(newname, self.mfix.ui.general.lineedit_keyword_run_name.text())
        mfxfile = os.path.join(self.rundir, '%s.mfx' % newname)
        self.assertTrue(os.path.exists(mfxfile))

    def test_run(self):
        self.skipTest("FIXME, where did the run button go?")
        if not self.find_exes(): # We should be able to mock this for testing
            self.skipTest("Only valid when executables are present")

        # before running
        self.assertTrue(self.mfix.ui.run.spinbox_mfix_executables.isVisibleTo(self.mfix.ui.run))
        self.assertFalse(self.mfix.ui.run.resume_mfix_button.isEnabled())
        self.assertTrue(self.mfix.ui.run.run_mfix_button.isEnabled())
        self.assertFalse(self.mfix.ui.run.stop_mfix_button.isEnabled())

        QTest.mouseClick(self.mfix.ui.toolbutton_run, QtCore.Qt.LeftButton)
        time.sleep(1)

        # during running
        self.assertFalse(self.mfix.ui.run.resume_mfix_button.isEnabled())
        self.assertFalse(self.mfix.ui.run.run_mfix_button.isEnabled())
        self.assertTrue(self.mfix.ui.run.stop_mfix_button.isEnabled())

        QTest.mouseClick(self.mfix.ui.run.stop_mfix_button, QtCore.Qt.LeftButton)

        # after running
        self.assertTrue(self.mfix.ui.run.resume_mfix_button.isEnabled())
        self.assertTrue(self.mfix.ui.run.run_mfix_button.isEnabled())
        self.assertFalse(self.mfix.ui.run.stop_mfix_button.isEnabled())

        QTest.mouseClick(self.mfix.ui.run.resume_mfix_button, QtCore.Qt.LeftButton)

        # after resuming
        self.assertFalse(self.mfix.ui.run.resume_mfix_button.isEnabled())
        self.assertFalse(self.mfix.ui.run.run_mfix_button.isEnabled())
        self.assertTrue(self.mfix.ui.run.stop_mfix_button.isEnabled())

        QTest.mouseClick(self.mfix.ui.run.stop_mfix_button, QtCore.Qt.LeftButton)

        logfile = os.path.join(self.rundir, 'DES_FB1.LOG')
        self.assertTrue(os.path.exists(logfile))

    def test_description_unicode(self):
        run_treeitem = None
        for item in self.mfix.ui.treewidget_model_navigation.findItems(
                    'general', QtCore.Qt.MatchContains | QtCore.Qt.MatchRecursive, 0):
            run_treeitem = item

        self.mfix.ui.treewidget_model_navigation.setCurrentItem(run_treeitem)

        description = self.mfix.ui.general.combobox_keyword_description.value
        self.mfix.ui.general.combobox_keyword_description.setFocus()
        QTest.keyClick(self.mfix.ui.general.combobox_keyword_description, QtCore.Qt.Key_Right)
        # FIXME get Qt to accept non-ASCII text
        # new_text = u'Παν語'
        new_text = u'some new text'
        QTest.keyClicks(self.mfix.ui.general.combobox_keyword_description, new_text)
        QTest.mouseClick(self.mfix.ui.toolbutton_save, QtCore.Qt.LeftButton)

        new_description = str(description + new_text)

        found = 0
        with open(os.path.join(self.rundir,'DES_FB1.mfx')) as ff:
            for line in ff.readlines():
                kv = line.split('=')
                if len(kv) > 1 and kv[0].strip()=='description':
                    self.assertEqual(kv[1].strip()[1:-1], new_description)
                    found += 1

        self.assertEqual(found, 1)

    def test_tstop(self):
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
        with open(os.path.join(self.rundir,'DES_FB1.mfx')) as ff:
            for line in ff.readlines():
                kv = line.split('=')
                if len(kv) > 1 and kv[0].strip()=='tstop':
                    self.assertAlmostEqual(float(kv[1]), new_tstop)
                    found += 1

        self.assertEqual(found, 1)


    def test_run_disabled_no_exe(self):
        self.skipTest("FIXME, where did the run button go?")
        if self.find_exes():
            self.skipTest("Only valid when executables are not present")
        self.assertFalse(self.mfix.ui.run.spinbox_mfix_executables.isVisibleTo(self.mfix.ui.run))
        self.assertFalse(self.mfix.ui.run.resume_mfix_button.isEnabled())
        self.assertFalse(self.mfix.ui.run.run_mfix_button.isEnabled())
        self.assertFalse(self.mfix.ui.run.stop_mfix_button.isEnabled())
