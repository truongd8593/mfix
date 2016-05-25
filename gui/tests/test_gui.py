# -*- coding: utf-8 -*-
"""
run with "nosetests -v"
or
python -m unittest discover
"""

import fnmatch
import glob
import os
import time
#from xvfbwrapper import Xvfb

from qtpy.QtTest import QTest
from qtpy import QtCore
from qtpy import QtWidgets

import logging

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
        patterns = ('mfix', 'pymfix')
        for paths in [glob.glob(os.path.join(self.mfix_home, 'build', '*', 'build-aux', pattern)) for pattern in patterns]:
            for path in paths:
                matches.append(path)
        for path in os.environ['PATH'].split(os.pathsep):
            for pattern in patterns:
                exe = os.path.join(path, pattern)
                if os.path.exists(exe):
                    matches.append(exe)

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
        def get_project_file(cls):
            return getattr(cls, 'project_file', None)

        def set_project_file(cls, value):
            cls.project_file = value

        gui.MfixGui.get_project_file = get_project_file
        gui.MfixGui.set_project_file = set_project_file

        log = logging.getLogger()
        log.root.setLevel(logging.WARN)
        self.mfix = gui.MfixGui(self.qapp)
        self.mfix.show()
        QTest.qWaitForWindowShown(self.mfix)

        self.mfix.get_open_filename = lambda: mfix_dat
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
        newname = 'DES_FB1_new_name'
        newpath = os.path.join(self.rundir, '%s.mfx' % newname)

        self.mfix.get_save_filename = lambda: newpath
        QtCore.QTimer.singleShot(100, dismiss)
        QTest.mouseClick(self.mfix.ui.toolbutton_save_as, QtCore.Qt.LeftButton)

        self.assertEqual(newname, self.mfix.ui.general.lineedit_keyword_run_name.text())
        self.assertTrue(os.path.exists(newpath))

    def test_run(self):
        if not self.find_exes(): # We should be able to mock this for testing
            self.skipTest("Only valid when executables are present")

        # before running
        self.assertTrue(self.mfix.ui.run.spinbox_mfix_executables.isVisibleTo(self.mfix.ui.run))
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.text() == "Run")

        # start run
        QTest.mouseClick(self.mfix.ui.toolbutton_run_stop, QtCore.Qt.LeftButton)
        time.sleep(1)

        # during running
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.text() == "Stop")

        # stop run
        QTest.mouseClick(self.mfix.ui.run.button_run_stop_mfix, QtCore.Qt.LeftButton)
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.text() == "Resume")

        # start resume
        QTest.mouseClick(self.mfix.ui.run.button_run_stop_mfix, QtCore.Qt.LeftButton)
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.text() == "Stop")

        # stop mfix
        QTest.mouseClick(self.mfix.ui.run.button_run_stop_mfix, QtCore.Qt.LeftButton)

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
        QTest.keyClick(self.mfix.ui.general.combobox_keyword_description, QtCore.Qt.Key_Enter)
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
        if self.find_exes():
            self.skipTest("Only valid when executables are not present")
        self.assertFalse(self.mfix.ui.toolbutton_run_stop.isEnabled())
        self.assertFalse(self.mfix.ui.toolbutton_reset_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.run.button_run_stop_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.run.button_reset_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.run.button_pause_mfix.isEnabled())
