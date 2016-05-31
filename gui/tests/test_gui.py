# -*- coding: utf-8 -*-
"""
run with "nosetests -v"
or
python -m unittest discover
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import glob
import os
#from xvfbwrapper import Xvfb

from qtpy.QtTest import QTest
from qtpy import QtCore
from qtpy.QtCore import Qt
from qtpy import QtWidgets

import logging

from tools.general import to_unicode_from_fs

from .helper_functions import TestQApplication
import gui

def dismiss():
    ''' dismiss modal QMessageBox '''
    for widget in QtWidgets.QApplication.instance().topLevelWidgets():
        if isinstance(widget, QtWidgets.QMessageBox):
            button = widget.button(QtWidgets.QMessageBox.Ok)
            if not button:
                button = widget.escapeButton()
            QTest.mouseClick(button, Qt.LeftButton)


class MfixGuiTests(TestQApplication):
    ''' unit tests for the GUI '''

    def find_exes(self): # All of mine live in Texas
        """find all mfix and pymfix executables"""
        matches = []
        names = 'mfix', 'mfix.exe', 'pymfix', 'pymfix.exe'
        for path in os.environ['PATH'].split(os.pathsep):
            if path=='' or path==os.path.curdir:
                continue
            for name in names:
                exe = os.path.join(path, name)
                if os.path.exists(exe):
                    matches.append(exe)
        for name in ('mfix', 'mfix.exe', 'pymfix', 'pymfix.exe'):
            wildcard = os.path.join(self.mfix_home, 'build',  name) # See open_project in gui.py
            for paths in glob.glob(wildcard):
                matches.extend(paths)

        return matches

    def setUp(self):
        ''' open FluidBed_DES for testing '''

        #self.xvfb = Xvfb(width=1280, height=720)
        #self.addCleanup(self.xvfb.stop)
        #self.xvfb.start()

        self.rundir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.mfix_home = os.path.dirname(self.rundir)
        self.rundir = os.path.join(self.mfix_home, 'tutorials', 'FluidBed_DES')
        self.runname = 'DES_FB1'
        mfix_dat = os.path.join(self.rundir, 'mfix.dat')

        patterns = [
            '*.LOG', '*.OUT', '*.RES', '*.SP?', '%s*' % self.runname, '*.mfx', 'MFIX.STOP',
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
        QTest.mouseClick(self.mfix.ui.toolbutton_open, Qt.LeftButton)

        self.assertEqual(self.runname, self.mfix.ui.general.lineedit_keyword_run_name.text())
        mfxfile = os.path.join(self.rundir, '%s.mfx' % self.runname)
        self.assertTrue(os.path.exists(mfxfile))

    def tearDown(self):
        patterns = [
            '*.LOG', '*.OUT', '*.RES', '*.SP?', '%s*' % self.runname, '*.mfx', 'MFIX.STOP',
            '*.pvd', '*.vtp', 'VTU_FRAME_INDEX.TXT']
        for paths in [glob.glob(os.path.join(self.rundir, n)) for n in patterns]:
            for path in paths:
                try:
                    os.remove(path)
                except OSError:
                    pass

        TestQApplication.tearDown(self)

    def get_tree_item(self, name):
        flags = Qt.MatchFixedString | Qt.MatchRecursive
        clist = self.mfix.ui.treewidget_model_navigation.findItems(name, flags, 0)
        assert len(clist) == 1
        return clist[0]

    def open_tree_item(self, name):
        self.mfix.ui.treewidget_model_navigation.setCurrentItem(self.get_tree_item(name))

    def test_save_as(self):
        #http://stackoverflow.com/questions/16536286/qt-ui-testing-how-to-simulate-a-click-on-a-qmenubar-item-using-qtest
        newname = '%s_new_name' % self.runname
        newpath = os.path.join(self.rundir, '%s.mfx' % newname)

        self.mfix.get_save_filename = lambda: newpath
        QtCore.QTimer.singleShot(100, dismiss)
        self.mfix.ui.action_save_as.trigger()

        self.assertEqual(newname, self.mfix.ui.general.lineedit_keyword_run_name.text())
        self.assertTrue(os.path.exists(newpath))

    def test_run(self):
        if not self.find_exes():
            self.skipTest("Only valid when executables are present (install mfix!)")

        # before running
        self.assertTrue(self.mfix.ui.run.combobox_mfix_executables.isVisibleTo(self.mfix.ui.run))
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.isEnabled())
        self.assertEqual(self.mfix.ui.run.button_run_stop_mfix.text(), "Run")

        # start run
        QTest.mouseClick(self.mfix.ui.toolbutton_run_stop_mfix, Qt.LeftButton)
        QTest.qWait(1000)

        # during running
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.isEnabled())
        self.assertEqual(self.mfix.ui.run.button_run_stop_mfix.text(), "Stop")

        # stop run
        QtCore.QTimer.singleShot(100, dismiss)
        QtCore.QTimer.singleShot(1000, dismiss)
        QTest.mouseClick(self.mfix.ui.run.button_run_stop_mfix, Qt.LeftButton)
        QTest.qWait(1000)
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.isEnabled())
        self.assertEqual(self.mfix.ui.run.button_run_stop_mfix.text(), "Resume")

        # start resume
        QTest.mouseClick(self.mfix.ui.run.button_run_stop_mfix, Qt.LeftButton)
        self.assertTrue(self.mfix.ui.run.button_run_stop_mfix.isEnabled())
        QTest.qWait(1000)
        self.assertEqual(self.mfix.ui.run.button_run_stop_mfix.text(),"Stop")

        # stop mfix
        QTest.mouseClick(self.mfix.ui.run.button_run_stop_mfix, Qt.LeftButton)

        logfile = os.path.join(self.rundir, '%s.LOG' % self.runname)
        self.assertTrue(os.path.exists(logfile))

    def test_description_ascii(self):
        self.open_tree_item('run')
        cb = self.mfix.ui.general.combobox_keyword_description
        description = cb.value
        cb.setFocus()
        QTest.keyClick(cb, Qt.Key_Right)
        new_text = u'some new text'
        QTest.keyClicks(cb, new_text)
        QTest.keyClick(cb, Qt.Key_Enter)
        QTest.mouseClick(self.mfix.ui.toolbutton_save, Qt.LeftButton)

        new_description = description + new_text

        found = 0
        with open(os.path.join(self.rundir,'%s.mfx' % self.runname)) as ff:
            for line in ff.readlines():
                kv = line.split('=')
                if len(kv) > 1 and kv[0].strip()=='description':
                    self.assertEqual(kv[1].strip()[1:-1], new_description)
                    found += 1

        self.assertEqual(found, 1)

    def test_description_unicode(self):
        self.open_tree_item('run')
        cb = self.mfix.ui.general.combobox_keyword_description
        description = cb.value
        cb.setFocus()
        QTest.keyClick(cb, Qt.Key_Right)

        new_text = u'mañana'
        ### KeyClicks won't take a unicode string - need to input ñ separately
        QTest.keyClicks(cb, "ma")
        QTest.keyClick(cb, Qt.Key_Ntilde)
        # how to input lower-case ñ ?
        QTest.keyClicks(cb, "ana")
        QTest.keyClick(cb, Qt.Key_Enter)
        QTest.mouseClick(self.mfix.ui.toolbutton_save, Qt.LeftButton)

        new_description = description + new_text

        found = 0
        with open(os.path.join(self.rundir,'%s.mfx' % self.runname)) as ff:
            for line in ff.readlines():
                line = to_unicode_from_fs(line)
                kv = line.split('=')
                if len(kv) > 1 and kv[0].strip()=='description':
                    self.assertEqual(kv[1].strip()[1:-1].lower(), new_description.lower())
                    # seeing 'Ñ' on py2,  'ñ' on py3.
                    found += 1

        self.assertEqual(found, 1)


    def test_tstop(self):
        self.open_tree_item('run')
        dt = self.mfix.ui.run.doublespinbox_keyword_dt.value()
        new_tstop = 5*dt
        self.mfix.ui.run.doublespinbox_keyword_tstop.setValue(new_tstop)
        QTest.mouseClick(self.mfix.ui.toolbutton_save, Qt.LeftButton)

        found = 0
        with open(os.path.join(self.rundir,'%s.mfx' % self.runname)) as ff:
            for line in ff.readlines():
                kv = line.split('=')
                if len(kv) > 1 and kv[0].strip() == 'tstop':
                    self.assertAlmostEqual(float(kv[1]), new_tstop)
                    found += 1

        self.assertEqual(found, 1)


    def test_run_disabled_no_exe(self):
        saved_path = os.environ['PATH']
        os.environ['PATH'] = ''
        self.mfix.reset()

        self.assertFalse(bool(self.find_exes()))
        self.assertFalse(self.mfix.ui.toolbutton_run_stop_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.toolbutton_reset_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.run.button_run_stop_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.run.button_reset_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.run.button_pause_mfix.isEnabled())
        os.environ["PATH"] = saved_path

    # TODO:  write more tests: tests for kill job, select different mfix exe,
    # close window with job running, failure to start job, etc
