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
from qtpy.QtCore import Qt, QTimer
from qtpy import QtWidgets

import logging
import errno

from mfixgui.tools.general import to_unicode_from_fs

from .helper_functions import TestQApplication, waitFor, waitForWindow
from mfixgui.mfixgui import MfixGui

class MfixGuiTests(TestQApplication):
    ''' unit tests for the GUI '''

    def __init__(self, *args, **kwargs):
        super(TestQApplication, self).__init__(*args, **kwargs)
        self.rundir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.mfix_home = os.path.dirname(self.rundir)
        self.rundir = os.path.join(self.mfix_home, 'tutorials', 'FluidBed_DES')
        self.runname = 'DES_FB1'

    def setUp(self):
        """open FluidBed_DES for testing"""
        #self.xvfb = Xvfb(width=1280, height=720)
        #self.addCleanup(self.xvfb.stop)
        #self.xvfb.start()
        log = logging.getLogger()
        log.root.setLevel(logging.INFO)

        mfix_dat = os.path.join(self.rundir, 'mfix.dat')

        patterns = ['*.LOG', '*.OUT', '*.RES', '*.SP?', self.runname+'*', '*.mfx', 'MFIX.STOP',
                    '*.pvd', '*.vtp', 'VTU_FRAME_INDEX.TXT']
        for paths in [glob.glob(os.path.join(self.rundir, pat)) for pat in patterns]:
            for path in paths:
                try:
                    os.remove(path)
                except OSError as e:
                    if e.errno != errno.ENOENT:
                        print(e)

        TestQApplication.setUp(self)
        def get_project_file(cls):
            return getattr(cls, 'project_file', None)

        def set_project_file(cls, value):
            cls.project_file = value

        MfixGui.get_project_file = get_project_file
        MfixGui.set_project_file = set_project_file

        self.mfix = MfixGui(self.qapp)
        self.mfix.show()
        self.assertTrue(waitForWindow(self.mfix), "main mfix app not open")

        self.mfix.get_open_filename = lambda: mfix_dat
        self.mfix.check_if_ok_to_write = lambda *args: True
        self.mfix.check_if_ok_to_rename = lambda *args: True
        self.mfix.check_if_ok_to_delete_files = lambda *args: True
        QTest.mouseClick(self.mfix.ui.toolbutton_open, Qt.LeftButton)

        # We will get a confirmer for auto-rename


        self.assertEqual(self.runname, self.mfix.project.get_value('run_name'))
        mfxfile = os.path.join(self.rundir, '%s.mfx' % self.runname)
        self.assertTrue(os.path.exists(mfxfile))
        self.mfix.save_project()
        self.assertTrue('*' not in self.mfix.windowTitle())

        # used to break out of waitFor state tests in test_run_mfix
        self.attempts = 0
        self.maxattempts = 100


    def tearDown(self):
        patterns = [
            '*.LOG', '*.OUT', '*.RES', '*.SP?', '%s*' % self.runname, '*.mfx', 'MFIX.STOP',
            '*.pvd', '*.vtp', 'VTU_FRAME_INDEX.TXT']
        for paths in [glob.glob(os.path.join(self.rundir, n)) for n in patterns]:
            for path in paths:
                try:
                    os.remove(path)
                except OSError as e:
                    if e.errno != errno.ENOENT:
                        print(e)

        # Forcing 'chk_batchq_end=True' set the 'unsaved' flag
        # so we get a popup - save file to avoid this
        self.mfix.save_project()

        TestQApplication.tearDown(self)

    def get_tree_item(self, name):
        flags = Qt.MatchFixedString | Qt.MatchRecursive
        clist = self.mfix.ui.treewidget_navigation.findItems(name, flags, 0)
        assert len(clist) == 1
        return clist[0]

    def open_tree_item(self, name):
        self.mfix.ui.treewidget_navigation.setCurrentItem(self.get_tree_item(name))

    def test_save_as(self):
        #http://stackoverflow.com/questions/16536286/qt-ui-testing-how-to-simulate-a-click-on-a-qmenubar-item-using-qtest
        newname = '%s_new_name' % self.runname
        newpath = os.path.join(self.rundir, '%s.mfx' % newname)

        self.mfix.get_save_filename = lambda: newpath
        self.mfix.ui.action_save_as.trigger()

        waitFor(500) # FIXME wait for a transition, not a hardcoded delay time
        self.assertEqual(newname, self.mfix.project.get_value('run_name'))
        self.assertTrue(os.path.exists(newpath))


    def test_run_mfix(self):

        self.skipTest("temporarily disabling to test the test harness")
        # avoid test hangs in the waitFor loops below
        def retry(msg, delay):
            print('retrying test of "%s" state' % msg)
            self.attempts += 1
            # timed waits don't work well when testing system is loaded
            if self.attempts >= self.maxattempts:
                self.fail('FAIL in %s test' % msg)
            waitFor(delay)

        #  FIXME:  The run dialog will get the exe from the ~/.config/MFIX file,
        #   need to control the QSettings for running tests
        mfix_exe = os.path.join(self.mfix_home, "mfix")

        # Don't run the test if default mfix binary doesn't exist
        if not (os.path.isfile(mfix_exe) and os.access(mfix_exe, os.X_OK)):
            self.skipTest("Only valid when %s is present" % mfix_exe)

        self.open_tree_item("run")
        while self.mfix.ui.run.button_run_mfix.text() != "Run":
            retry('runnable', 100)
        self.attempts = 0
        print("RUNNABLE")

        # Before running, button says 'Run'
        #self.assertTrue(cme.isVisibleTo(self.mfix.ui.run))
        self.assertEqual(self.mfix.ui.toolbutton_run_mfix.toolTip(), "Run MFIX")
        self.assertFalse(self.mfix.ui.run.button_stop_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.toolbutton_stop_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.run.button_run_mfix.isVisible())
        self.assertTrue(self.mfix.ui.toolbutton_run_mfix.isEnabled())

        # Open run dialog
        # FIXME: This will hang if run dialog doesn't find exe.
        # Need to dismiss the warning message box.
        QTest.mouseClick(self.mfix.ui.run.button_run_mfix, Qt.LeftButton)
        while not self.mfix.ui.run:
            retry('run dialog open', delay=100)
        self.attempts = 0

        rr = range(self.mfix.run_dialog.combobox_mfix_exe.count())
        self.assertTrue(mfix_exe in [self.mfix.run_dialog.combobox_mfix_exe.itemText(i) for i in rr])
        self.mfix.run_dialog.combobox_mfix_exe.setEditText(mfix_exe)
        self.mfix.run_dialog.handle_exe_change()
        self.mfix.run_dialog.mfix_exe = mfix_exe

        # Press OK in run dialog
        QTest.mouseClick(self.mfix.run_dialog.button_local_run, Qt.LeftButton)
        while self.mfix.ui.run.button_run_mfix.text() != "Unpause":
            retry('running/paused', delay=500)
        self.attempts = 0
        print("PAUSED")

        # Job is running/paused:
        #   run button disabled
        #   pause button disabled
        #   unpause button enabled
        #   stop button enabled
        #   reset button disabled
        self.assertTrue(self.mfix.ui.run.button_run_mfix.isVisible())
        self.assertTrue(self.mfix.ui.toolbutton_run_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.run.button_pause_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.toolbutton_pause_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.run.button_stop_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.toolbutton_stop_mfix.isEnabled())

        # Unpause
        QTest.mouseClick(self.mfix.ui.run.button_run_mfix, Qt.LeftButton)
        while self.mfix.ui.run.button_run_mfix.isEnabled():
            retry('unpaused', delay=500)
        self.attempts = 0
        print("UNPAUSED")

        # Job is running/UNpaused:
        #   run button disabled
        #   pause button enabled
        #   unpause button disabled
        #   stop button enabled
        #   reset button disabled
        self.assertTrue(self.mfix.ui.run.button_run_mfix.isVisible())
        self.assertFalse(self.mfix.ui.toolbutton_run_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.run.button_pause_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.toolbutton_pause_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.run.button_stop_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.toolbutton_stop_mfix.isEnabled())

        #self.assertEqual(self.mfix.ui.run.button_run_mfix.text(), "Pause")
        #self.assertEqual(self.mfix.ui.toolbutton_run_mfix.toolTip(), "Stop MFIX")

        # Stop run, button should say 'Resume'
        QTest.mouseClick(self.mfix.ui.run.button_stop_mfix, Qt.LeftButton)
        while self.mfix.ui.run.button_run_mfix.text() != "Resume":
            retry('stopped/resumable', delay=500)
        self.attempts = 0
        print("STOPPED")

        # Job is stopped/resumable:
        #   run button enabled with text 'Resume'
        #   pause button disabled
        #   unpause button invisible
        #   stop button disabled
        #   reset button enabled
        self.assertTrue(self.mfix.ui.run.button_run_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.toolbutton_run_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.run.button_stop_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.toolbutton_stop_mfix.isEnabled())

        self.assertEqual(self.mfix.ui.run.button_run_mfix.text(), "Resume")
        self.assertEqual(self.mfix.ui.toolbutton_run_mfix.toolTip(), "Resume previous MFIX run")

        # Resume
        # FIXME: isn't this the same ui state as running/paused?
        QTest.mouseClick(self.mfix.ui.run.button_run_mfix, Qt.LeftButton)
        waitFor(100)

        # Press OK in run dialog
        QTest.mouseClick(self.mfix.run_dialog.button_local_run, Qt.LeftButton)
        waitFor(5000)

        # Press OK to delete output files

        while self.mfix.ui.run.button_run_mfix.text() != "Unpause":
            retry('resumed/paused', delay=500)
        self.attempts = 0
        print("RESUMED")

        self.assertTrue(self.mfix.ui.run.button_run_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.toolbutton_run_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.run.button_pause_mfix.isEnabled())
        self.assertFalse(self.mfix.ui.toolbutton_pause_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.run.button_stop_mfix.isEnabled())
        self.assertTrue(self.mfix.ui.toolbutton_stop_mfix.isEnabled())
        waitFor(300)

        # Stop mfix, check for log.
        QTest.mouseClick(self.mfix.ui.run.button_stop_mfix, Qt.LeftButton)
        waitFor(100)

    def test_description_ascii(self):
        self.open_tree_item('run')
        cb = self.mfix.ui.model_setup.combobox_keyword_description
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
        cb = self.mfix.ui.model_setup.combobox_keyword_description
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
        dt = self.mfix.project.get_value('dt')
        new_tstop = 5*dt
        self.mfix.update_keyword('tstop', new_tstop)
        QTest.mouseClick(self.mfix.ui.toolbutton_save, Qt.LeftButton)

        found = 0
        with open(os.path.join(self.rundir,'%s.mfx' % self.runname)) as ff:
            for line in ff.readlines():
                kv = line.split('=')
                if len(kv) > 1 and kv[0].strip() == 'tstop':
                    self.assertAlmostEqual(float(kv[1]), new_tstop)
                    found += 1

        self.assertEqual(found, 1)


    #def test_force_kill(self):
    #def test_stop_mfix(self):
    # TODO:  write more tests: tests for kill job, select different mfix exe,
    # close window with job running, failure to start job, etc

class MfixGuiTests2(MfixGuiTests):

    def __init__(self, *args, **kwargs):
        super(MfixGuiTests, self).__init__(*args, **kwargs)
        self.rundir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.mfix_home = os.path.dirname(self.rundir)
        self.rundir = os.path.join(self.mfix_home, 'benchmarks', 'dem', 'mini-cfb')
        self.runname = 'MINI-CFB'
