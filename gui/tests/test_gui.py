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
from qtpy.QtCore import Qt, QTimer
from qtpy import QtWidgets

import logging
import errno

from tools.general import to_unicode_from_fs

from .helper_functions import TestQApplication
import gui

# TODO : replace all qWaits with a 'waitFor' function
# Note, qWaitForWindowShown is deprecated in qt5, which provides a better version,
#  including a timeout.

try:
    waitForWindow = QTest.qWaitForWindowActive
except NameError:
    waitForWindow = QTest.qWaitForWindowShown

class MfixGuiTests(TestQApplication):
    ''' unit tests for the GUI '''


    def click_ok(self):
        retry = 0
        while not (self.mfix.message_box and self.mfix.message_box.isVisible()) and retry < 100:
            QTest.qWait(10)
            retry += 1
        self.assertTrue(self.mfix.message_box and self.mfix.message_box.isVisible(), "message box not shown in 1s")
        button = self.mfix.message_box.button(QtWidgets.QMessageBox.Ok)
        if not button:
            button = self.mfix.message_box.escapeButton()
        QTest.mouseClick(button, Qt.LeftButton)

        retry = 0
        while (self.mfix.message_box and self.mfix.message_box.isVisible()) and retry < 100:
            QTest.qWait(10)
            retry += 1
        self.assertFalse(self.mfix.message_box.isVisible(), 'dialog box not closed within 1s')


    def setUp(self):
        """open FluidBed_DES for testing"""
        #self.xvfb = Xvfb(width=1280, height=720)
        #self.addCleanup(self.xvfb.stop)
        #self.xvfb.start()
        log = logging.getLogger()
        log.root.setLevel(logging.INFO)

        self.rundir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.mfix_home = os.path.dirname(self.rundir)
        self.rundir = os.path.join(self.mfix_home, 'tutorials', 'FluidBed_DES')
        self.runname = 'DES_FB1'
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

        gui.MfixGui.get_project_file = get_project_file
        gui.MfixGui.set_project_file = set_project_file

        self.mfix = gui.MfixGui(self.qapp)
        self.mfix.show()
        self.assertTrue(waitForWindow(self.mfix), "main mfix app not open")

        self.mfix.get_open_filename = lambda: mfix_dat
        QTimer.singleShot(500, self.click_ok)
        QTest.mouseClick(self.mfix.ui.toolbutton_open, Qt.LeftButton)

        # We will get a confirmer for auto-rename


        self.assertEqual(self.runname, self.mfix.project.run_name.value)
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

        QTest.qWait(500)
        self.assertEqual(newname, self.mfix.project.run_name.value)
        self.assertTrue(os.path.exists(newpath))


    def test_run_mfix(self):
        self.skipTest("TODO: run mfix test based on run_popup modal dialog inputs")
        #TODO: write similar test for pymfix
        mfix_exe = os.path.join(self.mfix_home, "mfix")
        cme = self.mfix.ui.run.combobox_mfix_exes
        if cme.findText(mfix_exe) < 0:
            self.skipTest("Only valid when %s is present" % mfix_exe)

        #  FIXME:  we're getting the exe from the ~/.config/MFIX file,
        #   need to control the QSettings for running tests, instead
        #   of doing this!
        cme.setCurrentText(mfix_exe)
        self.mfix.handle_select_exe()

        self.open_tree_item("run")

        runbuttons = (self.mfix.ui.run.button_run_mfix,
                      self.mfix.ui.toolbutton_run_mfix)
        stopbuttons = (self.mfix.ui.run.button_stop_mfix,
                       self.mfix.ui.toolbutton_stop_mfix)

        # For thoroughness, we could loop over stobutton[0] and stopbutton[1].
        #  But we trust that they are connected to the same slot

        # Before running, button says 'Run'
        self.assertTrue(cme.isVisibleTo(self.mfix.ui.run))
        self.assertTrue(all (b.isEnabled() for b in runbuttons))
        self.assertTrue(all (not b.isEnabled() for b in stopbuttons))
        self.assertEqual(runbuttons[0].text(), "Run")
        self.assertEqual(runbuttons[1].toolTip(), "Run MFIX")

        # Start run, run button disables, stop button enabled
        QTest.mouseClick(runbuttons[0], Qt.LeftButton)
        QTest.qWait(500)



        self.assertTrue(all (not b.isEnabled() for b in runbuttons))
        self.assertTrue(all (b.isEnabled() for b in stopbuttons))

        # No test for pause, which requires pymfix FIXME
        #self.assertEqual(runbuttons[0].text(), "Pause")
        #self.assertEqual(runbuttons[1].toolTip(), "Stop MFIX")

        # Stop run, button should say 'Resume'
        QTest.mouseClick(stopbuttons[0], Qt.LeftButton)
        QTest.qWait(100)
        self.assertTrue(all (b.isEnabled() for b in runbuttons))
        self.assertTrue(all (not b.isEnabled() for b in stopbuttons))
        self.assertEqual(runbuttons[0].text(), "Resume")
        self.assertEqual(runbuttons[1].toolTip(), "Resume previous MFIX run")

        # Resume
        QTest.mouseClick(runbuttons[0], Qt.LeftButton)
        QTest.qWait(100)
        self.assertTrue(all (not b.isEnabled() for b in runbuttons))
        self.assertTrue(all (b.isEnabled() for b in stopbuttons))
        QTest.qWait(300)

        # Stop mfix, check for log.
        QTest.mouseClick(stopbuttons[0], Qt.LeftButton)
        QTest.qWait(100)
        logfile = os.path.join(self.rundir, '%s.LOG' % self.runname)
        self.assertTrue(os.path.exists(logfile))

    def test_description_ascii(self):
        self.open_tree_item('run')
        cb = self.mfix.ui.model.combobox_keyword_description
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
        cb = self.mfix.ui.model.combobox_keyword_description
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


    #def test_force_kill(self):
    #def test_stop_mfix(self):
    # TODO:  write more tests: tests for kill job, select different mfix exe,
    # close window with job running, failure to start job, etc
