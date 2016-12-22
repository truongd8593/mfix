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
from mfixgui.gui import MfixGui

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
        self.mfix.confirm_write = lambda *args: True
        self.mfix.confirm_rename = lambda *args: True
        self.mfix.confirm_delete_files = lambda *args: True
        self.mfix.main_menu_animation_speed = 0
        self.mfix.animation_speed = 0
        QTest.mouseClick(self.mfix.ui.toolbutton_file, Qt.LeftButton)
        browse_item = self.mfix.ui.main_menu_list.item(2)
        browse_item_rect = self.mfix.ui.main_menu_list.visualItemRect(browse_item)
        waitFor(50)
        QTest.mouseClick(self.mfix.ui.main_menu_list.viewport(), Qt.LeftButton, Qt.NoModifier, browse_item_rect.center())
        waitFor(50)
        QTest.mouseClick(self.mfix.ui.main_menu_browse, Qt.LeftButton)

        # We will get a confirmer for auto-rename


        self.assertEqual(self.runname, self.mfix.project.get_value('run_name'))
        mfxfile = os.path.join(self.rundir, '%s.mfx' % self.runname)
        self.assertTrue(os.path.exists(mfxfile))
        self.mfix.save_project()
        self.assertTrue('*' not in self.mfix.windowTitle())


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
        self.mfix.save_as()

        waitFor(500) # FIXME wait for a transition, not a hardcoded delay time
        self.assertEqual(newname, self.mfix.project.get_value('run_name'))
        self.assertTrue(os.path.exists(newpath))

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

    def test_run_mfix(self):
        self.skipTest("test removed, see issues/186")
        # Need to write a new test.  The old test relied heavily
        # on the Run pane toolbuttons & their label text


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

# class MfixGuiTests2(MfixGuiTests):

    # def __init__(self, *args, **kwargs):
    #     super(MfixGuiTests, self).__init__(*args, **kwargs)
    #     self.rundir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    #     self.mfix_home = os.path.dirname(self.rundir)
    #     self.rundir = os.path.join(self.mfix_home, 'benchmarks', 'dem', 'mini-cfb')
    #     self.runname = 'MINI-CFB'
