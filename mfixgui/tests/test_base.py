# -*- coding: utf-8 -*-
"""
run with "nosetests -v"
or
python -m unittest discover
"""

from __future__ import print_function, absolute_import, unicode_literals, division

from qtpy.QtTest import QTest
from qtpy import QtCore
import math

from .helper_functions import TestQApplication, waitFor, waitForWindow
from mfixgui.widgets import base
from mfixgui.project import Equation
from mfixgui.constants import *


class BaseWidgetTests(TestQApplication):
    ''' unit tests for the GUI '''

    def setUp(self):
        TestQApplication.setUp(self)

        PARAMETER_DICT['test_param'] = 10

    def tearDown(self):
        if hasattr(self, 'widget'):
            self.widget.close()

    def test_lineedit_str(self):
        self.widget = base.LineEdit()
        self.widget.setdtype(str)
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, 'test string')
        waitFor(10)
        self.assertEqual(self.widget.value, 'test string')

    def test_lineedit_float(self):
        self.widget = base.LineEdit()
        self.widget.setdtype('dp')
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, '10.0')
        waitFor(10)
        self.assertEqual(self.widget.value, 10.0)

    def test_lineedit_float_eq(self):
        self.widget = base.LineEdit()
        self.widget.setdtype('dp')
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, '10.0*4')
        waitFor(10)
        self.assertIsInstance(self.widget.value, Equation)
        self.assertEqual(self.widget.value.eq, '10.0*4')

    def test_lineedit_float_eq_pi(self):
        self.widget = base.LineEdit()
        self.widget.setdtype('dp')
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, 'pi')
        waitFor(10)
        self.assertIsInstance(self.widget.value, Equation)
        self.assertEqual(self.widget.value.eq, 'pi')
        self.assertEqual(float(self.widget.value), math.pi)

    def test_lineedit_float_eq_param(self):
        self.widget = base.LineEdit()
        self.widget.setdtype('dp')
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, 'test_param')
        waitFor(10)
        self.assertIsInstance(self.widget.value, Equation)
        self.assertEqual(self.widget.value.eq, 'test_param')

#    def test_lineedit_float_range(self):
#        self.widget = base.LineEdit()
#        self.widget.setdtype('dp')
#        self.widget.show()
#        self.widget.setValInfo(1, 0)
#        waitForWindow(self.widget)
#        QTest.keyClicks(self.widget, '100')
#        waitFor(10)
#        with self.assertRaises(ValueError):
#            val = self.widget.value
#
#        QTest.keyClicks(self.widget, '-1')
#        waitFor(10)
#        with self.assertRaises(ValueError):
#            val = self.widget.value

    def test_lineedit_int(self):
        self.widget = base.LineEdit()
        self.widget.setdtype('i')
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, '10')
        waitFor(10)
        self.assertEqual(self.widget.value, 10)

    def test_lineedit_int_with_float(self):
        self.widget = base.LineEdit()
        self.widget.setdtype('i')
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, '10.0')
        waitFor(10)
        self.assertEqual(self.widget.value, 10)

    def test_lineedit_int_eq(self):
        self.widget = base.LineEdit()
        self.widget.setdtype('i')
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, '10.0*4')
        waitFor(10)
        self.assertIsInstance(self.widget.value, Equation)
        self.assertEqual(self.widget.value.eq, '10.0*4')

    def test_lineedit_int_eq_param(self):
        self.widget = base.LineEdit()
        self.widget.setdtype('i')
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, 'test_param')
        waitFor(10)
        self.assertIsInstance(self.widget.value, Equation)
        self.assertEqual(self.widget.value.eq, 'test_param')

    def test_lineedit_no_qcompleter(self):
        self.widget = base.LineEdit()
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, 'si')
        waitFor(10)
        self.assertFalse(self.widget._completer.popup().isVisible())

    def test_lineedit_qcompleter(self):
        self.widget = base.LineEdit()
        self.widget.allow_parameters = True
        self.widget.show()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, 'si')
        waitFor(10)
        self.assertTrue(self.widget._completer.popup().isVisible())

    def test_lineedit_qcompleter_completion_sin(self):
        self.widget = base.LineEdit()
        self.widget.allow_parameters = True
        self.widget.show()
        self.widget.setFocus()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, 'si')
        waitFor(10)
        QTest.keyClick(self.widget._completer.popup(), QtCore.Qt.Key_Enter)
        self.assertEqual(self.widget.text(), 'sin')

    def test_lineedit_qcompleter_completion_eq(self):
        self.widget = base.LineEdit()
        self.widget.allow_parameters = True
        self.widget.show()
        self.widget.setFocus()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, '10*3+4/co')
        waitFor(10)
        QTest.keyClick(self.widget._completer.popup(), QtCore.Qt.Key_Enter)
        self.assertEqual(self.widget.text(), '10*3+4/cos')

    def test_lineedit_qcompleter_completion_eq_middle(self):
        self.widget = base.LineEdit()
        self.widget.allow_parameters = True
        self.widget.show()
        self.widget.setFocus()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, '10*3+4/cos')
        waitFor(10)
        self.widget.setCursorPosition(5)
        QTest.keyClicks(self.widget, 'si')
        QTest.keyClick(self.widget._completer.popup(), QtCore.Qt.Key_Enter)
        self.assertEqual(self.widget.text(), '10*3+sin4/cos')

    def test_lineedit_create_parameter(self):
        self.widget = base.LineEdit()
        self.widget.allow_parameters = True
        self.widget.show()
        self.widget.setFocus()
        waitForWindow(self.widget)
        QTest.keyClicks(self.widget, '10.0')
        waitFor(10)
        self.widget.create_parameter()
        self.assertEqual(self.widget.text(), 'param')
