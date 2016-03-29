# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

import re

# local imports
from qtpy import QtGui, QtCore
from tools.mfixproject import KeyWord, Equation
from tools.general import to_text_string


class CommonBase(QtCore.QObject):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self):
        self.key = None
        self.defualtValue = None

    def emitUpdatedValue(self):
        self.value_updated.emit(self, {self.key: self.value}, None)

    def validator(self):
        return False

    def updateValue(self, key, newValue, args=None):
        pass

    @property
    def value(self):
        return None

    def setdtype(self, dtype=None):
        dtype = to_text_string(dtype).lower().strip()
        if dtype == to_text_string('i'):
            self.dtype = int
        elif dtype == to_text_string('dp'):
            self.dtype = float
        elif dtype == to_text_string('bool'):
            self.dtype = bool
        else:
            self.dtype = str

    def setValInfo(self, _max=None, _min=None, req=False):
        self._max = _max
        self._min = _min
        self.required = req

    def default(self, val=None):
        if val is not None:
            self.defualtValue = val

        if self.defualtValue is not None:
            self.updateValue(self.key, self.defualtValue, args=None)


class LineEdit(QtGui.QLineEdit, CommonBase):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtGui.QLineEdit.__init__(self, parent)
        CommonBase.__init__(self)

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.emitUpdatedValue)
        self.timer.setSingleShot(True)

        self.textEdited.connect(self.textEditedEvent)

        self.dtype = str

        self.regX_expression = re.compile('@\(([0-9.eEpiPI\+\-/*\(\))]+)\)')
        self.regX_mathOp = re.compile('([eEpiPI\+\-/*\^\(\)]+)')

    @property
    def value(self):
        if len(str(self.text())) == 0:
            return ''
        if self.dtype == str:
            return str(self.text())
        elif self.dtype == float:
            if self.regX_mathOp.findall(str(self.text())):
                return Equation(self.text())
            else:
                return float(str(self.text()))
        elif self.dtype == int:
            return int(float(str(self.text())))

    def updateValue(self, key, newValue, args=None):

        if newValue:
            if self.regX_expression.findall(str(newValue)):
                self.setText(self.regX_expression.findall(str(newValue))[0])
            else:
                self.setText(str(newValue).replace("'", '').replace('"', ''))
        else:
            self.setText('')

    def textEditedEvent(self, event):
        self.timer.stop()
        self.timer.start(100)

    def default(self, val=None):
        if val is not None:
            self.defualtValue = val

        if self.defualtValue is not None:
            self.updateValue(self.key, self.defualtValue, args=None)
        else:
            self.clear()


class CheckBox(QtGui.QCheckBox, CommonBase):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtGui.QCheckBox.__init__(self, parent)
        CommonBase.__init__(self)
        self.released.connect(self.emitUpdatedValue)

    @property
    def value(self):
        return bool(self.isChecked())

    def updateValue(self, key, newValue, args=None):
        if isinstance(newValue, KeyWord):
            newValue = newValue.value

        if newValue and isinstance(newValue, bool):
            self.setChecked(newValue)
        else:
            self.setChecked(False)


class ComboBox(QtGui.QComboBox, CommonBase):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtGui.QComboBox.__init__(self, parent)
        CommonBase.__init__(self)
        self.activated.connect(self.emitUpdatedValue)

        self.dtype = str

    @property
    def value(self):
        if self.dtype == int:
            return int(str(self.currentText()).split('-')[0].strip())
        elif self.dtype == bool:
            if str(self.currentText()).lower() == 'true':
                return True
            else:
                return False
        else:
            return str(self.currentText())

    def updateValue(self, key, newValue, args=None):
        if isinstance(newValue, KeyWord):
            newValue = newValue.value

        self.setCurrentText(newValue)

    def setCurrentText(self, newValue):
        for itm in range(self.count()):
            if self.dtype == str and str(newValue).lower() == str(self.itemText(itm)).lower():
                self.setCurrentIndex(itm)
                break
            elif self.dtype == int and int(newValue) == int(str(self.itemText(itm)).split('-')[0].strip()):
                self.setCurrentIndex(itm)
                break


class SpinBox(QtGui.QSpinBox, CommonBase):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtGui.QDoubleSpinBox.__init__(self, parent)
        CommonBase.__init__(self)

        self.valueChanged.connect(self.emitUpdatedValue)

        self.dtype = int

    def emitUpdatedValue(self):
        self.value_updated.emit(self, {self.key: self.value()}, None)

    def updateValue(self, key, newValue, args=None):
        if isinstance(newValue, KeyWord):
            newValue = newValue.value

        self.setValue(int(newValue))

    def setValInfo(self, _max=None, _min=None, req=False):
        if _max:
            self.setMaximum(int(_max))
        if _min:
            self.setMinimum(int(_min))


class DoubleSpinBox(QtGui.QDoubleSpinBox, CommonBase):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtGui.QDoubleSpinBox.__init__(self, parent)
        CommonBase.__init__(self)

        self.valueChanged.connect(self.emitUpdatedValue)

        self.dtype = float

    def emitUpdatedValue(self):
        self.value_updated.emit(self, {self.key: self.value()}, None)

    def updateValue(self, key, newValue, args=None):
        if isinstance(newValue, KeyWord):
            newValue = newValue.value

        self.setValue(float(newValue))

    def setValInfo(self, _max=None, _min=None, req=False):
        if _max:
            self.setMaximum(float(_max))
        if _min:
            self.setMinimum(float(_min))
