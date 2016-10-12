# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

import copy
from collections import OrderedDict
from qtpy import QtWidgets, QtCore, QtGui


# Note - some items moved from 'Widgets' to 'Core'
# depending on version of both qtpy and Qt.
# This cannot be determined just using the PYQT4
# and PYQT5 constants, hence this try/except, for
# portability.  With newer versions of qtpy this
# can be removed (symbols are in QtCore)

try:
    from qtpy.QtWidgets import QStringListModel
except ImportError:
    from qtpy.QtCore import QStringListModel

try:
    from qtpy.QtWidgets import QItemSelectionModel
except ImportError:
    from qtpy.QtCore import QItemSelectionModel


from tools.general import get_selected_row

import logging
log = logging.getLogger(__name__)

# optional imports
try:
    import numpy as np
except ImportError:
    np = None
    log.debug("can't import numpy")

try:
    import pandas as pd
except ImportError:
    pd = None
    log.debug("can't import pandas")

# local imports
from project import Keyword, Equation, FloatExp, make_FloatExp
from regexes import *
from constants import *
from tools.general import (to_text_string, get_icon, insert_append_action,
                           insert_append_separator, get_unique_string)
from tools.simpleeval import VALID_EXPRESION_NAMES
VALID_EXP_NAMES = VALID_EXPRESION_NAMES + SPECIAL_PARAMETERS


def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)

class RangeError(ValueError):
    pass

class BaseWidget(QtCore.QObject):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self):
        self.key = None
        self.default_value = None
        self.saved_value = None
        self.args = None
        self.min = None
        self.max = None
        self.required = None
        self.help_text = 'No help avaliable.'
        self.context_menu = None
        self.allow_parameters = False

    def extend_context_menu(self):
        menu = self.context_menu()
        first_default_action = menu.actions()

        if first_default_action:
            first_default_action = first_default_action[0]
        # help
        help_action = QtWidgets.QAction(
            get_icon('help.png'), 'Help', menu)
        help_action.triggered.connect(self.show_help_message)
        insert_append_action(menu, help_action, first_default_action)

        # create parameter
        if self.allow_parameters:
            create_param_action = QtWidgets.QAction(
                get_icon('functions.png'), 'Create Parameter', menu)
            create_param_action.triggered.connect(self.create_parameter)
            insert_append_action(menu, create_param_action, first_default_action)

        insert_append_separator(menu, first_default_action)

        return menu

    def contextMenuEvent(self, event):
        if self.context_menu is not None:
            menu = self.extend_context_menu()
            menu.exec_(event.globalPos())

    def show_help_message(self):
        message_box = QtWidgets.QMessageBox(self)
        if self.key is not None:
            key = ': ' + self.key
        else:
            key = ''
        message_box.setWindowTitle('Help' + key)
        message_box.setIcon(QtWidgets.QMessageBox.Information)

        # Text
        message_box.setText(self.help_text)

        message_box.addButton(QtWidgets.QMessageBox.Ok)
        message_box.exec_()

    def create_parameter(self):
        btn = QtWidgets.QMessageBox.Yes
        if isinstance(self.value, Equation):
            message_box = QtWidgets.QMessageBox(self)
            message_box.setWindowTitle('Warning')
            message_box.setIcon(QtWidgets.QMessageBox.Warning)

            # Text
            message_box.setText("Warning: Replace equation with parameter?")

            message_box.addButton(QtWidgets.QMessageBox.Yes)
            message_box.addButton(QtWidgets.QMessageBox.No)
            message_box.setDefaultButton(QtWidgets.QMessageBox.No)
            btn = message_box.exec_()

        if btn == QtWidgets.QMessageBox.Yes:
            name = get_unique_string('param', PARAMETER_DICT.keys())

            v = self.dtype(0)
            if self.value:
                v = self.dtype(self.value)
            PARAMETER_DICT[name] = v

            self.updateValue(None, name)
            self.emitUpdatedValue()

    def emitUpdatedValue(self):
        self.value_updated.emit(self, {self.key: self.value}, self.args)

    def setdtype(self, dtype=None):
        dtype = to_text_string(dtype).strip().lower()
        if dtype == 'i' or 'int' in dtype:
            self.dtype = int
        elif dtype in ('d', 'dp') or 'float' in dtype:
            self.dtype = float
        elif dtype == 'l' or  'bool' in dtype:
            self.dtype = bool
        elif dtype in ('c', 's') or 'str' in dtype:
            self.dtype = str
        else:
            raise TypeError(self.objectName(), dtype)

    def setValInfo(self, max=None, min=None, required=None):
        if max is not None:
            self.max = max
        if min is not None:
            self.min = min
        if required is not None:
            self.required = required

    def default(self, val=None):
        if val is not None:
            self.default_value = val

        if self.default_value is not None:
            self.updateValue(self.key, self.default_value, args=self.args)

        return self.default_value


class LineEdit(QtWidgets.QLineEdit, BaseWidget):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QLineEdit.__init__(self, parent)
        BaseWidget.__init__(self)
        self.textChanged.connect(self.mark_changed)
        self.editingFinished.connect(self.emitUpdatedValue)
        self.dtype = str
        self.text_changed_flag = False
        self.context_menu = self.createStandardContextMenu

        self._separators = ['*', '**', '/', '-', '+', ' ']
        self._completer_model = QStringListModel(sorted(list(PARAMETER_DICT.keys())))
        self._completer = QtWidgets.QCompleter()
        self._completer.setModel(self._completer_model)
        self._completer.setWidget(self)
        self._completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        self._completer.activated.connect(self._insertCompletion)
        self._keysToIgnore = [QtCore.Qt.Key_Enter,
                              QtCore.Qt.Key_Return,
                              QtCore.Qt.Key_Escape,
                              QtCore.Qt.Key_Tab] # shouldn't TAB take us to next field?

    @classmethod
    def report_value_error(self, text):
        print(text)

    def mark_changed(self):
        self.text_changed_flag = True

    def emitUpdatedValue(self):
        need_to_signal = self.text_changed_flag
        self.text_changed_flag = False
        if need_to_signal:
            value = self.value
            if value is not None:
                self.value_updated.emit(self, {self.key: value}, self.args)

    def check_range(self, val):
        # TODO instead of raising an exception, make this a boolean function
        if self.min is not None and val < self.min:
            raise RangeError("Value below allowed range:\n %s < %s" % (val, self.min))
        if self.max is not None and val > self.max:
            raise RangeError("Value above allowed range:\n %s > %s" % (val, self.max))

    @property
    def value(self):
        text = self.text().strip()
        parameters = VALID_EXP_NAMES + sorted(list(PARAMETER_DICT.keys()))

        if len(text) == 0:   # should we return None?
            return ''
        if self.dtype is str:
            return text
        elif self.dtype is float:
            if re_float.match(text) or re_int.match(text):
                try:
                    f = float(text)
                except ValueError as e: # Should not really happen, unless our regexes are bad
                    self.report_value_error(e)
                    return self.updateValue(None, self.saved_value)
                try:
                    self.check_range(f)
                    self.saved_value = f
                    return f
                except ValueError as e:
                    self.report_value_error(e)
                    return self.updateValue(None, self.saved_value)
            elif re_float_exp.match(text):
                try:
                    f = make_FloatExp(text)
                except ValueError as e:
                    self.report_value_error(e)
                    return self.updateValue(None, self.saved_value)
                try:
                    self.check_range(f)
                    self.saved_value = f
                    return f
                except ValueError as e:
                    self.report_value_error(e)
                    return self.updateValue(None, self.saved_value)
            elif re_math.search(text) or any(par in text for par in parameters):
                if text.startswith('@(') and text.endswith(')'):
                    text = text[2:-1]
                try:
                    eq = Equation(text, dtype=float)
                except ValueError as e:
                    self.report_value_error("Invalid equation: %s" % e)
                    return self.updateValue(None, self.saved_value)
                try:
                    f = float(eq)
                except ValueError as e:
                    self.report_value_error("Invalid equation: %s" % e)
                    return self.updateValue(None, self.saved_value)
                try:
                    self.check_range(f)
                    self.saved_value = eq
                    return eq
                except ValueError as e:
                    self.report_value_error(e)
                    return self.updateValue(None, self.saved_value)
            else:
                return self.updateValue(None, self.saved_value)

        elif self.dtype is int:
            if re_math.search(text) or any(par in text for par in parameters):
                # integer equations?  do we use this?
                try:
                    eq = Equation(text, dtype=int)
                except ValueError as e:
                    self.report_value_error("Invalid equation: %s" % e)
                    return self.updateValue(None, self.saved_value)
                try:
                    i = int(eq)
                except ValueError as e:
                    self.report_value_error("Invalid equation: %s" % e)
                    return self.updateValue(None, self.saved_value)
                try:
                    self.check_range(i)
                    self.saved_value = eq
                    return eq
                except ValueError as e:
                    self.report_value_error(e)
                    return self.updateValue(None, self.saved_value)
            else:
                try:
                    i = int(float(text))
                except ValueError as e:
                    self.report_value_error(e)
                    return self.updateValue(None, self.saved_value)
                try:
                    self.check_range(i)
                    self.saved_value = i
                    return i
                except ValueError as e:
                    self.report_value_error(e)
                    return self.updateValue(None, self.saved_value)
        else:
            raise TypeError(self.dtype)

    def updateValue(self, key, new_value, args=None):
        if new_value is None:
            self.setText('')
            self.saved_value = None
            return

        if new_value is not None:
            self.saved_value = new_value

        sval = to_text_string(new_value).strip()

        # Don't show @( ) for equations.  (Is this a good idea?)
        while sval.startswith("@(") and sval.endswith(")"):
            sval = sval[2:-1]

        self.setText(sval)
        return new_value

    def default(self, val=None):
        if BaseWidget.default(self,val) is None:
            self.clear()
            self.saved_value = None
            self.text_changed_flag = False

    #completer functions
    def _insertCompletion(self, completion):
        """
        This is the event handler for the QCompleter.activated(QString) signal,
        it is called when the user selects an item in the completer popup.
        """
        # TODO FIXME this auto-completes an empty string
        text_under = self.textUnderCursor()
        cur_text = self.text()
        i = self.cursorPosition()

        beg = cur_text[:i]
        end = cur_text[i:]

        if text_under:
            beg = rreplace(beg, text_under, completion, 1)
        else:
            beg += completion
        self.setText(beg+end)

        self.setCursorPosition(len(beg))

    def textUnderCursor(self):
        text = self.text()
        textUnderCursor = ''
        i = self.cursorPosition() - 1
        while i >= 0 and text[i] not in self._separators:
            textUnderCursor = text[i] + textUnderCursor
            i -= 1
        return textUnderCursor

    def keyPressEvent(self, event):
        if not self.allow_parameters:
            QtWidgets.QLineEdit.keyPressEvent(self, event)
            return

        if self._completer.popup().isVisible():
            if event.key() in self._keysToIgnore:
                event.ignore()
                return
        else:
            self._update_completion_list()
        QtWidgets.QLineEdit.keyPressEvent(self, event)
        completionPrefix = self.textUnderCursor()
        if completionPrefix != self._completer.completionPrefix():
            self._updateCompleterPopupItems(completionPrefix)
        if event.text() and completionPrefix:
            self._completer.complete()
        if not completionPrefix:
            self._update_completion_list(allow_blank=True)
            self._updateCompleterPopupItems('')

    def _update_completion_list(self, allow_blank=False):
        comp_list = copy.deepcopy(VALID_EXP_NAMES)
        if allow_blank:
            comp_list.insert(0, '')
        for key, value in PARAMETER_DICT.items():
            if key in comp_list:
                continue
            elif self.dtype == str and isinstance(value, str):
                comp_list.append(key)
            elif self.dtype in [int, float]:
                comp_list.append(key)
        comp_list.sort()

        self._completer_model.setStringList(comp_list)
        self._completer.setModel(self._completer_model)

    def _updateCompleterPopupItems(self, completionPrefix):
        """
        Filters the completer's popup items to only show items
        with the given prefix.
        """
        self._completer.setCompletionPrefix(completionPrefix)
        self._completer.popup().setCurrentIndex(
                self._completer.completionModel().index(0, 0))


class CheckBox(QtWidgets.QCheckBox, BaseWidget):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QCheckBox.__init__(self, parent)
        BaseWidget.__init__(self)
        # stateChanged:  called on both user interaction and programmatic change
        # clicked:  user interaction only
        self.clicked.connect(self.emitUpdatedValue)
        self.context_menu = QtWidgets.QMenu
        self.dtype = bool

    @property
    def value(self):
        if self.dtype == int:
            return int(self.isChecked())
        elif self.dtype == float:
            return float(self.isChecked())
        elif self.dtype == bool:
            return bool(self.isChecked())
        else:
            raise TypeError("Invalid dtype %s" % self.dtype)

    def updateValue(self, key, new_value, args=None):
        assert not isinstance(new_value, Keyword)  # value should not be keyword!
        self.setChecked(new_value)

    def default(self, val=None):
        if BaseWidget.default(self, val) is None:
            self.setChecked(False) #?


class ComboBox(QtWidgets.QComboBox, BaseWidget):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QComboBox.__init__(self, parent)
        BaseWidget.__init__(self)
        # activated: only on user setttings, not programmatic change
        self.activated.connect(self.emitUpdatedValue)
        #self.currentIndexChanged.connect(self.emitUpdatedValue)
        self.dtype = str
        self.is_pop_up = False
        self.context_menu = QtWidgets.QMenu

    @property
    def value(self):
        sval = to_text_string(self.currentText())
        if self.dtype == int: # Labeled string entries, just return the number
            return int(sval.split('-')[0].strip())
        elif self.dtype == bool:
            return sval.lower() == 'true'
        else:
            return sval

    def updateValue(self, key, new_value, args=None):
        assert not isinstance(new_value, Keyword) # value should not be kw
        if isinstance(new_value, int):
            self.setCurrentIndex(new_value)
        else:
            self.setCurrentText(new_value)

    def setCurrentText(self, new_value):
        for itm in range(self.count()):
            if self.dtype == str and to_text_string(new_value).lower() == to_text_string(self.itemText(itm)).lower():
                self.setCurrentIndex(itm)
                break
            elif self.dtype == int and int(new_value) == int(to_text_string(self.itemText(itm)).split('-')[0].strip()):
                self.setCurrentIndex(itm)
                break
        else:
            raise ValueError(new_value)

    def default(self, val=None):
        if BaseWidget.default(self,val) is None:
            self.setCurrentIndex(0) # ?


class SpinBox(QtWidgets.QSpinBox, BaseWidget):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QDoubleSpinBox.__init__(self, parent)
        BaseWidget.__init__(self)
        # Would be nice to distinguish user input from programmatic setting
        self.valueChanged.connect(self.emitUpdatedValue)
        self.dtype = int
        self.context_menu = QtWidgets.QMenu

    def emitUpdatedValue(self): # calls self.value() instead of using self.value
        self.value_updated.emit(self, {self.key: self.value()}, self.args)

    def updateValue(self, key, new_value, args=None):
        assert not isinstance(new_value, Keyword)
        self.setValue(int(new_value))

    def setValInfo(self, max=None, min=None, required=None):
        BaseWidget.setValInfo(self, max, min, required)
        if max:
            self.setMaximum(int(max))
        if min:
            self.setMinimum(int(min))

    def default(self, val=None):
        if BaseWidget.default(self,val) is None:
            self.setValue(0) #?

class DoubleSpinBox(QtWidgets.QDoubleSpinBox, BaseWidget):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QDoubleSpinBox.__init__(self, parent)
        BaseWidget.__init__(self)
        # Distinguish user input from programmatic setting
        #self.valueChanged.connect(self.emitUpdatedValue)
        self.editingFinished.connect(self.emitUpdatedValue)
        self.dtype = float
        self.context_menu = QtWidgets.QMenu

    def textFromValue(self, value):
        ret = repr(value)
        return ret

    def emitUpdatedValue(self):
        self.value_updated.emit(self, {self.key: self.value()}, self.args)

    def updateValue(self, key, new_value, args=None):
        assert not isinstance(new_value, Keyword)
        self.setValue(float(new_value))

    def setValInfo(self, max=None, min=None, required=None):
        BaseWidget.setValInfo(self, max, min, required)
        if max:
            self.setMaximum(float(max))
        if min:
            self.setMinimum(float(min))

    def default(self, val=None):
        if BaseWidget.default(self,val) is None:
            self.setValue(0.0) #?


# --- Table ---
class Table(QtWidgets.QTableView, BaseWidget):
    """A table view with built in dictionary and array table models. Custom
    delegates are also provided.

    Parameters
    ----------
    parent (QObject):
        parent of the widget (default None)
    dtype (type):
        the type of data, should be either `dict` for displaying
        `dict(dict())`, `list` for displaying `list(list())`, 2D
        `numpy.array`, or Pandas `DataFrame` (default None)
    columns (list):
        a list of column names, if `None`, hides column names (default [])
    rows (list):
        a list of row names, if `None`, hides row names  (default [])
    column_delegate (dict):
        a dictionary describing delegates for editing cells, column wise
        (default {})
    row_delegate (dict):
        a dictionary describing delegates for editing cells, row wise
        (default {})
    selection_behavior (str):
        a string describing the selection behavior. Either 'row', 'col', or
        'cell' for row selection, column selection, or single cell selection,
        respectively (default 'cell')
    multi_selection (bool):
        Either single selection (False), or multiple selections (True)
        (default False)

    Signals
    -------
    value_changed:
        emits the current value everytime the value is changed
    lost_focus:
        emits when the widget has lost focus
    new_selection:
        emits the from and to indices of a selection change. """

    value_changed = QtCore.Signal(object, object, object)
    lost_focus = QtCore.Signal(object)
    new_selection = QtCore.Signal(object, object)

    def __init__(self, parent=None, dtype=None, columns=[], rows=[],
                 column_delegate={}, row_delegate={}, selection_behavior='row',
                 multi_selection=False):

        QtWidgets.QTableView.__init__(self, parent)
        BaseWidget.__init__(self)

        self.dtype = dtype
        self.columns = columns
        self.rows = rows
        self.block_selection_change_event = False
        self.selection = []
        self.mouse_pos = None
        self._setModel()
        self.set_delegate(col=column_delegate, row=row_delegate)
        self.set_selection_model(selection_behavior, multi_selection)
        if columns is None:
            self.horizontalHeader().hide()
        if rows is None:
            self.verticalHeader().hide()

        # build context menu
        self.menu = QtWidgets.QMenu(self)
        applyAction = QtWidgets.QAction('Apply to Column', self)
        applyAction.triggered.connect(self.apply_val_to_column)
        self.menu.addAction(applyAction)

    def _setModel(self):

        # remove old model
        oldModel = self.model()
        if oldModel:
            oldModel.deleteLater()

        # Setup model
        if self.dtype in (dict, OrderedDict):
            model = DictTableModel(columns=self.columns, rows=self.rows)
        elif self.dtype in (list, tuple,
                            np.ndarray if np else None,
                            pd.DataFrame if pd else None):
            model = ArrayTableModel(columns=self.columns, rows=self.rows)
        else:
            model = None

        if model is not None:
            QtWidgets.QTableView.setModel(self, model)
    #        self.model().modelReset.connect(self.hideRows)
            self.model().value_updated.connect(self.value_changed.emit)
            self.model().modelAboutToBeReset.connect(self.save_selection)

            # Need a reference or it segfaults with pyside
            # http://stackoverflow.com/questions/19211430/pyside-segfault-when-using-qitemselectionmodel-with-qlistview
            selectModel = self.selectionModel()
            selectModel.selectionChanged.connect(self.selection_changed_event)

    def set_selection_model(self, behavior='row', multi=False):
        " set the selection model "

        if behavior == 'col' or behavior == 'column':
            self.setSelectionBehavior(
                QtWidgets.QAbstractItemView.SelectColumns)
        elif behavior == 'row':
            self.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        else:
            self.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectItems)

        if multi:
            self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        else:
            self.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)

        if self.selectionModel():
            selectionmodel = self.selectionModel()
            selectionmodel.selectionChanged.connect(
                self.selection_changed_event)

    def selection_changed_event(self):
        selection = self.selectionModel().selection().indexes()
        if not self.block_selection_change_event:
            if selection:
                # from
                from_ = selection[0]
                # to
                to = selection[-1]
                self.new_selection.emit([from_.row(), from_.column()],
                                        [to.row(), to.column()])
            else:
                self.new_selection.emit(None, None)

    def clear_selection(self):
        sel_model = self.selectionModel()
        sel_model.clearSelection()

    @property
    def value(self):
        if self.model():
            return self.model().datatable
        else:
            return None

    def set_value(self, value, block=None):
        if self.dtype != type(value):
            # FIXME, this is a problem if you try to mix dict and OrderedDict
            raise TypeError('Selected table model does not support type'
                            '{}'.format(type(value)))

        if self.model() is not None:
            self.model().update(value, block)

        # reset the selection
        self.block_selection_change_event = True
        select_model = self.selectionModel()
        for selection in self.selection:
            select_model.setCurrentIndex(
                selection,
                QItemSelectionModel.Select)
        self.block_selection_change_event = False

    def set_columns(self, cols):
        self.columns = cols
        self.model()._columns = cols

        if cols is None:
            self.horizontalHeader().hide()
        else:
            self.horizontalHeader().show()

    def set_rows(self, rows):
        self.rows = rows
        self.model()._rows = rows

        if rows is None:
            self.verticalHeader().hide()
        else:
            self.verticalHeader().show()

    @property
    def column_colors(self):
        return self.delegate.column_color_dict

    @column_colors.setter
    def column_colors(self, color_dict):
        self.delegate.column_color_dict = color_dict

    @property
    def row_colors(self):
        return self.delegate.row_color_dict

    @row_colors.setter
    def row_colors(self, color_dict):
        self.delegate.row_color_dict = color_dict

    def auto_update_rows(self, b):
        self.model().update_rows = b

    def show_horizontal_header(self, b):
        if b:
            self.horizontalHeader().show()
        else:
            self.horizontalHeader().hide()

    def show_vertical_header(self, b):
        if b:
            self.verticalHeader().show()
        else:
            self.verticalHeader().hide()

    def set_delegate(self, col, row):
        self.delegate = CustomDelegate(column_dict=col,
                                       row_dict=row)
        self.setItemDelegate(self.delegate)

    def fit_to_contents(self):
        self.resizeColumnsToContents()
        self.horizontalHeader().setStretchLastSection(True)
        self.resizeRowsToContents()

    def save_selection(self):
        self.selection = self.selectionModel().selection().indexes()

#    def hideRows(self):
#        data = self.model().datatable
#        if data:
#            for i in range(1, max(data.keys())+1):
#                if i not in data:
#                    self.setRowHidden(i-1, True)
#                else:
#                    self.setRowHidden(i-1, False)

    def contextMenuEvent(self, event):
        """Qt context menu over-ride"""
        self.mouse_pos = QtGui.QCursor.pos()
        # popup context menu
        self.menu.popup(self.mouse_pos)

    def get_clicked_cell(self):
        """get the cell of the contectMenuEvent, return row, col"""
        local_pos = self.viewport().mapFromGlobal(self.mouse_pos)
        return (self.rowAt(local_pos.y()), self.columnAt(local_pos.x()))

    def apply_val_to_column(self):
        row, column  = self.get_clicked_cell()
        value = self.model().data(col=column, row=row, role=QtCore.Qt.EditRole)
        self.model().apply_to_column(column, value)

    def current_rows(self):
        i = self.selectionModel().selectedRows()
        if i:
            return [ind.row() for ind in i]
        else:
            return []

    def current_columns(self):
        i = self.selectionModel().selectedColumns()
        if i:
            return [ind.column() for ind in i]
        else:
            return []

    def clear(self):
        self.model().update({}) # TODO: change based on dtype?

    def default(self):
        '''if there is a default value, set it, else clear'''
        if self.default_value is not None:
            self.model().update(copy.deepcopy(self.default_value))
        else:
            self.clear()


class CustomDelegate(QtWidgets.QStyledItemDelegate):
    def __init__(self, column_dict={}, row_dict={}, column_color_dict={},
                 row_color_dict={}):
        QtWidgets.QStyledItemDelegate.__init__(self)

        self.column_dict = column_dict
        self.row_dict = row_dict
        self.row_color_dict = row_color_dict
        self.column_color_dict = column_color_dict

    def set_column_widgets(self, column_dict):
        self.column_dict.update(column_dict)

    def set_row_widgets(self, row_dict):
        self.row_dict.update(row_dict)

    def createEditor(self, parent, option, index):
        if self.column_dict and index.column() in self.column_dict:
            widgetData = self.column_dict[index.column()]
        elif self.row_dict and index.row() in self.row_dict:
            widgetData = self.row_dict[index.row()]
        else:
            widgetData = {'widget': 'lineedit'}

        editor = None

        if widgetData['widget'] == 'combobox':
            editor = ComboBox(parent)
            if 'items' in widgetData:
                editor.addItems(widgetData['items'])
        elif widgetData['widget'] == 'checkbox':
            editor = CheckBox(parent)
        elif widgetData['widget'] == 'lineedit':
            editor = LineEdit(parent)
        elif widgetData['widget'] == 'spinbox':
            editor = SpinBox(parent)
        elif widgetData['widget'] == 'doublespinbox':
            editor = DoubleSpinBox(parent)

        if editor:
            if 'dtype' in widgetData:
                editor.setdtype(widgetData['dtype'])

            max = widgetData.get('max')
            min = widgetData.get('min')

            if hasattr(editor, 'setRange'):
                editor.setRange(min, max)

        return editor

    def setEditorData(self, widget, index):
        index.model().blockUpdate = True
        value = index.model().data(index, QtCore.Qt.EditRole)

        if widget.dtype is None:
            widget.dtype = type(value)
        widget.updateValue(None, value)

    def setModelData(self, widget, model, index):
        model.setData(index, widget.dtype(widget.value), QtCore.Qt.EditRole)
        model.blockUpdate = False

    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)
        # What is this for?

    def eventFilter(self, widget, event):
        if isinstance(widget, ComboBox):
            # print(event, event.type())
            # TODO: there is a bug here that doesn't return focus to the
            # tableview
            if event.type() == QtCore.QEvent.FocusOut:
                pop = widget.is_pop_up
                widget.is_pop_up = True  # When does this get cleared?
                return pop

            else:
                return QtWidgets.QStyledItemDelegate.eventFilter(
                    self, widget, event)

        else: # ???
            return QtWidgets.QStyledItemDelegate.eventFilter(
                self, widget, event)

    def paint(self, painter, option, index):

        if index.column() in self.column_dict and self.column_dict[index.column()]['widget'] == 'progressbar':

            progress = int(index.model().data(index, QtCore.Qt.EditRole))

            progressbar = QtWidgets.QStyleOptionProgressBar()
            progressbar.rect = option.rect
            progressbar.minimum = 0
            progressbar.maximum = 100
            progressbar.progress = progress
            progressbar.text = '{}%'.format(progress)
            progressbar.textVisible = True
            progressbar.textAlignment = QtCore.Qt.AlignCenter

            QtWidgets.QApplication.style().drawControl(
                QtWidgets.QStyle.CE_ProgressBar, progressbar, painter)
        else:
            QtWidgets.QStyledItemDelegate.paint(self, painter, option, index)


class DictTableModel(QtCore.QAbstractTableModel):
    """Table model that handles dict(dict()).

    Parameters
    ----------
    columns (list):
        a list of column names (used as dict index), suggested but not
        required. if not provided, defaults to keys of the dictionary
        (default [])
    rows (list)
        a list of row names (used as dict index), suggested but not required.
        if not provided, defualts to keys of the dictionary (default [])
    parent (QObject):
        parent of the model (default None) """

    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, columns=[], rows=[], parent=None, ):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.datatable = {}
        self._columns = columns
        self._rows = rows
        self.blockUpdate = False
        self.update_rows = False

    def update(self, data, block=None):
        if block is not None:
            self.blockUpdate = block

        if not self.blockUpdate and data is not None:
            self.beginResetModel()
            self.datatable = data
            self.endResetModel()

            if not self._columns:
                self._columns = []
                for value in self.datatable.values():
                    if len(value) > len(self._columns):
                        self._columns = value.keys()

            if not self._rows or self.update_rows:
                self._rows = list(self.datatable.keys())

    def setData(self, index=None, value=None, role=QtCore.Qt.EditRole,
                row=0, col=0):
        if role == QtCore.Qt.EditRole:
            if index is not None:
                row = index.row()
                col = index.column()

            if self._columns:
                col = self._columns[col]
            if self._rows:
                row = self._rows[row]

            if row not in self.datatable:
                self.datatable[row] = {}

            self.datatable[row][col] = value
            self.value_updated.emit(row, col, value)

    def apply_to_column(self, col, val):
        for i in range(self.rowCount()):
            self.setData(col=col, row=i,
                         value=copy.deepcopy(val))

    def rowCount(self, parent=QtCore.QModelIndex()):
        'return the row count'
        if self.datatable.keys():
            return len(self.datatable.keys())
        else:
            return 0

    def columnCount(self, parent=QtCore.QModelIndex()):
        'return the column count'

        if self._columns:
            cols = len(self._columns)
        else:
            cols = 0
            for value in self.datatable.values():
                cols = max(cols, len(value))

        return cols

    def data(self, index=None, role=QtCore.Qt.DisplayRole, row=0, col=0):
        if index is not None:
            row = index.row()
            col = index.column()

        if self._columns:
            col = self._columns[col]
        if self._rows:
            row = self._rows[row]

        if row in self.datatable and col in self.datatable[row]:
            value = self.datatable[row][col]
        else:
            value = None

        if role == QtCore.Qt.DisplayRole:
            if value is None:
                value = None
            elif isinstance(value, QtGui.QPixmap):
                value = None
            else:
                value = to_text_string(value)
            return value
        elif role == QtCore.Qt.EditRole:
            return value
        elif role == QtCore.Qt.BackgroundRole and hasattr(value, 'qcolor'):
            return value.qcolor
        elif role == QtCore.Qt.DecorationRole and isinstance(value,
                                                             QtGui.QPixmap):
            return value
        else:
            return None

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                if len(self._columns) > section:
                    return self._columns[section]
                else:
                    return section

            elif orientation == QtCore.Qt.Vertical:
                if self._rows and len(self._rows) > section:
                    return self._rows[section]
                else:
                    return section
        else:
            return None

    def flags(self, index):
        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable | \
            QtCore.Qt.ItemIsSelectable


class ArrayTableModel(QtCore.QAbstractTableModel):
    """Table model that handles the following data types:
        - list()
        - tuple()
        - list(list())
        - tuple(tuple())
        - list(tuple())
        - tuple(list())
        - 2D numpy.ndarray
        - pandas.DataFrame """

    value_updated = QtCore.Signal(object)

    def __init__(self, columns=[], rows=[], parent=None, ):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.datatable = []
        self._columns = columns
        self._rows = rows
        self.blockUpdate = False

    def update(self, data, block=None):
        if block is not None:
            self.blockUpdate = block
        if not self.blockUpdate:
            self.beginResetModel()
            self.datatable = data
            self.endResetModel()

            if pd and isinstance(self.datatable, pd.DataFrame):
                self._columns = list(self.datatable.columns)
                self._rows = list(self.datatable.index)

    def setData(self, index=None, value=None, role=QtCore.Qt.EditRole, col=0,
                row=0):
        if role == QtCore.Qt.EditRole:
            if index is not None:
                row = index.row()
                col = index.column()

            if pd and isinstance(self.datatable, pd.DataFrame):
                col = self.datatable.columns[col]
                self.datatable[col][row] = value
            else:
                self.datatable[row][col] = value

            self.value_updated.emit(self.datatable)

    def apply_to_column(self, col, val):
        for i in range(self.rowCount()):
            self.setData(col=col, row=i,
                         value=copy.deepcopy(val))

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.datatable)

    def columnCount(self, parent=QtCore.QModelIndex()):
        if (isinstance(self.datatable, (list, tuple)) and
                len(self.datatable) > 0):
            if isinstance(self.datatable[0], (list, tuple)):
                return len(self.datatable[0])
            else:
                return 1

        elif np and isinstance(self.datatable, np.ndarray):
            return int(self.datatable.shape[1])

        elif pd and isinstance(self.datatable, pd.DataFrame):
            return len(self.datatable.columns)

        else:
            return 0

    def data(self, index=None, role=QtCore.Qt.DisplayRole, row=0, col=0):
        if index is not None:
            row = index.row()
            col = index.column()

        if pd and isinstance(self.datatable, pd.DataFrame):
            col = self.datatable.columns[col]
            value = self.datatable[col][row]
        else:
            value = self.datatable[row][col]

        if role == QtCore.Qt.DisplayRole:
            if value is None:
                value = None
            else:
                value = to_text_string(value)
            return value
        elif role == QtCore.Qt.EditRole:
            return value
        elif role == QtCore.Qt.BackgroundRole and hasattr(value, 'qcolor'):
            return value.qcolor
        elif role == QtCore.Qt.DecorationRole and isinstance(value,
                                                             QtGui.QPixmap):
            return value
        else:
            return None

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                if self._columns and len(self._columns) > section:
                    return self._columns[section]
                else:
                    return section

            elif orientation == QtCore.Qt.Vertical:
                if self._rows and len(self._rows) > section:
                    return self._rows[section]
                else:
                    return section
        else:
            return None

    def flags(self, index):
        return (QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable |
                QtCore.Qt.ItemIsSelectable)
