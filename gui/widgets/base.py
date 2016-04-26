# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

import re
from collections import OrderedDict

# optional imports
try:
    import numpy as np
except ImportError:
    np = None

try:
    import pandas as pd
except ImportError:
    pd = None

# local imports
from qtpy import QtWidgets, QtCore
from tools.mfixproject import Keyword, Equation
from tools.general import to_text_string


class CommonBase(QtCore.QObject):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self):
        self.key = None
        self.defaultValue = None

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
            self.defaultValue = val

        if self.defaultValue is not None:
            self.updateValue(self.key, self.defaultValue, args=None)


class LineEdit(QtWidgets.QLineEdit, CommonBase):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QLineEdit.__init__(self, parent)
        CommonBase.__init__(self)

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.emitUpdatedValue)
        self.timer.setSingleShot(True)

        self.textEdited.connect(self.textEditedEvent)

        self.dtype = str

        self.regex_expression = re.compile('@\(([0-9.eEpiPI\+\-/*\(\))]+)\)')
        self.regex_mathOp = re.compile('([eEpiPI\+\-/*\^\(\)]+)')

    @property
    def value(self):
        if len(str(self.text())) == 0:
            return ''
        if self.dtype == str:
            return str(self.text())
        elif self.dtype == float:
            if self.regex_mathOp.findall(str(self.text())):
                return Equation(self.text())
            else:
                return float(str(self.text())) # TODO: validate input - cgw
        elif self.dtype == int:
            return int(float(str(self.text())))

    def updateValue(self, key, newValue, args=None):

        if newValue is not None:
            if self.regex_expression.findall(str(newValue)):
                self.setText(self.regex_expression.findall(str(newValue))[0])
            else:
                self.setText(str(newValue).replace("'", '').replace('"', ''))
        else:
            self.setText('')

    def textEditedEvent(self, event):
        self.timer.stop()
        self.timer.start(100)

    def default(self, val=None):
        if val is not None:
            self.defaultValue = val

        if self.defaultValue is not None:
            self.updateValue(self.key, self.defaultValue, args=None)
        else:
            self.clear()


class CheckBox(QtWidgets.QCheckBox, CommonBase):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QCheckBox.__init__(self, parent)
        CommonBase.__init__(self)
        self.released.connect(self.emitUpdatedValue)

    @property
    def value(self):
        return bool(self.isChecked())

    def updateValue(self, key, newValue, args=None):
        if isinstance(newValue, Keyword):
            newValue = newValue.value

        if newValue and isinstance(newValue, bool):
            self.setChecked(newValue)
        else:
            self.setChecked(False)


class ComboBox(QtWidgets.QComboBox, CommonBase):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QComboBox.__init__(self, parent)
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
        if isinstance(newValue, Keyword):
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


class SpinBox(QtWidgets.QSpinBox, CommonBase):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QDoubleSpinBox.__init__(self, parent)
        CommonBase.__init__(self)

        self.valueChanged.connect(self.emitUpdatedValue)

        self.dtype = int

    def emitUpdatedValue(self):
        self.value_updated.emit(self, {self.key: self.value()}, None)

    def updateValue(self, key, newValue, args=None):
        if isinstance(newValue, Keyword):
            newValue = newValue.value

        self.setValue(int(newValue))

    def setValInfo(self, _max=None, _min=None, req=False):
        if _max:
            self.setMaximum(int(_max))
        if _min:
            self.setMinimum(int(_min))


class DoubleSpinBox(QtWidgets.QDoubleSpinBox, CommonBase):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QDoubleSpinBox.__init__(self, parent)
        CommonBase.__init__(self)

        self.valueChanged.connect(self.emitUpdatedValue)

        self.dtype = float

    def emitUpdatedValue(self):
        self.value_updated.emit(self, {self.key: self.value()}, None)

    def updateValue(self, key, newValue, args=None):
        if isinstance(newValue, Keyword):
            newValue = newValue.value

        self.setValue(float(newValue))

    def setValInfo(self, _max=None, _min=None, req=False):
        if _max:
            self.setMaximum(float(_max))
        if _min:
            self.setMinimum(float(_min))


# --- Table ---
class Table(QtWidgets.QTableView, CommonBase):
    '''
    A table view with built in dictionary and array table models. Custom
    delegates are also provided.

    Parameters
    ----------
    parent (QObject):
        parent of the widget (default None)
    dtype (str):
        the type of data, should be either `dict` for displaying
        `dict(dict())`, `list` for displaying `list(list())`, 2D
        `numpy.array`, or Pandas `DataFrame` (default None)
    columns (list):
        a list of column names, if `None`, hides column names (default [])
    rows (list):
        a list of row names, if `None`, hides row names  (default [])
    columnDelegate (dict):
        a dictionary describing delegates for editing cells, column wise
        (default {})
    rowDelegate (dict):
        a dictionary describing delegates for editing cells, row wise
        (default {})
    selectionBehavior (str):
        a string describing the selection behavior. Either 'row', 'col', or
        'cell' for row selection, column selection, or single cell selection,
        respectively (default 'cell')
    selectionMode (str):
        a string describing the selection mode. Either 'single', or 'multi' for
        single selection or multiple selections (default 'single')

    Signals
    -------
    value_changed:
        emits the current value everytime the value is changed
    lost_focus:
        emits when the widget has lost focus
    new_selection:
        emits the from and to indices of a slection change.
    '''
    value_changed = QtCore.Signal(object)
    lost_focus = QtCore.Signal(object)
    new_selection = QtCore.Signal(object, object)

    def __init__(self, parent=None, dtype=None, columns=[], rows=[],
                 column_delegate={}, row_delegate={}, selection_behavior='row',
                 multi_selection=False):

        QtWidgets.QTableView.__init__(self, parent)
        CommonBase.__init__(self)

        self.dtype = dtype
        self.columns = columns
        self.rows = rows
        self.block_selection_change_event = False

        self._setModel()

        if columns is None:
            self.horizontalHeader().hide()
        if rows is None:
            self.verticalHeader().hide()

        self.set_delegate(col=column_delegate, row=row_delegate)

        self.set_selection_model(selection_behavior, multi_selection)

    def _setModel(self):

        # remove old model
        oldModel = self.model()
        if oldModel:
            oldModel.deleteLater()

        # Setup model
        if self.dtype in [dict, OrderedDict]:
            model = DictTableModel(columns=self.columns, rows=self.rows)
        elif self.dtype in [list, tuple,
                            np.ndarray if np else None,
                            pd.DataFrame if pd else None]:
            model = ArrayTableModel(columns=self.columns, rows=self.rows)
        else:
            model = None

        if model is not None:
            QtWidgets.QTableView.setModel(self, model)
    #        self.model().modelReset.connect(self.hideRows)
            self.model().value_updated.connect(
                lambda: self.value_changed.emit(self.value))
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

        if multi == 'multi':
            self.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
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

    @property
    def value(self):
        if self.model():
            return self.model().datatable
        else:
            return None

    def set_value(self, value):
        if self.dtype != type(value):
            raise TypeError('Selected table model does not support type'
                            '{}'.format(type(value)))

        if self.model() is not None:
            self.model().update(value)

        # reset the selection
        self.block_selection_change_event = True
        select_model = self.selectionModel()
        for selection in self.selection:
            select_model.setCurrentIndex(
                selection,
                QtWidgets.QItemSelectionModel.Select)
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

        self.delegate = CustomDelegate(columnDict=col,
                                       rowDict=row)
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

        # build context menu
        self.menu = QtWidgets.QMenu(self)
        applyAction = QtWidgets.QAction('Apply to Column', self)
        applyAction.triggered.connect(self.apply_val_to_column)
        self.menu.addAction(applyAction)

        # popup context menu
        self.menu.popup(QtWidgets.QCursor.pos())

    def apply_val_to_column(self):
        i = self.selectionModel().selection().indexes()[-1]
        value = self.model().data(i, QtCore.Qt.EditRole)
        self.model().apply_to_column(i.column(), value)

    def current_row(self):
        i = self.selectionModel().selection().indexes()
        if i:
            return i[-1].row()
        else:
            return -1

    def current_column(self):
        i = self.selectionModel().selection().indexes()
        if i:
            return i[-1].column()
        else:
            return -1

    def clear(self):
        self.model().update({})


class CustomDelegate(QtWidgets.QStyledItemDelegate):
    def __init__(self, columnDict={}, rowDict={}):
        QtWidgets.QStyledItemDelegate.__init__(self)

        self.columnDict = columnDict
        self.rowDict = rowDict

    def setColumnWidgets(self, columnW):
        self.columnDict.update(columnW)

    def setRowWidgets(self, rowW):
        self.rowDict.update(rowW)

    def createEditor(self, parent, option, index):
        if self.columnDict and index.column() in self.columnDict:
            widgetData = self.columnDict[index.column()]
        elif self.rowDict and index.row() in self.rowDict:
            widgetData = self.rowDict[index.row()]
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
                editor.setDtype(widgetData['dtype'])

            if 'max' in widgetData:
                _max = widgetData['max']
            else:
                _max = None

            if 'min' in widgetData:
                _min = widgetData['min']
            else:
                _min = None

            if hasattr(editor, 'setRange'):
                editor.setRange(_min, _max)

        return editor

    def setEditorData(self, widget, index):
        index.model().blockUpdate = True
        value = index.model().data(index, QtCore.Qt.EditRole)

        if widget.dtype is None:
            widget.dtype = type(value)
        widget.setValue(value)

    def setModelData(self, widget, model, index):
        model.setData(index, widget.dtype(widget.value), QtCore.Qt.EditRole)
        model.blockUpdate = False

    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)

    def eventFilter(self, widget, event):

        if isinstance(widget, ComboBox):
            # print(event, event.type())
            # TODO: there is a bug here that doesn't return focus to the
            # tableview
            if event.type() == QtCore.QEvent.FocusOut:
                if widget.isPopup:
                    return True
                else:
                    widget.isPopup = True
                    return False
            else:
                return QtWidgets.QStyledItemDelegate.eventFilter(self,
                                                                 widget,
                                                                 event)

        else:
            return QtWidgets.QStyledItemDelegate.eventFilter(self,
                                                             widget,
                                                             event)

        return False


class DictTableModel(QtCore.QAbstractTableModel):
    '''
    Table model that handles dict(dict()).

    Parameters
    ----------
    columns (list):
        a list of column names (used as dict index), suggested but not
        required. if not provided, defualts to keys of the dictionary
        (default [])
    rows (list)
        a list of row names (used as dict index), suggested but not required.
        if not provided, defualts to keys of the dictionary (default [])
    parent (QObject):
        parent of the model (default None)
    '''
    value_updated = QtCore.Signal(object)

    def __init__(self, columns=[], rows=[], parent=None, ):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.datatable = {}
        self._columns = columns
        self._rows = rows
        self.blockUpdate = False
        self.update_rows = False

    def update(self, data):
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
                col=None, row=None):
        if role == QtCore.Qt.EditRole:
            if index is not None:
                i = index.row()
                j = index.column()
            elif col is not None and row is not None:
                i = row
                j = col

            if self._columns:
                j = self._columns[j]
            if self._rows:
                i = self._rows[i]

            if i not in self.datatable:
                self.datatable[i] = {}

            self.datatable[i][j] = value
            self.value_updated.emit(self.datatable)

    def apply_to_column(self, col, val):
        for i in range(self.rowCount()):
            self.setData(col=col, row=i, value=val)

    def rowCount(self, parent=QtCore.QModelIndex()):
        if self.datatable.keys():
            return len(self.datatable.keys())
        else:
            return 0

    def columnCount(self, parent=QtCore.QModelIndex()):
        maxKeys = 0

        for value in self.datatable.values():
            maxKeys = max(maxKeys, len(value))

        return maxKeys

    def data(self, index, role=QtCore.Qt.DisplayRole):
        i = index.row()
        j = index.column()

        if self._columns:
            j = self._columns[j]
        if self._rows:
            i = self._rows[i]

        if i in self.datatable and j in self.datatable[i]:
            value = self.datatable[i][j]
        else:
            value = None

        if role == QtCore.Qt.DisplayRole:
            if value is None:
                value = None
            else:
                value = str(value)
            return value
        elif role == QtCore.Qt.EditRole:
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
    '''
    Table model that handles the following data types:
        - list()
        - tuple()
        - list(list())
        - tuple(tuple())
        - list(tuple())
        - tuple(list())
        - 2D numpy.ndarray
        - pandas.DataFrame
    '''
    value_updated = QtCore.Signal(object)

    def __init__(self, columns=[], rows=[], parent=None, ):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.datatable = []
        self._columns = columns
        self._rows = rows
        self.blockUpdate = False

    def update(self, data):
        if not self.blockUpdate:
            self.beginResetModel()
            self.datatable = data
            self.endResetModel()

            if pd and isinstance(self.datatable, pd.DataFrame):
                self._columns = list(self.datatable.columns)
                self._rows = list(self.datatable.index)

    def setData(self, index=None, value=None, role=QtCore.Qt.EditRole,
                col=None, row=None):
        if role == QtCore.Qt.EditRole:
            if index is not None:
                i = index.row()
                j = index.column()
            elif col is not None and row is not None:
                i = row
                j = col

            if pd and isinstance(self.datatable, pd.DataFrame):
                j = self.datatable.columns[j]
                self.datatable[j][i] = value
            else:
                self.datatable[i][j] = value

            self.value_updated.emit(self.datatable)

    def applyToColumn(self, col, val):
        for i in range(self.rowCount()):
            self.setData(col=col, row=i, value=val)

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.datatable)

    def columnCount(self, parent=QtCore.QModelIndex()):
        if (isinstance(self.datatable, list) or
                isinstance(self.datatable, tuple)) and \
                len(self.datatable) > 0:
            if isinstance(self.datatable[0], list) or \
                    isinstance(self.datatable[0], tuple):
                return len(self.datatable[0])
            else:
                return 1

        elif np and isinstance(self.datatable, np.ndarray):
            return int(self.datatable.shape[1])

        elif pd and isinstance(self.datatable, pd.DataFrame):
            return len(self.datatable.columns)

        else:
            return 0

    def data(self, index, role=QtCore.Qt.DisplayRole):
        i = index.row()
        j = index.column()

        if pd and isinstance(self.datatable, pd.DataFrame):
            j = self.datatable.columns[j]
            value = self.datatable[j][i]
        else:
            value = self.datatable[i][j]

        if role == QtCore.Qt.DisplayRole:
            if value is None:
                value = None
            else:
                value = str(value)
            return value
        elif role == QtCore.Qt.EditRole:
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
        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable | \
            QtCore.Qt.ItemIsSelectable
