#!/usr/bin/env python
"""Parameter Dialog and management of the PARAMETER_DICT constant"""

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict
import copy

from qtpy import QtWidgets, QtCore

from mfixgui.constants import (
    PARAMETER_DICT,
    SPECIAL_PARAMETERS,
)
from mfixgui.project import Equation
from mfixgui.regexes import RE_MATH
from mfixgui.tools.general import (
    get_icon,
    get_unique_string,
)
from mfixgui.tools.util import (
    parse_key_with_args,
)
from mfixgui.tools.simpleeval import DEFAULT_FUNCTIONS, DEFAULT_NAMES
from mfixgui.widgets.base import Table

PROTECTED_NAMES = list(DEFAULT_FUNCTIONS.keys()) + list(DEFAULT_NAMES.keys()) + SPECIAL_PARAMETERS

class ParameterDialog(QtWidgets.QDialog):

    def __init__(self, parent):
        QtWidgets.QDialog.__init__(self, parent)

        self.data_old = {}

        self.setWindowIcon(get_icon('mfix.png'))
        self.setWindowTitle('Parameters')

        self.grid_layout = QtWidgets.QGridLayout(self)
        self.grid_layout.setContentsMargins(5, 5, 5, 5)

        # --- tool buttons ---
        self.button_bar = QtWidgets.QWidget(self)
        self.button_bar_layout = QtWidgets.QHBoxLayout(self.button_bar)
        self.button_bar_layout.setContentsMargins(0, 0, 0, 0)
        self.button_bar.setLayout(self.button_bar_layout)
        self.button_bar.setGeometry(QtCore.QRect(0, 0, 300, 300))
        self.grid_layout.addWidget(self.button_bar, 0, 0)

        self.toolbutton_add = QtWidgets.QToolButton()
        self.toolbutton_add.pressed.connect(self.new_parameter)
        self.toolbutton_add.setIcon(get_icon('add.png'))

        self.toolbutton_remove = QtWidgets.QToolButton()
        self.toolbutton_remove.pressed.connect(self.remove_parameter)
        self.toolbutton_remove.setIcon(get_icon('remove.png'))
        self.toolbutton_remove.setEnabled(False)

        self.toolbutton_copy = QtWidgets.QToolButton()
        self.toolbutton_copy.pressed.connect(self.copy_parameter)
        self.toolbutton_copy.setIcon(get_icon('copy.png'))
        self.toolbutton_copy.setEnabled(False)

        for widget in [self.toolbutton_add, self.toolbutton_remove,
                       self.toolbutton_copy]:
            self.button_bar_layout.addWidget(widget)
            widget.setAutoRaise(True)

        self.button_bar_layout.addStretch()

        # --- table ---
        self.table = Table(
            parent=self,
            dtype=OrderedDict,
            selection_behavior='row',
            multi_selection='multi',
            columns=['parameter', 'type', 'value'],
            column_delegate={
                0: {
                    'widget': 'lineedit',
                },
                1: {
                    'widget': 'combobox',
                    'items':  ['integer', 'float', 'string'],
                },
                2: {
                    'widget': 'lineedit',
                },
            }
        )
        self.table.show_vertical_header(False)
        self.table.auto_update_rows(True)
        self.table.new_selection.connect(self.table_clicked)
        self.table.value_changed.connect(self.parameter_changed)

        self.grid_layout.addWidget(self.table, 1, 0)

        # --- apply/close ---
        self.pushbutton_close = QtWidgets.QPushButton('Close')
        self.pushbutton_close.pressed.connect(self.close)
        self.grid_layout.addWidget(self.pushbutton_close, 2, 0)

    def update_table(self, data):
        self.table.set_value(data)
        self.data_old = copy.deepcopy(data)

    def table_clicked(self):
        rows = self.table.current_rows()

        if rows:
            self.toolbutton_remove.setEnabled(True)
            self.toolbutton_copy.setEnabled(True)
        else:
            self.toolbutton_remove.setEnabled(False)
            self.toolbutton_copy.setEnabled(False)

    def new_parameter(self):
        data = self.table.value
        new_name = get_unique_string(
            'new', [val['parameter'] for val in data.values()])

        new_index = 0
        if len(data.keys()) > 0:
            new_index = max(data.keys())+1

        data[new_index] = {'parameter': new_name, 'type': 'float',
                           'value': 0.0}

        self.update_table(data)

    def remove_parameter(self):
        rows = self.table.current_rows()

        if rows:
            data = self.table.value
            cant_remove = []
            for row in rows:
                index = list(data.keys())[row]
                name = data[index]['parameter']
                if name in self.parent().project.parameter_key_map:
                    cant_remove.append(name)
                else:
                    data.pop(index)
            self.update_table(data)
            if cant_remove:
                txt = ('The following parameters are being used: <b>{}</b>.'
                       ' Please remove reference before deleting.')
                self.parent().message(title='Error',
                                      text=txt.format(', '.join(cant_remove)))

        self.toolbutton_remove.setDown(False)

    def copy_parameter(self):
        rows = self.table.current_rows()

        if rows:
            data = self.table.value

            for row in rows:
                index = list(data.keys())[row]
                name = data[index]['parameter']

                new_name = get_unique_string(
                    name, [val['parameter'] for val in data.values()])

                new_index = 0
                if len(data.keys()) > 0:
                    new_index = max(data.keys())+1

                data[new_index] = {'parameter': new_name,
                                   'type': copy.deepcopy(data[index]['type']),
                                   'value': copy.deepcopy(data[index]['value'])}
            self.update_table(data)

        self.toolbutton_remove.setDown(False)

    def load_parameters(self):
        new_param_dict = OrderedDict()
        for i, key in enumerate(PARAMETER_DICT.keys()):
            dtype = 'string'
            value = PARAMETER_DICT[key]
            if isinstance(value, float):
                dtype = 'float'
            elif isinstance(value, int):
                dtype = 'integer'
            elif isinstance(value, Equation):
                dtype = {int:'integer', float:'float'}[value.dtype]
                value = value.eq
            new_param_dict[i] = {'parameter': key, 'type': dtype,
                                 'value': value}

        self.update_table(new_param_dict)
        self.table.fit_to_contents()

    @property
    def parameters(self):
        param_dict = OrderedDict()
        data = self.table.value
        param_names = [val['parameter'] for val in data.values()]
        for key, value in data.items():
            par_value = str(value['value'])
            if value['type'] in ['float', 'integer']:
                if RE_MATH.search(par_value) or any(par in par_value for par in param_names):
                    par_value = Equation(par_value)
            elif value['type'] == 'float':
                par_value = float(value['value'])
            elif value['type'] == 'integer':
                par_value = int(value['value'])
            param_dict[value['parameter']] = par_value
        return param_dict

    def get_parameters(self):
        self.load_parameters()
        self.changed_parameters = set()
        self.exec_()
        self.update_parameter_dict(self.parameters)
        self.table.clear_selection()
        return self.changed_parameters

    def update_parameter_dict(self, data):
        PARAMETER_DICT.update(data)

        for key in set(PARAMETER_DICT.keys()) - set(data.keys()):
            PARAMETER_DICT.pop(key)

    def parameter_changed(self, row, col, value):
        """parameter changed"""
        data = self.table.value

        # check value
        if col == 'value':
            self.changed_parameters.add(data[row]['parameter'])
            old_value = self.data_old[row][col]
            new_value = self.check_value(value, old_value, self.data_old[row]['type'])
            data[row][col] = new_value
            self.update_table(data)

        # check name
        elif col == 'parameter':
            old_name = self.data_old[row][col]
            new_name = self.check_name(value, old_name)

            data[row][col] = new_name
            self.update_table(data)
            self.table.clear_selection()

            if new_name != old_name:
                self.change_parameter_name(old_name, new_name)

    def check_value(self, value, old_value, dtype):
        if dtype == 'float':
            try:
                value = float(value)
            except:
                self.parent().message(title='Error', text='The value: <b>{}</b> is not a valid float'.format(value))
                value = old_value
        elif dtype == 'integer':
            try:
                value = int(value)
            except:
                self.parent().message(title='Error', text='The value: <b>{}</b> is not a valid integer'.format(value))
                value = old_value
        return value

    def check_name(self, name, old_name):
        """check the parameter name"""
        data = self.table.value
        param_names = [val['parameter'] for val in data.values()]
        param_names.remove(name)

        if name in PROTECTED_NAMES or name in self.parent().keyword_doc:
            self.parent().message(title='Error', text='The parameter name: <b>{}</b> is protected and cannot be used'.format(name))
            return old_name
        elif name in param_names:
            return get_unique_string(name, param_names)
        else:
            return name

    def change_parameter_name(self, old_name, new_name):
        """a parameter was renamed, update equations"""

        # update project keywords
        proj = self.parent().project
        p_map = proj.parameter_key_map

        if old_name in p_map:
            for keyword in p_map[old_name]:
                key, args = parse_key_with_args(keyword)
                eq = proj.get_value(key, default=None, args=args)
                if eq is not None and isinstance(eq, Equation):
                    eq.eq = eq.eq.replace(old_name, new_name)

            p_map[new_name] = p_map.pop(old_name)
            self.changed_parameters.add(new_name)

        # update regions
        regions = self.parent().ui.regions
        p_map = regions.parameter_key_map
        data = regions.tablewidget_regions.value

        if old_name in p_map:
            for keyword in p_map[old_name]:
                name, key = keyword.split(',')
                if 'to' in key or 'from' in key:
                    item = key.split('_')
                    index = ['x', 'y', 'z'].index(item[1])
                    val = data[name][item[0]][index]
                elif 'filter' in key:
                    item = key.split('_')
                    index = ['x', 'y', 'z'].index(item[1])
                    val = data[name][item[0]][index]
                else:
                    val = data[name][key]

                val.eq = val.eq.replace(old_name, new_name)
            p_map[new_name] = p_map.pop(old_name)
            regions.update_region_parameters()

        # update geometry
        geo = self.parent().vtkwidget
        p_map = geo.parameter_key_map
        data = geo.geometrydict

        if old_name in p_map:
            for keyword in p_map[old_name]:
                name, key = keyword.split(',')
                val = data[name][key]

                val.eq = val.eq.replace(old_name, new_name)
            p_map[new_name] = p_map.pop(old_name)
            geo.selected_geometry_changed()
