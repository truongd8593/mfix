from collections import OrderedDict
import json
from qtpy import QtWidgets, QtCore
import copy

from base import Table
from tools.general import get_icon, get_unique_string
from constants import *

class ParameterDialog(QtWidgets.QDialog):

    def __init__(self, parent):
        QtWidgets.QDialog.__init__(self, parent)

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
            selection_behavior='cell',
            columns=['parameter', 'type', 'value'],
            column_delegate={0: {'widget': 'lineedit',
                                 },
                             1: {'widget': 'combobox',
                                 'items':  ['integer', 'float', 'string'],
                                 },
                             2: {'widget': 'lineedit',
                                 },
                             }
            )
        self.table.show_vertical_header(False)
        self.table.auto_update_rows(True)
        self.table.new_selection.connect(self.table_clicked)

        self.grid_layout.addWidget(self.table, 1, 0)

        # --- apply/close ---
        self.pushbutton_close = QtWidgets.QPushButton('Close')
        self.pushbutton_close.pressed.connect(self.close)
        self.grid_layout.addWidget(self.pushbutton_close, 2, 0)

    def table_clicked(self):
        row = self.table.current_row()

        if row >= 0:
            self.toolbutton_remove.setEnabled(True)
            self.toolbutton_copy.setEnabled(True)
        else:
            self.toolbutton_remove.setEnabled(False)
            self.toolbutton_copy.setEnabled(False)

    def new_parameter(self):
        param = self.table.value
        new_name = get_unique_string(
            'new', [val['parameter'] for val in param.values()])

        param[len(param)] = {'parameter': new_name, 'type': 'float',
                             'value': 0.0}

        self.table.set_value(param)

    def remove_parameter(self):
        row = self.table.current_row()

        if row >= 0:
            data = self.table.value
            name = list(data.keys())[row]
            data.pop(name)
            self.table.set_value(data)

    def copy_parameter(self):
        row = self.table.current_row()

        if row >= 0:
            data = self.table.value
            name = data[row]['parameter']

            new_name = get_unique_string(
                name, [val['parameter'] for val in data.values()])

            data[len(data)] = {'parameter': new_name,
                               'type': copy.deepcopy(data[row]['type']),
                               'value': copy.deepcopy(data[row]['value'])}
            self.table.set_value(data)

    def load_parameters(self):
        new_param_dict = OrderedDict()
        for i, key in enumerate(PARAMETER_DICT.keys()):
            dtype = 'string'
            if isinstance(PARAMETER_DICT[key], float):
                dtype = 'float'
            elif isinstance(PARAMETER_DICT[key], int):
                dtype = 'integer'
            new_param_dict[i] = {'parameter': key, 'type': dtype,
                                 'value': PARAMETER_DICT[key]}

        self.table.set_value(new_param_dict)
        self.table.fit_to_contents()

    @property
    def parameters(self):
        param_dict = OrderedDict()
        for key, value in self.table.value.items():
            par_value = str(value['value'])
            if value['type'] == 'float':
                par_value = float(value['value'])
            elif value['type'] == 'integer':
                par_value = int(value['value'])
            param_dict[value['parameter']] = par_value
        return param_dict

    def get_parameters(self):
        self.load_parameters()
        self.exec_()
        self.update_parameter_dict(self.parameters)

    def parameters_from_str(self, string):
        """load parameter data from a saved string"""
        loaded_data = json.loads(string)

        if 'order' not in loaded_data:
            return

        data = OrderedDict()
        for par in loaded_data['order']:
            data[par] = loaded_data['parameters'][par]

        self.update_parameter_dict(data)

    def parameters_to_str(self):
        """convert parameter data to a string for saving"""
        data = {'order': list(PARAMETER_DICT.keys()),
                'parameters': PARAMETER_DICT
                }
        return json.dumps(data)

    def update_parameter_dict(self, data):
        PARAMETER_DICT.update(data)

        for key in set(PARAMETER_DICT.keys()) - set(data.keys()):
            PARAMETER_DICT.pop(key)