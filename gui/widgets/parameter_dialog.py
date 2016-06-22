
from collections import OrderedDict
from qtpy import QtWidgets, QtCore

from base import Table
from tools.general import get_icon, get_unique_string


class ParameterDialog(QtWidgets.QDialog):

    def __init__(self, parent, param_dict):
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
#        self.toolbutton_remove.pressed.connect(self.reset_view)
        self.toolbutton_remove.setIcon(get_icon('remove.png'))

        for widget in [self.toolbutton_add, self.toolbutton_remove]:
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
                                 'dtype':  'dp'
                                 },
                             }
            )
        self.table.show_vertical_header(False)
        self.table.auto_update_rows(True)

        self.grid_layout.addWidget(self.table, 1, 0)

        new_param_dict = OrderedDict()
        for i, key in enumerate(param_dict.keys()):
            new_param_dict[i] = {'parameter': key, 'type': 'float',
                                 'value': param_dict[key]}

        self.table.set_value(new_param_dict)
        self.table.fit_to_contents()

        # --- apply/close ---

    def new_parameter(self):
        param = self.table.value
        new_name = get_unique_string(
            'new', [val['parameter'] for val in param.values()])

        param[len(param)] = {'parameter': new_name, 'type': 'float',
                             'value': 0.0}

        self.table.set_value(param)

    @property
    def parameters(self):
        param_dict = OrderedDict()
        for key, value in self.table.value.items():
            param_dict[value['parameter']] = value['value']
        return param_dict

    @staticmethod
    def get_parameters(parameters, parent=None):

        dialog = ParameterDialog(parent, parameters)
        result = dialog.exec_()
        return dialog.parameters
