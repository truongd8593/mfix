# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

# impots
from qtpy import QtWidgets, QtCore
import copy

# local imports
from .base import Table
from gui.project import Keyword


class LinearEquationTable(QtWidgets.QWidget):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

        self.setObjectName('LinearEquationTable')

        # variables
        self.name_to_key = {
            'Scheme':           'discretize',
            'Method':           'leq_method',
            'Tolerance':        'leq_tol',
            'Iterations':       'leq_it',
            'Sweep':            'leq_sweep',
            'Preconditioner':   'leq_pc',
            'Under Relaxation': 'ur_fac',
            }
        # Build inverse dictionary
        self.key_to_name = dict((v, k) for (k, v) in self.name_to_key.items())

        # Move to constants.py
        self.rows = ['Gas Pressure', 'Volume Fraction', 'U', 'V', 'W', 'Energy',
                     'Mass Fraction', 'Granular  Temp', 'K-ε, Scalar', 'Diffusion']
        (GAS_PRESSURE, VOLUME_FRACTION, U, V, W, ENERGY,
         MASS_FRACTION, GRANULAR_TEMP, K_E, DIFFUSION) = self.rows
        # build default dictionary
        self.solverdict = {}
        for key in self.rows:
            self.solverdict[key] = {
                'Scheme':          0,
                'Method':          2,
                'Tolerance':       1E-4,
                'Sweep':           'RSRS',
                'Preconditioner':  'LINE',
                'Iterations':       15,
                'Under Relaxation': 0.8,
                }

        self.solverdict[GAS_PRESSURE]['Iterations'] = 20
        self.solverdict[VOLUME_FRACTION]['Iterations'] = 20
        self.solverdict[VOLUME_FRACTION]['Under Relaxation'] = 0.5
        self.solverdict[MASS_FRACTION]['Under Relaxation'] = 1.0
        self.solverdict[GRANULAR_TEMP]['Under Relaxation'] = 0.5

        for key in [U, V, W]:
            self.solverdict[key]['Iterations'] = 5
            self.solverdict[key]['Under Relaxation'] = 0.5

        # Solver Table
        self.table = Table(
            parent=self,
            dtype=dict,
            selection_behavior='cell',
            columns=['Scheme', 'Method', 'Tolerance', 'Iterations', 'Sweep',
                     'Preconditioner', 'Under Relaxation'],
            rows=self.rows,
            column_delegate={0: {'widget': 'combobox',
                                 'items':  ['0 - First-order upwinding',
                                            '1 - First-order upwinding (down-wind factors)',
                                            '2 - Superbee',
                                            '3 - Smart',
                                            '4 - Ultra-Quick',
                                            '5 - Quickest',
                                            '6 - MUSCL',
                                            '7 - van Leer',
                                            '8 - minmod',
                                            '9 - Central'],
                                 'dtype':  'i',
                                 },
                             1: {'widget': 'combobox',
                                 'items':  ['0 - SOR', '2 - BiCGSTAB',
                                            '3 - GMRES', '5 - CG'],
                                 'dtype':  'i',
                                 },
                             2: {'widget': 'lineedit',
                                 'dtype':  'dp'
                                 },
                             3: {'widget': 'lineedit',
                                 'dtype':  'i'},
                             4: {'widget': 'combobox',
                                 'items':  ['RSRS', 'ISIS', 'JSJS', 'KSKS',
                                            'ASAS']
                                 },
                             5: {'widget': 'combobox',
                                 'items':  ['NONE', 'LINE', 'DIAG']
                                 },
                             6: {'widget': 'lineedit',
                                 'dtype':  'dp',
                                 'max':    1.0,
                                 'min':    0.0
                                 },
                             }
            )

        self.table.value_changed.connect(self.model_edited)

        # set the default value
        self.table.default_value = copy.deepcopy(self.solverdict)

        self.grid_layout = QtWidgets.QGridLayout(self)
        self.grid_layout.addWidget(self.table, 0, 0)
        self.grid_layout.setSpacing(0)
        self.grid_layout.setContentsMargins(0, 0, 0, 0)

        self.table.set_value(self.solverdict)
        self.table.fit_to_contents()

    def model_edited(self, row, column, value):

        key = self.name_to_key[column]
        self.value_updated.emit(self,
                                {key: value},
                                # 1-based indexing
                                1 + self.rows.index(row))

    def updateValue(self, key, value, args):
        if isinstance(value, Keyword):
            value = value.value
        name = self.key_to_name[key]
        row_num = args[0] - 1         # 1-based indexing
        self.solverdict[self.rows[row_num]][name] = value
        self.table.update()
