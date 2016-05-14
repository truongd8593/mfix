# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

# local imports
from qtpy import QtWidgets, QtCore
from .base import Table
from project import Keyword


class LinearEquationTable(QtWidgets.QWidget):
    value_updated = QtCore.Signal(object, object, object)

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

        self.setObjectName('LinearEquationTable')

        # variables
        self.name_to_key = {
            'Method':           'leq_method',
            'Tolerance':        'leq_tol',
            'Iterations':       'leq_it',
            'Sweep':            'leq_sweep',
            'Preconditioner':   'leq_pc',
            'Under Relaxation': 'ur_fac',
            }
        # Build inverse dictionary
        self.key_to_name = dict((v, k) for (k, v) in self.name_to_key.items())

        self.rows = ['Gas Press', 'Vol Frac', 'U', 'V', 'W', 'Energy',
                     'Mass Frac', 'Gran. T', 'K-Ep, Scl.']

        self.solverdict = {
            'Gas Press': {'Method':          2,
                          'Tolerance':       1E-4,
                          'Iterations':      20,
                          'Sweep':           'RSRS',
                          'Preconditioner':  'LINE',
                          'Under Relaxation': 0.8,
                          },
            'Vol Frac': {'Method':           2,
                         'Tolerance':        1E-4,
                         'Iterations':       20,
                         'Sweep':            'RSRS',
                         'Preconditioner':   'LINE',
                         'Under Relaxation': 0.5,
                         },
            'U': {'Method':           2,
                  'Tolerance':        1E-4,
                  'Iterations':       5,
                  'Sweep':            'RSRS',
                  'Preconditioner':   'LINE',
                  'Under Relaxation': 0.5,
                  },
            'V': {'Method':           2,
                  'Tolerance':        1E-4,
                  'Iterations':       5,
                  'Sweep':            'RSRS',
                  'Preconditioner':   'LINE',
                  'Under Relaxation': 0.5,
                  },
            'W': {'Method':           2,
                  'Tolerance':        1E-4,
                  'Iterations':       5,
                  'Sweep':            'RSRS',
                  'Preconditioner':   'LINE',
                  'Under Relaxation': 0.5,
                  },
            'Energy': {'Method':           2,
                       'Tolerance':        1E-4,
                       'Iterations':       15,
                       'Sweep':            'RSRS',
                       'Preconditioner':   'LINE',
                       'Under Relaxation': 0.8,
                       },
            'Mass Frac': {'Method':           2,
                          'Tolerance':        1E-4,
                          'Iterations':       15,
                          'Sweep':            'RSRS',
                          'Preconditioner':   'LINE',
                          'Under Relaxation': 1.0,
                          },
            'Gran. T': {'Method':           2,
                        'Tolerance':        1E-4,
                        'Iterations':       15,
                        'Sweep':            'RSRS',
                        'Preconditioner':   'LINE',
                        'Under Relaxation': 0.5,
                        },
            'K-Ep, Scl.': {'Method':           2,
                           'Tolerance':        1E-4,
                           'Iterations':       15,
                           'Sweep':            'RSRS',
                           'Preconditioner':   'LINE',
                           'Under Relaxation': 0.8,
                           },
            }

        # Solver Table
        self.table = Table(
            parent=self,
            dtype=dict,
            selection_behavior='cell',
            columns=['Method', 'Tolerance', 'Iterations', 'Sweep',
                     'Preconditioner', 'Under Relaxation'],
            rows=self.rows,
            column_delegate={0: {'widget': 'combobox',
                                 'items':  ['0 - SOR', '2 - BiCGSTAB',
                                            '3 - GMRES', '5 - CG'],
                                 'dtype':  'i',
                                 },
                             1: {'widget': 'lineedit',
                                 'dtype':  'dp'
                                 },
                             2: {'widget': 'lineedit',
                                 'dtype':  'i'},
                             3: {'widget': 'combobox',
                                 'items':  ['RSRS', 'ISIS', 'JSJS', 'KSKS',
                                            'ASAS']
                                 },
                             4: {'widget': 'combobox',
                                 'items':  ['NONE', 'LINE', 'DIAG']
                                 },
                             5: {'widget': 'lineedit',
                                 'dtype':  'dp',
                                 'max':    1.0,
                                 'min':    0.0
                                 },
                             }
            )

        self.table.value_changed.connect(self.model_edited)

        self.grid_layout = QtWidgets.QGridLayout(self)
        self.grid_layout.addWidget(self.table, 0, 0)
        self.grid_layout.setSpacing(0)
        self.grid_layout.setContentsMargins(0, 0, 0, 0)

        self.table.set_value(self.solverdict)
        self.table.fit_to_contents()

    def model_edited(self, row, column, value):

        key = self.name_to_key[column]
        self.value_updated.emit(self,
                                {key: Keyword(key, value)},
                                # 1-based indexing
                                1 + self.rows.index(row))

    def updateValue(self, key, value, args):
        if isinstance(value, Keyword):
            value = value.value
        name = self.key_to_name[key]
        row_num = args[0] - 1         # 1-based indexing
        self.solverdict[self.rows[row_num]][name] = value
        self.table.update()
