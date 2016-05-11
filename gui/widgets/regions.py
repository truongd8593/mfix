# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

import os
import copy
from collections import OrderedDict
from qtpy import QtWidgets

# TODO: add pyside? There is an issue to add this to qtpy:
# https://github.com/spyder-ide/qtpy/issues/16
# TODO: cache ui file to a py file, use the ui file if newer, else py file
try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic

# local imports
from tools.general import (get_unique_string, widget_iter, CellColor)


class RegionsWidget(QtWidgets.QWidget):

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)


        uifiles = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'uifiles')
        uic.loadUi(os.path.join(uifiles, 'regions.ui'), self)

        self.combobox_regions_type.addItems(['point', 'XY-plane',
                                             'XZ-plane', 'YZ-plane', 'box',
                                             'STL',
                                             ])

        self.toolbutton_region_add.pressed.connect(
            self.new_region)
        self.toolbutton_region_delete.pressed.connect(
            self.delete_region)
        self.toolbutton_region_copy.pressed.connect(
            self.copy_region)
        self.toolbutton_color.pressed.connect(self.change_color)

        tablewidget = self.tablewidget_regions
        tablewidget.dtype = OrderedDict
        tablewidget._setModel()
        tablewidget.set_columns(['visibility', 'color', 'type', 'from', 'to'])
        tablewidget.show_vertical_header(True)
        tablewidget.set_value(OrderedDict())
        tablewidget.auto_update_rows(True)
        tablewidget.set_selection_model('cell', multi=False)
        tablewidget.new_selection.connect(self.update_region_parameters)
#        tablewidget.value_changed.connect(self.region_value_changed)

        for widget in widget_iter(self.groupbox_region_parameters):
            if hasattr(widget, 'value_updated'):
                widget.value_updated.connect(self.region_value_changed)

                # example <name>: lineedit_regions_to_x
                name = str(widget.objectName())
                if '_to_' in name:
                    widget.key = '_'.join(name.split('_')[-2:])
                    widget.dtype = float
                elif '_from_' in str(widget.objectName()):
                    widget.key = '_'.join(name.split('_')[-2:])
                    widget.dtype = float
                elif 'name' in name:
                    widget.key = 'name'
                    widget.dtype = str
                elif 'type' in name:
                    widget.key = 'type'
                    widget.dtype = str
                elif 'stl_shape' in name:
                    widget.key = 'stl_shape'
                    widget.dtype = str
                elif 'slice' in name:
                    widget.key = 'slice'
                    widget.dtype = bool
                elif 'filter' in name:
                    widget.key = '_'.join(name.split('_')[-2:])
                    widget.dtype = float
                elif 'deviation_angle' in name:
                    widget.key = 'deviation_angle'
                    widget.dtype = float

    def new_region(self):
        'create a new region'

        data = self.tablewidget_regions.value

        name = get_unique_string('new', list(data.keys()))

        data[name] = {'type':            'box',
                      'from':            [0, 0, 0],
                      'to':              [0, 0, 0],
                      'color':           CellColor([1, 0, 0]),
                      'slice':           True,
                      'filter':          [0, 0, 0],
                      'stl_shape':       'box',
                      'deviation_angle': 10,
                      'visibility':      True,
                      }

        self.vtkwidget.new_region(name, data[name])

        self.tablewidget_regions.set_value(data)
        self.tablewidget_regions.fit_to_contents()

    def delete_region(self):
        'remove the currently selected region'

        row = self.tablewidget_regions.current_row()

        if row >= 0:
            data = self.tablewidget_regions.value
            name = list(data.keys())[row]
            data.pop(name)
            self.tablewidget_regions.set_value(data)
            self.vtkwidget.delete_region(name)

            if data:
                self.groupbox_region_parameters.setEnabled(True)
            else:
                self.groupbox_region_parameters.setEnabled(False)

    def copy_region(self):
        'copy the currently selected region'

        row = self.tablewidget_regions.current_row()

        if row >= 0:
            data = self.tablewidget_regions.value
            name = list(data.keys())[row]
            new_region = copy.deepcopy(data[name])

            new_name = get_unique_string(name, list(data.keys()))
            data[new_name] = new_region

            self.tablewidget_regions.set_value(data)

            self.vtkwidget.new_region(new_name, data[new_name])

    def update_region_parameters(self):
        'a new region was selected, update region widgets'

        row = self.tablewidget_regions.current_row()

        if row >= 0:
            # enable groupbox
            self.groupbox_region_parameters.setEnabled(True)

            data = self.tablewidget_regions.value
            name = list(data.keys())[row]

            # enable widgets
            self.enable_disable_widgets(name)

            # color
            self.toolbutton_color.setStyleSheet(
                "QToolButton{{ background: rgb({},{},{});}}".format(
                    *data[name]['color'].color_int))

            self.lineedit_regions_name.updateValue(None, name)
            self.combobox_regions_type.updateValue(
                None,
                data[name]['type'])

            for widget, value in zip([self.lineedit_regions_from_x,
                                      self.lineedit_regions_from_y,
                                      self.lineedit_regions_from_z],
                                     data[name]['from']
                                     ):
                widget.updateValue(None, value)

            for widget, value in zip([self.lineedit_regions_to_x,
                                      self.lineedit_regions_to_y,
                                      self.lineedit_regions_to_z],
                                     data[name]['to']
                                     ):
                widget.updateValue(None, value)

            # stl
            self.combobox_stl_shape.updateValue(
                None,
                data[name]['stl_shape'])

            self.checkbox_slice_facets.updateValue(
                None,
                data[name]['slice'])

            self.lineedit_deviation_angle.updateValue(
                None,
                data[name]['deviation_angle'])

            for widget, value in zip([self.lineedit_filter_x,
                                      self.lineedit_filter_y,
                                      self.lineedit_filter_z],
                                     data[name]['filter']
                                     ):
                widget.updateValue(None, value)

        else:
            self.groupbox_region_parameters.setEnabled(False)

    def region_value_changed(self, widget, value, args):
        'one of the region wigets values changed, update'
        row = self.tablewidget_regions.current_row()
        data = self.tablewidget_regions.value
        name = list(data.keys())[row]
        key = value.keys()[0]

        if 'to' in key or 'from' in key:
            item = key.split('_')
            index = ['x', 'y', 'z'].index(item[1])
            data[name][item[0]][index] = value.values()[0]

            self.vtkwidget.update_region(name, data[name])

        elif 'name' in key and name != value.values()[0]:

            new_name = get_unique_string(value.values()[0], list(data.keys()))
            data = OrderedDict([(new_name, v) if k == name else (k, v) for
                                k, v in data.items()])

            self.vtkwidget.change_region_name(name, new_name)

        elif 'type' in key:
            data[name]['type'] = value.values()[0]

            self.vtkwidget.change_region_type(name, data[name])
            self.enable_disable_widgets(name)

        elif 'stl_shape' in key:
            data[name]['stl_shape'] = value.values()[0]

        elif 'slice' in key:
            data[name]['slice'] = value.values()[0]

        elif 'filter' in key:
            item = key.split('_')
            index = ['x', 'y', 'z'].index(item[1])
            data[name][item[0]][index] = value.values()[0]

        elif 'deviation_angle' in key:
            data[name]['deviation_angle'] = value.values()[0]

        self.tablewidget_regions.set_value(data)

    def enable_disable_widgets(self, name):

        data = self.tablewidget_regions.value
        # enable stl widgets
        self.groupbox_stl.setEnabled(data[name]['type'] == 'STL')

        if data[name]['type'] == 'point':
            self.lineedit_regions_from_x.setEnabled(True)
            self.lineedit_regions_to_x.setEnabled(False)
            self.lineedit_regions_from_y.setEnabled(True)
            self.lineedit_regions_to_y.setEnabled(False)
            self.lineedit_regions_from_z.setEnabled(True)
            self.lineedit_regions_to_z.setEnabled(False)
        elif data[name]['type'] == 'XY-plane':
            self.lineedit_regions_from_x.setEnabled(True)
            self.lineedit_regions_to_x.setEnabled(True)
            self.lineedit_regions_from_y.setEnabled(True)
            self.lineedit_regions_to_y.setEnabled(True)
            self.lineedit_regions_from_z.setEnabled(True)
            self.lineedit_regions_to_z.setEnabled(False)
        elif data[name]['type'] == 'XZ-plane':
            self.lineedit_regions_from_x.setEnabled(True)
            self.lineedit_regions_to_x.setEnabled(True)
            self.lineedit_regions_from_y.setEnabled(True)
            self.lineedit_regions_to_y.setEnabled(False)
            self.lineedit_regions_from_z.setEnabled(True)
            self.lineedit_regions_to_z.setEnabled(True)
        elif data[name]['type'] == 'YZ-plane':
            self.lineedit_regions_from_x.setEnabled(True)
            self.lineedit_regions_to_x.setEnabled(False)
            self.lineedit_regions_from_y.setEnabled(True)
            self.lineedit_regions_to_y.setEnabled(True)
            self.lineedit_regions_from_z.setEnabled(True)
            self.lineedit_regions_to_z.setEnabled(True)
        else:
            self.lineedit_regions_to_x.setEnabled(True)
            self.lineedit_regions_from_x.setEnabled(True)
            self.lineedit_regions_to_x.setEnabled(True)
            self.lineedit_regions_from_y.setEnabled(True)
            self.lineedit_regions_to_y.setEnabled(True)
            self.lineedit_regions_from_z.setEnabled(True)
            self.lineedit_regions_to_z.setEnabled(True)

    def change_color(self):
        color = QtWidgets.QColorDialog.getColor()

        if color.isValid():
            row = self.tablewidget_regions.current_row()
            data = self.tablewidget_regions.value
            name = list(data.keys())[row]

            data[name]['color'].color = color.getRgbF()[:-1]

            self.toolbutton_color.setStyleSheet(
                "QToolButton{{ background: rgb({},{},{});}}".format(
                    *data[name]['color'].color_int))

            self.vtkwidget.change_region_color(name, data[name]['color'])
