# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

import os
import json
import copy
from collections import OrderedDict
from qtpy import QtWidgets, QtGui, QtCore

# TODO: add pyside? There is an issue to add this to qtpy:
# https://github.com/spyder-ide/qtpy/issues/16
# TODO: cache ui file to a py file, use the ui file if newer, else py file
try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic

# local imports
from tools.general import (get_unique_string, widget_iter, CellColor,
                           get_image_path)


class RegionsWidget(QtWidgets.QWidget):

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

        uifiles = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                               'uifiles')
        uic.loadUi(os.path.join(uifiles, 'regions.ui'), self)

        self.extent_lineedits = [self.lineedit_regions_from_x,
                                 self.lineedit_regions_to_x,
                                 self.lineedit_regions_from_y,
                                 self.lineedit_regions_to_y,
                                 self.lineedit_regions_from_z,
                                 self.lineedit_regions_to_z]

        self.combobox_regions_type.addItems(['point', 'XY-plane',
                                             'XZ-plane', 'YZ-plane', 'box',
                                             'STL',
                                             ])

        self.toolbutton_region_add.pressed.connect(self.new_region)
        self.toolbutton_region_delete.pressed.connect(self.delete_region)
        self.toolbutton_region_copy.pressed.connect(self.copy_region)
        self.toolbutton_color.pressed.connect(self.change_color)

        tablewidget = self.tablewidget_regions
        tablewidget.dtype = OrderedDict
        tablewidget._setModel() # Should be in __init__
        tablewidget.set_value(OrderedDict())
        tablewidget.set_columns(['visible', 'color', 'type', 'from', 'to'])
        tablewidget.show_vertical_header(True)
        tablewidget.auto_update_rows(True)
        tablewidget.set_selection_model('cell', multi=False)
        tablewidget.new_selection.connect(self.update_region_parameters)
        tablewidget.clicked.connect(self.cell_clicked)
        tablewidget.default_value = OrderedDict()
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

    def cell_clicked(self, index):
        if index.column() == 0:
            data = self.tablewidget_regions.value
            name = list(data.keys())[index.row()]

            if data[name]['visibility']:
                image = QtGui.QPixmap(get_image_path(
                    'visibilityofftransparent.png'))
            else:
                image = QtGui.QPixmap(get_image_path('visibility.png'))

            vis = data[name]['visibility'] = not data[name]['visibility']
            self.vtkwidget.change_region_visibility(name, vis)

            image = image.scaled(16, 16, QtCore.Qt.KeepAspectRatio,
                                 QtCore.Qt.SmoothTransformation)
            data[name]['visible'] = image

            self.tablewidget_regions.set_value(data)

    def new_region(self, name='new', extents=[[0, 0, 0], [0, 0, 0]],
                   rtype='box'):
        """create a new region"""

        data = self.tablewidget_regions.value

        name = get_unique_string(name, list(data.keys()))

        image = QtGui.QPixmap(get_image_path('visibility.png'))
        image = image.scaled(16, 16, QtCore.Qt.KeepAspectRatio,
                             QtCore.Qt.SmoothTransformation)

        data[name] = {'type':            rtype,
                      'from':            extents[0],
                      'to':              extents[1],
                      'color':           CellColor([1, 0, 0]),
                      'slice':           True,
                      'filter':          [0, 0, 0],
                      'stl_shape':       'box',
                      'deviation_angle': 10,
                      'visibility':      True,
                      'visible':         image }


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

            for widget, value in zip(self.extent_lineedits[::2],
                                     data[name]['from']
                                     ):
                widget.updateValue(None, value)

            for widget, value in zip(self.extent_lineedits[1::2],
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

        enable_list = [True]*6
        if data[name]['type'] == 'point':
            enable_list[1::2] = [False]*3
        elif data[name]['type'] == 'XY-plane':
            enable_list[5] = False
        elif data[name]['type'] == 'XZ-plane':
            enable_list[3] = False
        elif data[name]['type'] == 'YZ-plane':
            enable_list[1] = False

        for widget, enable in zip(self.extent_lineedits, enable_list):
            widget.setEnabled(enable)

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

    def regions_to_str(self):
        """ convert regions data to a string for saving """
        data = {'order': list(self.tablewidget_regions.value.keys()),
                'regions': {}
                }
        for region in data['order']:
            region = data['regions'][region] = {}
            for key in self.tablewidget_regions.value[region].keys():
                if key == 'visible':
                    pass
                elif key == 'color':
                    region[key] = self.tablewidget_regions.value[region][key].color
                else:
                    region[key] = self.tablewidget_regions.value[region][key]
        return json.dumps(data)

    def regions_from_str(self, string):
        """ load regions data from a saved string """
        loaded_data = json.loads(string)

        if 'order' not in loaded_data:
            return

        data = OrderedDict()
        for region in loaded_data['order']:
            data[region] = {}

            for key in loaded_data['regions'][region].keys():

                if key == 'visibility':
                    # create the image
                    if loaded_data['regions'][region]['visibility']:
                        image = QtGui.QPixmap(get_image_path('visibility.png'))
                    else:
                        image = QtGui.QPixmap(get_image_path(
                            'visibilityofftransparent.png'))
                    image = image.scaled(16, 16, QtCore.Qt.KeepAspectRatio,
                                         QtCore.Qt.SmoothTransformation)

                    data[region]['visible'] = image
                    data[region][key] = loaded_data['regions'][region][key]

                elif key == 'color':
                    data[region][key] = CellColor(
                        loaded_data['regions'][region]['color'])
                else:
                    data[region][key] = loaded_data['regions'][region][key]

            self.vtkwidget.new_region(region, data[region])

        self.tablewidget_regions.set_value(data)
        self.tablewidget_regions.fit_to_contents()

    def clear(self):
        """ delete all regions from the widget and vtk """

        for region in self.tablewidget_regions.value.keys():
            self.vtkwidget.delete_region(region)

        self.tablewidget_regions.set_value(OrderedDict())

    def extract_regions(self, proj):
        """ extract regions from IC, BC, PS """

        for condtype, conds in [('ic', proj.ics), ('bc', proj.bcs),
                                ('is', proj.iss), ('ps', proj.pss)]:
            for cond in conds:
                extents = []
                for key in ['{}_x_w', '{}_x_e', '{}_y_s', '{}_y_n', '{}_z_b',
                            '{}_z_t']:
                    key = key.format(condtype)
                    if key in cond:
                        extents.append(float(cond[key]))
                    else:
                        extents.append(0.0)

                # if extents are not already used by a region, add it
                if not self.check_extents_in_regions(extents):

                    extents = [extents[::2], extents[1::2]]

                    if 'bc_type' in cond and \
                            str(cond['bc_type']).startswith('cg'):
                        rtype = 'stl'
                    else:
                        rtype = self.get_region_type(extents)

                    name = '{}_{}'.format(condtype.upper(), cond.ind)
                    self.new_region(name, extents, rtype)

    def check_extents_in_regions(self, extents):
        """ check to see if the extents are already in a region """
        region_dict = self.tablewidget_regions.value
        for key in region_dict.keys():
            region_extent = [None]*6
            region_extent[::2] = region_dict[key]['from']
            region_extent[1::2] = region_dict[key]['to']

            if extents == region_extent:
                return True
        return False

    def get_region_type(self, extents):
        """ given the extents, guess the region type """

        rtype = 'box'
        if extents[0] == extents[1]:
            rtype = 'point'

        for r, f, t in zip(['YZ-plane', 'XZ-plane', 'XY-plane'],
                           extents[0], extents[1]):
            if f == t:
                rtype = r
                break

        return rtype
