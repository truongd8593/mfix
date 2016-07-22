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
from qtpy import uic

# local imports
from tools.general import (get_unique_string, widget_iter, CellColor,
                           get_image_path)


class RegionsWidget(QtWidgets.QWidget):

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.parent = parent

        uifiles = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                               'uifiles')
        uic.loadUi(os.path.join(uifiles, 'regions.ui'), self)

        def get_visibility_image(visible):
            image = QtGui.QPixmap(get_image_path(
                'visibility.png' if visible
                else 'visibilityofftransparent.png'))
            return image.scaled(16, 16, QtCore.Qt.KeepAspectRatio,
                                QtCore.Qt.SmoothTransformation)

        # Cache pixmaps
        self.visibility_image = {True:get_visibility_image(True),
                                 False:get_visibility_image(False)}

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
        self.toolbutton_region_delete.setEnabled(False) #Need a selection
        self.toolbutton_region_copy.pressed.connect(self.copy_region)
        self.toolbutton_region_copy.setEnabled(False) #Need a selection
        self.toolbutton_color.pressed.connect(self.change_color)

        tablewidget = self.tablewidget_regions
        tablewidget.dtype = OrderedDict
        tablewidget._setModel() # FIXME: Should be in __init__
        tablewidget.set_value(OrderedDict())
        tablewidget.set_columns(['visible', 'color', 'type'])
        tablewidget.show_vertical_header(True)
        tablewidget.auto_update_rows(True)
        tablewidget.new_selection.connect(self.update_region_parameters)
        tablewidget.clicked.connect(self.cell_clicked)
        tablewidget.default_value = OrderedDict()
        tablewidget.value_changed.connect(self.table_value_changed)
        self.inhibit_toggle = True

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


    def get_visibility_image(self, visible=True):
        return self.visibility_image[visible]


    def cell_clicked(self, index):
        if self.inhibit_toggle: # Don't toggle visibility on a row selection event
            self.inhibit_toggle = False
            return
        self.inhibit_toggle = False
        if index.column() == 0:
            data = self.tablewidget_regions.value
            name = list(data.keys())[index.row()]

            vis = data[name]['visibility'] = not data[name]['visibility']
            self.vtkwidget.change_region_visibility(name, vis)

            data[name]['visible'] = self.get_visibility_image(vis)

            self.tablewidget_regions.set_value(data)

    def new_region(self, name='new', extents=[[0, 0, 0], [0, 0, 0]],
                   rtype='box', defer_render=False):
        """create a new region"""

        data = self.tablewidget_regions.value
        name = get_unique_string(name, list(data.keys()))
        image = self.get_visibility_image()
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


        self.vtkwidget.new_region(name, data[name], defer_render=defer_render)
        self.tablewidget_regions.set_value(data)
        if defer_render:
            return

        self.tablewidget_regions.fit_to_contents()
        self.tablewidget_regions.selectRow(len(data)-1) # Select new row

    def delete_region(self):
        'remove the currently selected region'
        rows = self.tablewidget_regions.current_rows()

        if rows:
            data = self.tablewidget_regions.value

            for row in rows:
                name = list(data.keys())[row]
                data.pop(name)
                self.vtkwidget.delete_region(name)

            self.tablewidget_regions.set_value(data)

        nrows = len(data)
        if rows[-1] == nrows: # We deleted the last row,
            #https://mfix.netl.doe.gov/gitlab/develop/mfix/issues/99
            if nrows > 0:
                self.tablewidget_regions.selectRow(nrows-1)

        enabled = (nrows > 0) # Any left?
        self.groupbox_region_parameters.setEnabled(enabled)
        self.toolbutton_region_delete.setEnabled(enabled)
        self.toolbutton_region_copy.setEnabled(enabled)
        self.update_region_parameters()


    def copy_region(self):
        'copy the currently selected region'
        rows = self.tablewidget_regions.current_rows()

        if rows:
            data = self.tablewidget_regions.value
            
            for row in rows:
                name = list(data.keys())[row]
                new_region = copy.deepcopy(data[name])

                new_name = get_unique_string(name, list(data.keys()))
                data[new_name] = new_region
                data[new_name]['visible'] = self.get_visibility_image(
                    data[new_name]['visibility'])                    
                self.vtkwidget.new_region(new_name, data[new_name])

            self.tablewidget_regions.set_value(data)
            self.tablewidget_regions.fit_to_contents()
    
    def update_region_parameters(self):
        'a new region was selected, update region widgets'
        rows = self.tablewidget_regions.current_rows()

        self.inhibit_toggle = True
        enabled = bool(rows)
        self.toolbutton_region_delete.setEnabled(enabled)
        self.toolbutton_region_copy.setEnabled(enabled)
        self.groupbox_region_parameters.setEnabled(enabled)

        if enabled:
            data = self.tablewidget_regions.value
            name = list(data.keys())[rows[-1]]
            data = data[name]
            # enable widgets
            self.enable_disable_widgets(name)
        else:
            data = {'filter': [0, 0, 0],
                    'to': [0, 0, 0],
                    'deviation_angle': 10,
                    'slice': True,
                    'from': [0, 0, 0],
                    'color': CellColor(),
                    'stl_shape': 'box',
                    'type': 'box',
                    'visibility': True}
            name = ''

        # color
        self.toolbutton_color.setStyleSheet(
            "QToolButton{{ background: rgb({},{},{});}}".format(
                *data['color'].color_int))

        self.lineedit_regions_name.updateValue(None, name)
        self.combobox_regions_type.updateValue(None, data['type'])

        for widget, value in zip(self.extent_lineedits[::2], data['from']):
            widget.updateValue(None, value)

        for widget, value in zip(self.extent_lineedits[1::2], data['to']):
            widget.updateValue(None, value)

        # stl
        self.combobox_stl_shape.updateValue(None, data['stl_shape'])

        self.checkbox_slice_facets.updateValue(None, data['slice'])

        self.lineedit_deviation_angle.updateValue(None, data['deviation_angle'])

        for widget, value in zip([self.lineedit_filter_x,
                                  self.lineedit_filter_y,
                                  self.lineedit_filter_z],
                                 data['filter']):
            widget.updateValue(None, value)

    def region_value_changed(self, widget, value, args):
        'one of the region wigets values changed, update'
        rows = self.tablewidget_regions.current_rows()
        data = self.tablewidget_regions.value
        name = list(data.keys())[rows[-1]]
        key = value.keys()[0]

        if 'to' in key or 'from' in key:
            item = key.split('_')
            index = ['x', 'y', 'z'].index(item[1])
            data[name][item[0]][index] = value.values()[0]

            self.vtkwidget.update_region(name, data[name])

        elif 'name' in key and name != value.values()[0]:
            # Don't allow rename of a region in active use (ICs, etc)
            if self.check_region_in_use(name):
                self.parent.message(text="Region %s is in use" % name)
                widget.update_value(name)
                return

            new_name = get_unique_string(value.values()[0], list(data.keys()))
            data = OrderedDict(((new_name, v) if k == name else (k, v) for
                                (k, v) in data.items()))

            self.vtkwidget.change_region_name(name, new_name)

        elif 'type' in key:
            data[name]['type'] = value.values()[0]

            self.vtkwidget.change_region_type(name, data[name])
            self.enable_disable_widgets(name)

        elif 'stl_shape' in key:
            data[name]['stl_shape'] = value.values()[0]
            self.vtkwidget.update_region(name, data[name])

        elif 'slice' in key:
            data[name]['slice'] = value.values()[0]
            self.vtkwidget.update_region(name, data[name])

        elif 'filter' in key:
            item = key.split('_')
            index = ['x', 'y', 'z'].index(item[1])
            data[name][item[0]][index] = value.values()[0]
            self.vtkwidget.update_region(name, data[name])

        elif 'deviation_angle' in key:
            data[name]['deviation_angle'] = value.values()[0]
            self.vtkwidget.update_region(name, data[name])

        self.tablewidget_regions.set_value(data)

    def table_value_changed(self, name, key, value):
        data = self.tablewidget_regions.value

        if key == 'type':
            self.vtkwidget.change_region_type(name, data[name])
        elif key == 'color':
            self.vtkwidget.change_region_color(name, data[name]['color'])

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
            rows = self.tablewidget_regions.current_rows()
            data = self.tablewidget_regions.value
            name = list(data.keys())[rows[-1]]

            data[name]['color'].color = color.getRgbF()[:-1]

            self.toolbutton_color.setStyleSheet(
                "QToolButton{{ background: rgb({},{},{});}}".format(
                    *data[name]['color'].color_int))

            self.vtkwidget.change_region_color(name, data[name]['color'])

    def check_region_in_use(self, name):
        # TODO: Check initial conditions, boundary conditions, etc
        return False

    def regions_to_str(self):
        """ convert regions data to a string for saving """
        data = {'order': list(self.tablewidget_regions.value.keys()),
                'regions': {}
                }
        for region in data['order']:
            data['regions'][region] = {}
            for key in self.tablewidget_regions.value[region].keys():
                if key == 'visible':
                    pass
                elif key == 'color':
                    data['regions'][region][key] = self.tablewidget_regions.value[region][key].color
                else:
                    data['regions'][region][key] = self.tablewidget_regions.value[region][key]
        return json.dumps(data)

    def regions_from_str(self, string):
        """ load regions data from a saved string """
        loaded_data = json.loads(string) # Order of dict has been lost

        if 'order' not in loaded_data:
            return

        data = OrderedDict() # Rebuild data dict

        for region in loaded_data['order']:
            # Copy dictionary entry to ordered dict
            region_data = data[region] = loaded_data['regions'][region]
            if 'visibility' in region_data:
                # Create pixmap for 'visible' column
                region_data['visible'] = self.get_visibility_image(region_data['visibility'])
            if 'color' in region_data:
                # Convert to a CellColor object
                region_data['color'] = CellColor(region_data['color'])

            self.vtkwidget.new_region(region, region_data)

        self.tablewidget_regions.set_value(data)
        self.tablewidget_regions.fit_to_contents()

    def extract_regions(self, proj, defer_render=False):
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

                    if ('bc_type' in cond and
                            str(cond['bc_type']).startswith('cg')):
                        rtype = 'stl'
                    else:
                        rtype = self.get_region_type(extents)

                    name = '{}_{}'.format(condtype.upper(), cond.ind)
                    self.new_region(name, extents, rtype, defer_render=defer_render)

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
        else:
            for r, f, t in zip(['YZ-plane', 'XZ-plane', 'XY-plane'],
                               extents[0], extents[1]):
                if f == t:
                    rtype = r
                    break
        return rtype

    def get_region_dict(self):
        """ return region dict, for use by clients"""
        region_dict = self.tablewidget_regions.value
        return copy.deepcopy(region_dict) # Allow clients to modify dict
