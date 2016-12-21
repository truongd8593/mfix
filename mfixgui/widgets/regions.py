# -*- coding: utf-8 -*-
#!/usr/bin/env python

from __future__ import print_function, absolute_import, unicode_literals, division

import os
import glob
import copy
from collections import OrderedDict
from qtpy import QtWidgets, QtGui, QtCore

from qtpy import uic

# local imports
from mfixgui.tools.general import (get_unique_string, widget_iter, CellColor, get_pixmap)

from mfixgui.project import Equation, ExtendedJSON

DEFAULT_REGION_DATA = {
    'filter': [0, 0, 0],
    'to': [0, 0, 0],
    'deviation_angle': 10,
    'slice': True,
    'from': [0, 0, 0],
    'color': CellColor([1, 0, 0]),
    'stl_shape': 'box',
    'type': 'box',
    'visibility': True}

def safe_int(obj, default=None):
    try:
        return int(obj)
    except:
        return default

def clean_region_dict(region_dict):
    """remove qt objects and values that are equal to the default"""
    clean_dict = {}
    for key, value in region_dict.items():
        if key in ['visible']:
            pass
        elif key == 'color':
            clean_dict['color'] = region_dict['color'].color_hex
        elif isinstance(value, list) and any(isinstance(v, Equation) for v in value):
            clean_dict[key] = value
        elif value != DEFAULT_REGION_DATA[key]:
            clean_dict[key] = value
    return clean_dict


class RegionsWidget(QtWidgets.QWidget):

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.parent = parent
        self.parameter_key_map = {}

        uifiles = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                               'uifiles')
        uic.loadUi(os.path.join(uifiles, 'regions.ui'), self)

        self.extent_lineedits = [self.lineedit_regions_from_x,
                                 self.lineedit_regions_to_x,
                                 self.lineedit_regions_from_y,
                                 self.lineedit_regions_to_y,
                                 self.lineedit_regions_from_z,
                                 self.lineedit_regions_to_z]
        for ext in self.extent_lineedits:
            ext.allow_parameters = True
            ext.dtype = float
            ext.help_text = 'Physical coordinates describing the bounds of the region.'

        self.combobox_regions_type.addItems([
            'point', 'XY-plane', 'XZ-plane', 'YZ-plane', 'box', 'STL', ])

        self.combobox_regions_type.help_text = 'The type of region.'
        self.lineedit_regions_name.help_text = 'Name of the region. Used througout the gui to reference the region'
        self.combobox_stl_shape.help_text = 'Shape to be used to select facets.'
        self.checkbox_slice_facets.help_text = 'Slice facets if the facet is on the selection boundary.'
        for wid in [self.lineedit_filter_x, self.lineedit_filter_y, self.lineedit_filter_z]:
            wid.allow_parameters = True
            wid.help_text = 'Vector to filter facets with. If the facet '\
                            'normal is the same as the vector, then the facet'\
                            ' will be added to the region, else discarded.'
        self.lineedit_deviation_angle.allow_parameters = True
        self.lineedit_deviation_angle.help_text = 'Angle to provide a'\
            'tolerence to the filtering of facets based on the facet normal.'

        self.toolbutton_region_add.clicked.connect(self.new_region)
        self.toolbutton_region_delete.clicked.connect(self.delete_region)
        self.toolbutton_region_delete.setEnabled(False) #Need a selection
        self.toolbutton_region_copy.clicked.connect(self.copy_region)
        self.toolbutton_region_copy.setEnabled(False) #Need a selection
        self.toolbutton_color.clicked.connect(self.change_color)

        tablewidget = self.tablewidget_regions
        tablewidget.dtype = OrderedDict
        tablewidget._setModel()
        tablewidget.set_selection_model()
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
        self.error = self.parent.error
        self.warning = self.warn = self.parent.warn


    def reset_regions(self):
        self.tablewidget_regions.value.clear()
        self.parameter_key_map = {}


    def get_visibility_image(self, visible=True):
        return get_pixmap('visibility.png' if visible else 'visibilityofftransparent.png',
                          16, 16)


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


    def new_region(self, name=None, extents=None, rtype=None, defer_update=False):
        """create a new region"""
        # This is used both as a signal callback and an API function,
        # so there's some complexity with default args/
        if name in (None, True, False): # 'clicked' signal arguments
            name =  'R_1' # shorter than 'region', better than 'new'
        # Would be nice to name a box 'box', plane 'plane', etc but
        #  they all start as 'box' and then are changed to new type

        data = self.tablewidget_regions.value
        name = get_unique_string(name, list(data.keys()))

        reg_dat = data[name] = copy.deepcopy(DEFAULT_REGION_DATA)
        if rtype is not None and extents is not None:
            reg_dat['type'] = rtype
            reg_dat['from'] = extents[0]
            reg_dat['to'] = extents[1]
        reg_dat['visible'] = self.get_visibility_image(True)
        reg_dat['color'].rand()

        self.vtkwidget.new_region(name, data[name])
        self.tablewidget_regions.set_value(data)
        if defer_update:
            return

        self.tablewidget_regions.fit_to_contents()
        self.tablewidget_regions.selectRow(len(data)-1) # Select new row
        self.parent.set_unsaved_flag()
        self.parent.update_nav_tree() # Enable/disable ICs/BCs etc


    def delete_region(self):
        'remove the currently selected region'
        rows = self.tablewidget_regions.current_rows()

        if rows:
            data = self.tablewidget_regions.value
            for row in rows:
                name = list(data.keys())[row]
                if self.check_region_in_use(name):
                    self.parent.message(text="Region %s is in use" % name)
                    return
                deleted_region = data.pop(name)
                self.vtkwidget.delete_region(name)
                self.remove_from_parameter_map(name, deleted_region)

            self.tablewidget_regions.set_value(data)
            self.vtkwidget.render()
            self.parent.set_unsaved_flag()
            self.parent.update_nav_tree()

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
                # update parameters before the new region is added to the dictionary
                self.update_parameter_name(None,new_name, new_region)
                data[new_name] = new_region
                data[new_name]['visible'] = self.get_visibility_image(
                    data[new_name]['visibility'])
                self.vtkwidget.new_region(new_name, data[new_name])

            self.tablewidget_regions.set_value(data)
            self.tablewidget_regions.fit_to_contents()
            self.tablewidget_regions.selectRow(len(data)-1)
            self.parent.set_unsaved_flag()
            self.parent.update_nav_tree()


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
            data = copy.deepcopy(DEFAULT_REGION_DATA)
            name = ''

        # color
        self.toolbutton_color.setStyleSheet(
            "QToolButton{{ background: rgb({},{},{});}}".format(
                *data['color'].color_int))

        self.lineedit_regions_name.updateValue(None, name)
        self.combobox_regions_type.updateValue(None, data['type'])

        # Don't change type of in-use region
        #self.combobox_regions_type.setEnabled(not self.check_region_in_use(name))
        # TODO (Should we just disable deletion and renaming too, rather than
        #  intercept and warn?)

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


    def region_value_changed(self, widget, value, args, name=None,
                             update_param=True):
        'one of the region widgets values changed, update'
        rows = self.tablewidget_regions.current_rows()
        data = self.tablewidget_regions.value
        if name is None:
            name = list(data.keys())[rows[-1]]
        elif name not in data:
            self.error("region %s not defined" % name)
            return
        key = list(value.keys())[0]
        row_data = data[name]

        if self.check_region_in_use(name):
            if key == 'type':
                # Should we just disable the widgets for in-use regions?
                self.parent.message(text="Region %s is in use, cannot change type" % name)
                widget.setCurrentText(row_data.get('type'))
                return

            elif key == 'name':
                new_name = value.get('name')
                if new_name:
                    self.parent.change_region_name(name, new_name)
                else:
                    self.parent.error('invalid value %s' % value)
                    return

        self.parent.set_unsaved_flag()

        if 'to' in key or 'from' in key:
            item = key.split('_')
            index = ['x', 'y', 'z'].index(item[1])
            val = list(value.values())[0]

            # check for equation and special parameters max, min
            if isinstance(val, Equation):
                used = val.get_used_parameters()
                if 'min' in used or 'max' in used:
                    val.eq = val.eq.replace('min', item[1]+'min').replace('max', item[1]+'max')
                    widget.updateValue(None, val)
            if update_param:
                self.update_parameter_map(value[key], name, key)

            # update data dict
            data[name][item[0]][index] = val

            # check for and update point and plane extents
            shape = data[name]['type']
            if shape in ('plane', 'point') and item[1] not in shape.lower():
                self.update_parameter_map(val, name, 'to_'+str(item[1]))
                data[name]['to'][index] = val
                self.extent_lineedits[index*2+1].updateValue(None, val)

            # propagate values
            for update in (self.vtkwidget.update_region,
                           self.parent.update_region):
                update(name, data[name])

        elif 'name' in key and name != list(value.values())[0]:
            new_name = get_unique_string(list(value.values())[0], list(data.keys()))
            data = OrderedDict(((new_name, v) if k == name else (k, v) for
                                (k, v) in data.items()))

            self.vtkwidget.change_region_name(name, new_name)
            # TODO FIXME fit table to contents

            # update parameter map
            self.update_parameter_name(name, new_name, data[new_name])

        elif 'type' in key:
            shape = list(value.values())[0]
            data[name]['type'] = shape
            self.vtkwidget.change_region_type(name, data[name])
            self.enable_disable_widgets(name)

            # check for plane, update extents
            if 'plane' in shape:
                shape = shape.lower()
                index = [d not in shape for d in 'xyz'].index(True)
                f_value = data[name]['from'][index]
                self.update_parameter_map(f_value, name, 'to_'+'xyz'[index])
                data[name]['to'][index] = f_value
                self.extent_lineedits[index*2+1].updateValue(None, f_value)
            # check for point, update extents
            elif shape == 'point':
                for index in (0,1,2):
                    f_value = data[name]['from'][index]
                    data[name]['to'][index] = f_value
                    self.update_parameter_map(f_value, name, 'to_'+'xyz'[index])
                    self.extent_lineedits[index*2+1].updateValue(None, f_value)

        elif 'stl_shape' in key:
            data[name]['stl_shape'] = list(value.values())[0]
            self.vtkwidget.update_region(name, data[name])

        elif 'slice' in key:
            data[name]['slice'] = list(value.values())[0]
            self.vtkwidget.update_region(name, data[name])

        elif 'filter' in key:
            item = key.split('_')
            index = ['x', 'y', 'z'].index(item[1])
            data[name][item[0]][index] = list(value.values())[0]
            self.vtkwidget.update_region(name, data[name])

            if update_param:
                self.update_parameter_map(value[key], name, key)

        elif 'deviation_angle' in key:
            data[name]['deviation_angle'] = list(value.values())[0]
            self.vtkwidget.update_region(name, data[name])

            if update_param:
                self.update_parameter_map(value[key], name, key)

        self.tablewidget_regions.set_value(data)

        if key == 'type':
            self.parent.update_nav_tree() # ICs/BCs availability depends on region types


    def table_value_changed(self, name, key, value):
        # When is this called?
        data = self.tablewidget_regions.value
        if key == 'type':
            self.vtkwidget.change_region_type(name, data[name])
        elif key == 'color':
            self.vtkwidget.change_region_color(name, data[name]['color'])
        self.parent.set_unsaved_flag()


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

            data[name]['color'].color = color

            self.toolbutton_color.setStyleSheet(
                "QToolButton{{ background: rgb({},{},{});}}".format(
                    *data[name]['color'].color_int))

            self.parent.set_unsaved_flag()
            self.vtkwidget.change_region_color(name, data[name]['color'])


    def check_region_in_use(self, name):
        return self.parent.check_region_in_use(name)


    def regions_to_str(self):
        """ convert regions data to a string for saving """
        data = {'order': list(self.tablewidget_regions.value.keys()),
                'regions': {}
                }
        for region in data['order']:
            data['regions'][region] = clean_region_dict(self.tablewidget_regions.value[region])
        return ExtendedJSON.dumps(data)


    def regions_from_str(self, string):
        """ load regions data from a saved string """
        loaded_data = ExtendedJSON.loads(string) # Order of dict has been lost

        if 'order' not in loaded_data:
            return

        data = OrderedDict() # Rebuild data dict

        for region in loaded_data['order']:
            # Copy dictionary entry to ordered dict
            region_data = data[region] = copy.deepcopy(DEFAULT_REGION_DATA)
            region_data.update(loaded_data['regions'][region])
            if 'visibility' in region_data:
                # Create pixmap for 'visible' column
                region_data['visible'] = self.get_visibility_image(region_data['visibility'])
            if 'color' in region_data:
                # Convert to a CellColor object
                region_data['color'] = CellColor(region_data['color'])

            #build parameter map
            for key, value in region_data.items():
                if isinstance(value, list):
                    for v, k in zip(value, ['x', 'y', 'z']):
                        self.update_parameter_map(v, region, '_'.join([key,k]))
                else:
                    self.update_parameter_map(value, region, key)

            self.vtkwidget.new_region(region, region_data)

        self.tablewidget_regions.set_value(data)
        self.tablewidget_regions.fit_to_contents()


    def extract_regions(self, proj, proj_dir=None):
        """ extract regions from IC, BC, PS, IS """
        if self.tablewidget_regions.value:
            # We assume regions_dict has been initialized correctly
            # from mfix_gui_comments.
            return

        stl_files = []
        stl_nums = []
        if proj_dir:
            # look for geometry_#####.stl files
            stl_files = glob.glob(os.path.join(proj_dir, 'geometry_*.stl'))
            # extract numbers
            stl_nums = [safe_int(f.split('.')[0].split('_')[-1]) for f in stl_files]

        for condtype, conds in (('ic_', proj.ics), ('bc_', proj.bcs),
                                ('is_', proj.iss), ('ps_', proj.pss)):
            for cond in conds:
                extents = []
                extents_keys = False
                for key in ('x_w', 'x_e', 'y_s', 'y_n', 'z_b', 'z_t'):
                    key = condtype + key
                    if key in cond:
                        extents.append(float(cond[key]))
                        extents_keys = True
                    else:
                        extents.append(0.0)

                # reformat extents
                extents = [extents[::2], extents[1::2]]

                # infer region type from extents
                rtype = self.get_region_type(extents)

                # create a name
                name = condtype.upper() + str(cond.ind) # "BC_1"

                add = False
                # handle CG_* regions
                if ('bc_type' in cond and cond['bc_type'].value.lower().startswith('cg')):
                    rtype = 'STL'
                    add = True
                    # single stl bc, assume region fills domain
                    if cond.ind == proj.get_value('stl_bc_id'):
                        ext = [Equation(s) for s in ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']]
                        extents = [ext[::2], ext[1::2]]
                        for key, value in [('from', ext[::2]), ('to', ext[1::2])]:
                            for v, k in zip(value, ['x', 'y', 'z']):
                                self.update_parameter_map(v, name, '_'.join([key,k]))
                    else:
                        if cond.ind in stl_nums:
                            extents = self.vtkwidget.get_stl_extents(stl_files[stl_nums.index(cond.ind)])
                            # reformat extents
                            extents = [extents[::2], extents[1::2]]



                # if extents are not already used by a region, add it
                elif not self.check_extents_in_regions(extents) and extents_keys:
                    add = True
                elif not extents_keys:
                    self.warn('{} does not have extents defined and is not a cartesian grid, ignoring'.format(name))
                else:
                    self.warn('could not infer region from {}'.format(name))

                if add:
                    self.new_region(name, extents, rtype, defer_update=True)

    def check_extents_in_regions(self, extents):
        """ check to see if the extents are already in a region """
        region_dict = self.tablewidget_regions.value
        for key in region_dict.keys():
            region = region_dict.get(key)
            region_extent = [region['from'], region['to']]
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
        """return region dict, for use by clients"""
        region_dict = self.tablewidget_regions.value
        return copy.deepcopy(region_dict)  # Allow clients to modify dict

    def get_value(self, name, key):
        """given a region name and value key, return the value"""

        data = self.tablewidget_regions.value.get(name)
        # safe to assume key is always in data?
        if data is None:
            val = None
        elif 'to' in key or 'from' in key or 'filter' in key:
            item = key.split('_')
            index = ['x', 'y', 'z'].index(item[1])
            val = data[item[0]][index]
        else:
            val = data[key]
        return val

    def update_parameter_map(self, new_value, name, key):
        """update the mapping of parameters and keywords"""
        name_key = ','.join([name, key])
        # new params
        new_params = []
        if isinstance(new_value, Equation):
            new_params = new_value.get_used_parameters()

        # old params
        old_value = self.get_value(name, key)

        old_params = []
        if isinstance(old_value, Equation):
            old_params = old_value.get_used_parameters()

        add = set(new_params)-set(old_params)
        for param in add:
            self.add_namekey(param, name_key)

        remove = set(old_params)-set(new_params)
        for param in remove:
            self.remove_namekey(param, name_key)

    def add_namekey(self, param, name_key):
        """add the name_key to the parameter list"""
        if param not in self.parameter_key_map:
            self.parameter_key_map[param] = set()
        self.parameter_key_map[param].add(name_key)

    def remove_namekey(self, param, name_key):
        """remove the name_key from the parameter list"""
        self.parameter_key_map[param].remove(name_key)
        if len(self.parameter_key_map[param]) == 0:
            self.parameter_key_map.pop(param)

    def update_parameters(self, params):
        """parameters have changed, update regions"""
        for param in params:
            if param in self.parameter_key_map:
                for var in self.parameter_key_map[param]:
                    name, key = var.split(',')
                    self.region_value_changed(
                        None, {key: self.get_value(name, key)}, None,
                        name=name, update_param=False)

    def remove_from_parameter_map(self, name, del_region):
        """a region was deleted, make sure to remove from parameter map"""
        for key, value in del_region.items():
            name_key = ','.join([name, key])
            if isinstance(value, list):
                for xyz, item in zip(['x', 'y', 'z'], value):
                    self._remove_key('_'.join([name_key, xyz]), item)
            else:
                self._remove_key(key, value)

    def update_parameter_name(self, old_name, new_name, region):
        """a region name changed, update map"""
        if old_name:
            self.remove_from_parameter_map(old_name, region)
        for k, v in region.items():
            if isinstance(v, list):
                for xyz, item in zip(['x', 'y', 'z'], v):
                    self.update_parameter_map(item, new_name, '_'.join([k, xyz]))
            else:
                self.update_parameter_map(v, new_name, k)

    def _remove_key(self, name_key, value):
        if not isinstance(value, Equation): return

        for param in value.get_used_parameters():
            keys = self.parameter_key_map.get(param, None)
            if keys:
                keys.remove(name_key)
                if len(self.parameter_key_map[param]) == 0:
                    self.parameter_key_map.pop(param)
            else:
                self.warn("couldn't remove {} from {}".format(param, keys))
