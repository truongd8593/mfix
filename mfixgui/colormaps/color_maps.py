# -*- coding: utf-8 -*-
'''
This script will read the *.rgb files and produce pngs from them to be
displayed in the gui.
'''
from __future__ import print_function, absolute_import, unicode_literals, division
import glob
import os
from collections import OrderedDict

try:
    import vtk
except ImportError:
    vtk = None

SCRIPT_DIRECTORY = os.path.abspath(os.path.dirname(__file__))
COLOR_MAP_EXT = '.rgb'
COLOR_MAPS = glob.glob(os.path.join(SCRIPT_DIRECTORY, '*' + COLOR_MAP_EXT))

def read_rgb(path):
    """
    Given a file path, extract the rgb values.
    """
    rgb_list = []
    with open(path) as f:
        for line in f:
            rgb = line.split()
            if len(rgb) == 3:
                try:
                    rgb_list.append([float(v) for v in rgb])
                except:
                    pass
    return rgb_list

def get_color_map_dict():
    '''return a dictionary of color maps'''
    color_dict = OrderedDict()
    for f in COLOR_MAPS:
        name = os.path.basename(f).replace(COLOR_MAP_EXT, '').lower()
        color_dict[name] = read_rgb(f)
    return color_dict

def build_vtk_lookup_tables():
    '''build and return lookup tables for vtk'''
    LUTs = OrderedDict()
    for name, colors in get_color_map_dict().items():
        lookup_table = vtk.vtkLookupTable()
        lookup_table.SetNumberOfTableValues(len(colors))

        for i, c in enumerate(colors):
            lookup_table.SetTableValue(i, c[0], c[1], c[2], 1.0)

        LUTs[name] = lookup_table
    return LUTs

if __name__ == "__main__":
    '''Use matplotlib to generate previews of the color maps and save them as
    png files so that they can be displayed in the gui'''
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    def plot_color_gradients(name, cmap, gradient):
        fig, axes = plt.subplots(1, figsize=(1, 0.25))
        fig.subplots_adjust(top=1, bottom=0, left=0, right=1)
        axes.imshow(gradient, aspect='auto', cmap=cmap)
        # Turn off *all* ticks & spines, not just the ones with colormaps.
        axes.set_axis_off()
        fig.savefig(name+'.png', dpi=5)

    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    for name, colors in get_color_map_dict().items():
        cmap = LinearSegmentedColormap.from_list(name, colors)
        plot_color_gradients(name, cmap, gradient)
