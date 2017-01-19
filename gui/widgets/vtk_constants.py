"""
Constants for vtk
"""
from collections import OrderedDict
import vtk
# Qt imports
from qtpy import QtGui

DEFAULT_MESH_NAME = 'mesh.vtu'

CELL_TYPE_ENUM = {
    0:  'empty_cell',
    1:  'vertex',
    2:  'poly_vertex',
    3:  'line',
    4:  'poly_line',
    5:  'triangle',
    6:  'triangle_strip',
    7:  'polygon',
    8:  'pixel',
    9:  'quad',
    10: 'tetra',
    11: 'voxel',
    12: 'hexahedron',
    13: 'wedge',
    14: 'pyramid',
    15: 'pentagonal_prism',
    16: 'hexagonal_prism',
    21: 'quadratic_edge',
    22: 'quadratic_triangle',
    23: 'quadratic_quad',
    24: 'quadratic_tetra',
    25: 'quadratic_hexahedron',
    26: 'quadratic_wedge',
    27: 'quadratic_pyramid',
    28: 'biquadratic_quad',
    29: 'triquadratic_hexahedron',
    30: 'quadratic_linear_quad',
    31: 'quadratic_linear_wedge',
    32: 'biquadratic_quadratic_wedge',
    33: 'biquadratic_quadratic_hexahedron',
    34: 'biquadratic_traingle',
    35: 'cubic_line',
    36: 'quadratic_polygon',
    }

# Primitives
PRIMITIVE_DICT = OrderedDict([
    ('sphere',   vtk.vtkSphereSource),
    ('box',      vtk.vtkCubeSource),
    ('cylinder', vtk.vtkCylinderSource),
    ('cone',     vtk.vtkConeSource),
    ])

DEFAULT_PRIMITIVE_PARAMS = {
    'centerx':         0.0,
    'centery':         0.0,
    'centerz':         0.0,
    'rotationx':       0.0,
    'rotationy':       0.0,
    'rotationz':       0.0,
    'radius':          1.0,
    'directionx':      1.0,
    'directiony':      0.0,
    'directionz':      0.0,
    'lengthx':         1.0,
    'lengthy':         1.0,
    'lengthz':         1.0,
    'height':          1.0,
    'resolution':      10,
    'thetaresolution': 10,
    'phiresolution':   10,
    'visible':         True,
    'geo_type':        'primitive',
    'type':            '',
    }

# Parametrics
PARAMETRIC_DICT = OrderedDict([
    ('torus',           vtk.vtkParametricTorus),
    ('boy',             vtk.vtkParametricBoy),
    ('conic_spiral',    vtk.vtkParametricConicSpiral),
    ('cross_cap',       vtk.vtkParametricCrossCap),
    ('dini',            vtk.vtkParametricDini),
    ('ellipsoid',       vtk.vtkParametricEllipsoid),
    ('enneper',         vtk.vtkParametricEnneper),
    ('figure_8_klein',  vtk.vtkParametricFigure8Klein),
    ('klein',           vtk.vtkParametricKlein),
    ('mobius',          vtk.vtkParametricMobius),
    ('random_hills',    vtk.vtkParametricRandomHills),
    ('roman',           vtk.vtkParametricRoman),
    ('super_ellipsoid', vtk.vtkParametricSuperEllipsoid),
    ('super_toroid',    vtk.vtkParametricSuperToroid),
    ])

DEFAULT_PARAMETRIC_PARAMS = {
    'translationx':       0.0,
    'translationy':       0.0,
    'translationz':       0.0,
    'centerx':            0.0,
    'centery':            0.0,
    'centerz':            0.0,
    'rotationx':          0.0,
    'rotationy':          0.0,
    'rotationz':          0.0,
    'radius':             1.0,
    'radiusx':            1.0,
    'radiusy':            1.0,
    'radiusz':            1.0,
    'ringradius':         1.0,
    'crosssectionradius': 0.5,
    'zscale':             0.125,
    'ascale':             0.8,
    'bscale':             0.2,
    'bfunc':              1.0,
    'cfunc':              0.1,
    'nfunc':              2.0,
    'nhills':             30,
    'variancex':          2.5,
    'scalex':             0.3,
    'variancey':          2.5,
    'scaley':             0.3,
    'amplitude':          2.0,
    'scaleamplitude':     0.3,
    'allowrandom':        True,
    'n1':                 1.0,
    'n2':                 1.0,
    'visible':            True,
    'geo_type':           'parametric',
    'type':               ''}

# Filters
FILTER_DICT = OrderedDict([
    ('clean',                 vtk.vtkCleanPolyData),
    ('fill_holes',            vtk.vtkFillHolesFilter),
    ('triangle',              vtk.vtkTriangleFilter),
    ('decimate',              vtk.vtkDecimatePro),
    ('quadric_decimation',    vtk.vtkQuadricDecimation),
    ('quadric_clustering',    vtk.vtkQuadricClustering),
    ('linear_subdivision',    vtk.vtkLinearSubdivisionFilter),
    ('loop_subdivision',      vtk.vtkLoopSubdivisionFilter),
    ('butterfly_subdivision', vtk.vtkButterflySubdivisionFilter),
    ('smooth',                vtk.vtkSmoothPolyDataFilter),
    ('windowed_sinc',         vtk.vtkWindowedSincPolyDataFilter),
    ('reverse_sense',        vtk.vtkReverseSense),
    ])

DEFAULT_FILTER_PARAMS = {
    'linestopoints':        True,
    'polystolines':         True,
    'stripstopolys':        True,
    'maximumholesize':      1.0,
    'processvertices':      True,
    'processlines':         True,
    'targetreduction':      0.2,
    'preservetopology':     True,
    'splitmesh':            True,
    'deletevertices':       False,
    'divisionsx':           10,
    'divisionsy':           10,
    'divisionsz':           10,
    'autoadjustdivisions':  True,
    'visible':              True,
    'relaxation':           0.01,
    'iterations':           20,
    'boundarysmoothing':    True,
    'featureangle':         45.0,
    'featureedgesmoothing': False,
    'edgeangle':            15.0,
    'passband':             0.1,
    'manifoldsmoothing':    False,
    'normalize':            False,
    'reversecells':         False,
    'reversenormals':       True,
    'geo_type':             'filter',
    'type':                 ''}

DEFAULT_BOOLEAN_PARAMS = {
    'children': [],
    'visible':  True,
    'type':     '',
    'geo_type': 'boolean',
    }

DEFAULT_STL_PARAMS = {
    'type':            'stl',
    'filename':        None,
    'centerx':         0.0,
    'centery':         0.0,
    'centerz':         0.0,
    'rotationx':       0.0,
    'rotationy':       0.0,
    'rotationz':       0.0,
    'translationx':    0.0,
    'translationy':    0.0,
    'translationz':    0.0,
    'visible':         True,
    'geo_type':        'stl'}

DEFAULT_PARAMS = {
    'primitive':   DEFAULT_PRIMITIVE_PARAMS,
    'parametric':  DEFAULT_PARAMETRIC_PARAMS,
    'filter':      DEFAULT_FILTER_PARAMS,
    'boolean':     DEFAULT_BOOLEAN_PARAMS,
    'stl':         DEFAULT_STL_PARAMS,
}


DEFAULT_VISUAL_PROPS = {
    'mesh':            {'color': QtGui.QColor(244,  67,  40), 'visible':True, 'opacity':1, 'rep':'wire'},
    'background_mesh': {'color': QtGui.QColor(100, 182, 247), 'visible':True, 'opacity':1, 'rep':'wire'},
    'geometry':        {'color': QtGui.QColor(224, 224, 224), 'visible':True, 'opacity':1, 'rep':'wire'},
    'regions':         {'color': QtGui.QColor(224, 224, 224), 'visible':True, 'opacity':0.5, 'rep':'solid'},
    }

# add edge color
for geo in list(DEFAULT_VISUAL_PROPS.keys()):
    DEFAULT_VISUAL_PROPS[geo]['edge'] = DEFAULT_VISUAL_PROPS[geo]['color'].darker()
