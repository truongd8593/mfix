# Constants
from collections import OrderedDict

# Solver types
# must match combobox_solver in model.ui
SINGLE, TFM, DEM, PIC, HYBRID = range(5)
CONSTANT, AIR, UDF = 0, 1, 2  # model types
VARIABLE = MIXTURE = OTHER = AIR  # continuum, etc

SPECIAL_PARAMETERS = ['min', 'max']
for c in ['x', 'y', 'z']:
    SPECIAL_PARAMETERS.extend([c+'min', c+'max'])

PARAMETER_DICT = OrderedDict([(key, 0.0) for key in SPECIAL_PARAMETERS])
