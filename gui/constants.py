# Constants
import math
from collections import OrderedDict

# Solver types
# must match combobox_solver in model.ui
SINGLE, TFM, DEM, PIC, HYBRID = range(5)
CONSTANT, AIR, UDF = 0, 1, 2  # model types
VARIABLE = MIXTURE = OTHER = AIR  # continuum, etc

PARAMETER_DICT = OrderedDict([
    ('pi', math.pi),
    ('e', math.e),
])
