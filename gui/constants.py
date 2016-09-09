# Constants
from collections import OrderedDict

# Solver types
# must match combobox_solver in model.ui
SINGLE, TFM, DEM, PIC, HYBRID = range(5)
CONSTANT, AIR, UDF = 0, 1, 2  # model types
VARIABLE = MIXTURE = OTHER = AIR  # continuum, etc

BC_TYPES = ['MI', 'PO', 'NSW', 'FSW', 'PSW', 'PI', 'MO']

BC_NAMES = ['Mass Inflow', 'Pressure Outflow', 'No Slip Wall',
            'Free Slip Wall', 'Partial Slip Wall',
            'Pressure Inflow', 'Mass Outflow']

(MASS_INFLOW, PRESSURE_OUTFLOW,
 NO_SLIP_WALL, FREE_SLIP_WALL, PARTIAL_SLIP_WALL,
 PRESSURE_INFLOW, MASS_OUTFLOW) = range(7)

(NO_FLUX, SPECIFIED_TEMPERATURE, SPECIFIED_FLUX, CONVECTIVE_FLUX) = range(4)

DEFAULT_BC_TYPE = NO_SLIP_WALL



IS_NAMES = ['Impermeable', 'X-Axis Impermeable', 'Y-Axis Impermeable', 'Z-Axis Impermeable',
            'Semi-permeable', 'X-Axis semi-permeable','Y-Axis semi-permeable','Z-Axis semi-permeable']

IS_TYPES = ['IMPERMEABLE', 'X_IMPERMEABLE', 'Y_IMPERMEABLE', 'Z_IMPERMEABLE',
            'SEMIPERMEABLE', 'X_SEMIPERMEABLE', 'Y_SEMIPERMEABLE', 'Z_SEMIPERMEABLE']

(IMPERMEABLE, X_IMPERMEABLE, Y_IMPERMEABLE, Z_IMPERMEABLE,
 SEMIPERMEABLE, X_SEMIPERMEABLE, Y_SEMIPERMEABLE, Z_SEMIPERMEABLE) = range(8)

DEFAULT_IS_TYPE = IMPERMEABLE



SPECIAL_PARAMETERS = ['min', 'max']
for c in ['x', 'y', 'z']:
    SPECIAL_PARAMETERS.extend([c+'min', c+'max'])

PARAMETER_DICT = OrderedDict([(key, 0.0) for key in SPECIAL_PARAMETERS])
