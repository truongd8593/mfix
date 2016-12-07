# Constants
from collections import OrderedDict



# Solver types
# must match combobox_solver in model.ui
SINGLE, TFM, DEM, PIC, HYBRID = range(5)
CONSTANT, AIR, UDF = 0, 1, 2  # model types
VARIABLE = MIXTURE = OTHER = AIR  # continuum, etc

DRAG_TYPES = ['SYAM_OBRIEN', # (DEFAULT)
              'BVK',
              'GIDASPOW',
              'GIDASPOW_BLEND',
              'GIDASPOW_PCF',
              'GIDASPOW_BLEND_PCF',
              'HYS',
              'KOCH_HILL',
              'KOCH_HILL_PCF',
              'WEN_YU',
              'WEN_YU_PCF',
              'USER_DRAG']

DEFAULT_DRAG_TYPE = 'SYAM_OBRIEN'

TURBULENCE_MODELS = ['NONE', 'MIXING_LENGTH', 'K_EPSILON']
DEFAULT_TURBULENCE_MODEL = 'NONE'

SUBGRID_TYPES = ['NONE', 'IGCI', 'MILIOLI']
DEFAULT_SUBGRID_TYPE = 'NONE'


KT_TYPES = ['ALGEBRAIC', 'LUN_1984', 'IA_NONEP', 'SIMONIN',
            'AHMADI', 'GD_99', 'GTSH', 'GHD']

DEFAULT_KT_TYPE = 'ALGEBRAIC'

FRICTION_MODELS = ['SCHAEFFER', 'SRIVASTAVA', 'NONE']

RDF_TYPES = ['LEBOWITZ', 'LEBOWITZ', #sic
             'MANSOORI', 'MODIFIED_LEBOWITZ', 'MODIFIED_MANSOORI']

DEFAULT_RDF_TYPE='LEBOWITZ'

BLENDING_FUNCTIONS = ['NONE', 'TANH_BLEND', 'SIGM_BLEND']
DEFAULT_BLENDING_FUNCTION = 'NONE'

BC_TYPES = ['MI', 'PO', 'NSW', 'FSW', 'PSW', 'PI', 'MO', 'CYCLIC'] # 'CYCLIC' is not really a bc_type

BC_NAMES = ['Mass Inflow', 'Pressure Outflow', 'No Slip Wall',
            'Free Slip Wall', 'Partial Slip Wall',
            'Pressure Inflow', 'Mass Outflow',
            'Cyclic Boundary']

(MASS_INFLOW, PRESSURE_OUTFLOW,
 NO_SLIP_WALL, FREE_SLIP_WALL, PARTIAL_SLIP_WALL,
 PRESSURE_INFLOW, MASS_OUTFLOW, CYCLIC_BOUNDARY) = range(8)

(NO_FLUX, SPECIFIED_TEMPERATURE, SPECIFIED_FLUX, CONVECTIVE_FLUX) = range(4)

DEFAULT_BC_TYPE = 'NSW'


IS_NAMES = ['Impermeable', 'X-Axis Impermeable', 'Y-Axis Impermeable', 'Z-Axis Impermeable',
            'Semi-permeable', 'X-Axis semi-permeable','Y-Axis semi-permeable','Z-Axis semi-permeable']

IS_TYPES = ['IMPERMEABLE', 'X_IMPERMEABLE', 'Y_IMPERMEABLE', 'Z_IMPERMEABLE',
            'SEMIPERMEABLE', 'X_SEMIPERMEABLE', 'Y_SEMIPERMEABLE', 'Z_SEMIPERMEABLE']

(IMPERMEABLE, X_IMPERMEABLE, Y_IMPERMEABLE, Z_IMPERMEABLE,
 SEMIPERMEABLE, X_SEMIPERMEABLE, Y_SEMIPERMEABLE, Z_SEMIPERMEABLE) = range(8)

DEFAULT_IS_TYPE = 'IMPERMEABLE'

# ./model/param_mod.f:67:      INTEGER, PARAMETER :: DIM_M = 10 # max # of solids phases
DIM_M = 10
#model/param_mod.f:      INTEGER, PARAMETER :: DIM_EQS = 10
DIM_EQS = 10


SPECIAL_PARAMETERS = ['min', 'max']
for c in ['x', 'y', 'z']:
    SPECIAL_PARAMETERS.extend([c+'min', c+'max'])

PARAMETER_DICT = OrderedDict([(key, 0.0) for key in SPECIAL_PARAMETERS])

MAX_RECENT_PROJECTS = 20
