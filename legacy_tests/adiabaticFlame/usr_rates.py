# !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
# !                                                                      !
# !  Module name: USR_RATES                                              !
# !                                                                      !
# !  Purpose:                                                            !
# !                                                                      !
# !  Author: J.Musser                                   Date: 10-Oct-12  !
# !                                                                      !
# !  Comments:                                                           !
# !                                                                      !
# !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

from mfix import constant as CONSTANT
from mfix import fldvar as FLDVAR
from mfix import physprop as PHYSPROP

def usr_rates(IJK, RATES):
    IJK = IJK - 1

    for line in open('species.inc').readlines():
        if '::' in line:
            exec(line.split('::')[1].strip())

    # ! Reaction specific variables:
    # !`````````````````````````````````````````````````````````````````````//
    #       DOUBLE PRECISION c_O2    ! Oxygen concentration mol/cm^3
    #       DOUBLE PRECISION c_CH4   ! Methane concentration mol/cm^3

    # ! CH4_Comb:  CH4 + 2O2 --> CO2 + 2H2O  (mol/cm^3.s)
    # !---------------------------------------------------------------------//
    # ! Note: The CH4 combustion rate is artificial and used for the
    # ! adiabatic flame test case.

    c_O2  = FLDVAR.ro_g[IJK]*FLDVAR.x_g[IJK,O2-1]/PHYSPROP.mw_g[O2-1]
    c_CH4 = FLDVAR.ro_g[IJK]*FLDVAR.x_g[IJK,CH4-1]/PHYSPROP.mw_g[CH4-1]

    RATES[CH4_Comb-1] = CONSTANT.c[1-1] * FLDVAR.ep_g[IJK] * c_O2 * c_CH4
