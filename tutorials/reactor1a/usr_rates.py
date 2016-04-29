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

    for line in open('species.inc').readlines():
        if '::' in line:
            exec(line.split('::')[1].strip())

    # ! Reaction specific variables:
    # !`````````````````````````````````````````````````````````````````````//
    #       DOUBLE PRECISION c_A    ! Species A concentration mol/cm^3

    # ! AtoR:  A --> 3.0*R         (mol/cm^3.s)
    # !---------------------------------------------------------------------//

    c_A = FLDVAR.ro_g[IJK-1]*FLDVAR.x_g[IJK-1,A-1]/PHYSPROP.mw_g[A-1]

    RATES[AtoR-1] = FLDVAR.ep_g[IJK-1] * CONSTANT.c[1-1] * (c_A)**CONSTANT.c[2-1]
