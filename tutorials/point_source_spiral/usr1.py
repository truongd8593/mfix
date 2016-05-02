# !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
# !                                                                      C
# !  Module name: USR1                                                   C
# !  Purpose: This routine is called from the time loop and is           C
# !           user-definable.  The user may insert code in this routine  C
# !           or call appropriate user defined subroutines.  This        C
# !           can be used for setting or checking errors in quantities   C
# !           that vary with time.  This routine is not called from an   C
# !           IJK loop, hence all indices are undefined.                 C               C
# !                                                                      C
# !  Author:                                            Date: dd-mmm-yy  C
# !  Reviewer:                                          Date: dd-mmm-yy  C
# !                                                                      C
# !  Revision Number:                                                    C
# !  Purpose:                                                            C
# !  Author:                                            Date: dd-mmm-yy  C
# !  Reviewer:                                          Date: dd-mmm-yy  C
# !                                                                      C
# !  Literature/Document References:                                     C
# !                                                                      C
# !  Variables referenced:                                               C
# !  Variables modified:                                                 C
# !                                                                      C
# !  Local variables:                                                    C
# !                                                                      C
# !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

from math import cos, sin, pi
from mfix import physprop as PHYSPROP
from mfix import ps as PS
from mfix import run as RUN

def usr1():

    # ! One revolution over a three second period.
    lRad = (2.0e0*pi*RUN.time)/3.0e0

    # ! Calculate the normalized velocity components.
    lV = (0.12e0/0.15e0)
    lU = (0.09e0/0.15e0)*cos(lRad)
    lW = (0.09e0/0.15e0)*sin(lRad)

    # ! Update the gas phase components.
    PS.ps_v_g[1-1] = lV
    PS.ps_u_g[1-1] = lU
    PS.ps_w_g[1-1] = lW

    # ! Update the solids phase components.
    for M in range(1, PHYSPROP.mmax+1):
        PS.ps_v_s[1-1, M-1] = lV
        PS.ps_u_s[1-1, M-1] = lU
        PS.ps_w_s[1-1, M-1] = lW
