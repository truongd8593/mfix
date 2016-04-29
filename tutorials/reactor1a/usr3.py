# !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
# !                                                                      C
# !  Module name: USR3                                                   C
# !  Purpose: This routine is called after the time loop ends and is
# !           user-definable.  The user may insert code in this routine
# !           or call appropriate user defined subroutines.
# !           This routine is not called from an IJK loop, hence
# !           all indices are undefined.                                 C
# !                                                                      C
# !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

from mfix import compar as COMPAR
from mfix import fldvar as FLDVAR
from mfix import geometry as GEOMETRY

def funijk(i,j,k):
    return COMPAR.ijk_array_of[i,j,k]

def usr3():

    IJK1 = funijk(2, 1, 1) - 1
    IJK2 = funijk(2, GEOMETRY.jmax1, 1) - 1
    conv =  1. -  ( (FLDVAR.x_g[IJK2, 0] * FLDVAR.rop_g[IJK2] * FLDVAR.v_g[IJK2])
                    /(FLDVAR.x_g[IJK1, 0] * FLDVAR.rop_g[IJK1] * FLDVAR.v_g[IJK1])
    )
    ff = open('POST_Conversion.dat', 'w+')
    ff.write('Conversion = %s' % conv)
    ff.close()
