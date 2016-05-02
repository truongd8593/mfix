# !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
# !                                                                      C
# !  Module name: USR3                                                   C
# !  Purpose: This routine is called after the time loop ends and is
# !           user-definable.  The user may insert code in this routine
# !           or call appropriate user defined subroutines.
# !           This routine is not called from an IJK loop, hence
# !           all indices are undefined.                                 C
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

from mfix import fldvar as FLDVAR

def usr3():

      dat = open('POST_Aflame.dat')
      dat.write('Adiabatic Flame Temperature = ', FLDVAR.t_g[5])
      dat.write('P_g = ', FLDVAR.p_g[5])
      dat.write('CH4=', FLDVAR.x_g[5, 1], ' O2=', FLDVAR.x_g[5, 2])
      dat.write('CO2=', FLDVAR.x_g[5, 3], ' H2O=', FLDVAR.x_g[5, 4])
      dat.write('N2 =', FLDVAR.x_g[5, 5])
      dat.close()
