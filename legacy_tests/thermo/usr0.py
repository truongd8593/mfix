# !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
# !                                                                      C
# !  Module name: USR0                                                   C
# !  Purpose: This routine is called before the time loop starts and is  C
# !           user-definable.  The user may insert code in this routine  C
# !           or call appropriate user defined subroutines.  This        C
# !           can be used for setting constants and checking errors in   C
# !           data.  This routine is not called from an IJK loop, hence  C
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

from mfix import rxn_com as RXN_COM

def usr0():
    # TYPE(REACTION_BLOCK), POINTER :: RxN

    # ! Open the POST data file and write out the recation data table.
    dat = open('POST_Thermo.dat')
    for L in range(1, NO_OF_RXNS+1):
        # RxN => Reaction(L)
        RXN_COM.write_rxn_summary(RxN, SPECIES_ALIAS_g[:], SPECIES_ALIAS_s[:,:], False, 678)
    dat.close()
