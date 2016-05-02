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
def usr_rates(IJK, RATES):

    # !-----------------------------------------------
    #       INCLUDE 'species.inc'
    #       INCLUDE 'usrnlst.inc'

    # ! This case is designed so that there is one fluid cell for each
    # ! reaction. The first fluid cell start with and IJK value of 50.
    # ! Therefore, subtract of 49 to have a one-to-one map between the
    # ! rates array and fluid cells.
    # !
    # ! This is overly complicated but done so that the table generated in
    # ! usr3 for this test case can be populated without directly modifying
    # ! any real source code.
    L = IJK - 49
    RATES[L-1] = 1
