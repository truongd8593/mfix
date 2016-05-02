# !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
# !                                                                      !
# !  Module name: USR3                                                   !
# !  Author: J. Musser                                  Date: 31-MAR-14  !
# !                                                                      !
# !  Purpose: Write out the heat of reactions calculated in RATES0.      !
# !                                                                      !
# !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
def usr3():

    # ! Error index
    IER = 0
    # ! total heat of reaction.
    # DOUBLE PRECISION :: tHORg, tHORs(1:3)

    # ! Initialize the 'sum' variables for the gas and solids phase heates
    # ! of reaction.
    tHORg = 0.0d0
    tHORs = 0.0d0

    # ! Calculate user defined reaction rates.
    mfix.rrates0(IER)

    # ! Open the POST file again to append the heat of reaction data.
    lUnit = open('POST_Thermo.dat','a')

    # ! Write the table header.
    lUnit.write('Reaction','HORg', 'HORs1', 'HORs2', 'HORs3')

    # ! This case has ten fulid cells and ten reactions. Each fluid cell has
    # ! a reaction rate of "ONE" for a single reaction and "ZERO" for all
    # ! other reactions. Note that the first 49 IJK values are for ghost
    # ! cells and therefore 49 is added to the reaction loop counter to
    # ! map between the rates array and the fluid cells.
    for LC in range(1, NO_OF_RXNS+1):
        IJK = LC + 49

        # ! Write out the calculated heat of reaction for reaction LC.
        lUnit.write(Reaction[LC]%Name[1:18],HOR_g[IJK], HOR_s[IJK,1:3])
        tHORg = tHORg + HOR_g(IJK)
        tHORs = tHORs + HOR_s(IJK,1:3)

    lUnit.write('Total', tHORg, tHORs[1:3])
    lUnit.close()
