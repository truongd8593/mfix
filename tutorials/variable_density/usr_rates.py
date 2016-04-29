# !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
# !                                                                      !
# !  Module name: USR_RATES                                              !
# !                                                                      !
# !  Purpose: Hook for user defined reaction rates.                      !
# !                                                                      !
# !  Author: J.Musser                                   Date: 10-Oct-12  !
# !                                                                      !
# !  Comments: Write reaction rates in units of moles/sec.cm^3 (cgs) or  !
# !  kmoles/sec.m^3 (SI). Units should match those specified in the data !
# !  file.
# !                                                                      !
# !  Example reaction: Methane combustion                                !
# !                                                                      !
# !  mfix.dat input:                                                     !
# !``````````````````````````````````````````````````````````````````````!
# !    @(RXNS)                                                           !
# !      CH4_Comb { chem_eq = "CH4 + 2.0*O2 --> CO2 + 2.0*H2O" }         !
# !    @(END)                                                            !
# !``````````````````````````````````````````````````````````````````````!
# !                                                                      !
# !  usr_rates.f input:                                                  !
# !``````````````````````````````````````````````````````````````````````!
# !    c_O2  = (RO_g(IJK)*X_g(IJK,O2)/MW_g(O2))                          !
# !    c_CH4 = (RO_g(IJK)*X_g(IJK,CH4)/MW_g(CH4))                        !
# !    RATES(CH4_Comb) = 2.0d5 * EP_g(IJK) * c_O2 * c_CH4                !
# !``````````````````````````````````````````````````````````````````````!
# !  * Species alias and reaction names given in the data file can be    !
# !    used in reference to the reaction index in RATES and a species    !
# !    index in gas/solids phase variables.                              !
# !                                                                      !
# !  * Additional information is provided in section 5.11 of the code    !
# !    Readme.                                                           !
# !                                                                      !
# !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
import math
from mfix import fldvar as FLDVAR
from mfix import physprop as PHYSPROP
from mfix import toleranc as TOLERANC

def usr_rates(IJK, RATES):

    print ("RATES==",len(RATES),RATES)

    for line in open('species.inc').readlines():
        if '::' in line:
            exec(line.split('::')[1].strip())

    IJK = IJK - 1

    # ! Reaction specific variables:
    # !`````````````````````````````````````````````````````````````````````//
    # ! Bounded phase temperatures (K)
    # DOUBLE PRECISION :: TgX   ! gas phase
    # DOUBLE PRECISION :: TsX   ! solids phase 1
    # DOUBLE PRECISION :: TgsX  ! Average of gas/solids 1

    # ! Molar concentration of solids species (mol/cm^3)
    # DOUBLE PRECISION :: c_Biomass
    # DOUBLE PRECISION :: c_Moisture

    # ! Bound the gas and solids phase temperatures.
    TMAX = TOLERANC.tmax
    MW_s = PHYSPROP.mw_s
    X_s = FLDVAR.x_s
    T_g = FLDVAR.t_g
    T_s = FLDVAR.t_s
    ROP_g = FLDVAR.rop_g
    ROP_s = FLDVAR.rop_s

    TgX  = min(TMAX, T_g[IJK])
    TsX  = min(TMAX, T_s[IJK,1-1])

    # ! Compute the gas/solids average bounded temperature.
    TgsX = (TgX + TsX) / 2

    # ! Compute concentration of moisture (gmole/cm^3)
    c_Moisture = ROP_s[IJK,1-1]*X_s[IJK,1-1,Moisture-1]/MW_s[1-1,Moisture-1]

    # ! Compute concentration of carbon (gmole/cm^3)
    c_Biomass = ROP_s[IJK,1-1]*X_s[IJK,1-1,Biomass-1]/MW_s[1-1,Biomass-1]


    # ! Reaction rates:
    # !`````````````````````````````````````````````````````````````````````//

    RATES[Drying-1]    = 5.13e6 * math.exp(-1.058e4/TsX) * c_Moisture

    RATES[Pyrolysis-1] = 1.30e8 * math.exp(-1.688e4/TsX) * c_Biomass

    RATES[Tarring-1]   = 2.00e8 * math.exp(-1.602e4/TsX) * c_Biomass
    RATES[Charring-1]  = 1.08e7 * math.exp(-1.460e4/TsX) * c_Biomass

    # !      RATES(:) = ZERO
