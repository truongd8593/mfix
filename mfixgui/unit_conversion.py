#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division

""" Factors for converting CGS to SI models, etc
    converting CGS to SI models"""

cgs_to_SI = {
#<start generated code>
#<end generated code>
}

if (__name__ == '__main__'):
    import os
    import re
    from tools.general import SCRIPT_DIRECTORY
    from tools.namelistparser import buildKeywordDoc

    doc = buildKeywordDoc(os.path.join(SCRIPT_DIRECTORY, os.pardir, 'model'))
    keys = doc.keys()
    keys.sort()

    pat = re.compile(r'\[(.*) in SI.*\]')


    unit_to_SI = { 'm':    1e-2,  # from cm
                   'kg':   1e-3,  # from g
                   'J':    4.184, # from cal
                   'kmol': 1e-3,  # from mol
                   'N':    1e-5,  # from dyne
                   'Pa':   0.1,   # from barye
                   'K':    1,
                   's':    1
    }

    # Try to get units from keyword doc. Add to predefined factors above

    for key in keys[:]:
        if key in cgs_to_SI:
            continue
        if 'gravity' in key:
            factor = 1e-2 # cm -> m
            cgs_to_SI[key] = factor
            continue
        entry = doc[key]
        if not isinstance(entry, dict):
            keys.remove(key)
            del doc[key]
            continue
        dtype = entry['dtype']
        if dtype !=  'DP':
            keys.remove(key)
            continue
        desc = entry['description']
        dlower = desc.lower()

        if 'fraction' in dlower:
            cgs_to_SI[key] = 1
            continue

        match = pat.search(desc)

        if match:
            SI_unit = match.group(1)
            if SI_unit.startswith('(') and SI_unit.endswith(')'):
                SI_unit = SI_unit[1:-1]
            if '/' in SI_unit:
                num, denom = SI_unit.split('/')
            else:
                num, denom = SI_unit, ''
            if num.startswith('(') and num.endswith(')'):
                num = num[1:-1]
            if denom.startswith('(') and denom.endswith(')'):
                denom = denom[1:-1]
            factor = 1.0
            for term in num.split('.'):
                if '^' in term:
                    u, pow = term.split('^')
                    pow = int(pow)
                else:
                    u = term
                    pow = 1
                factor *= unit_to_SI[u]**pow
            if denom:
                for term in denom.split('.'):
                    if '^' in term:
                        u, pow = term.split('^')
                        pow = int(pow)
                    else:
                        u = term
                        pow = 1
                    factor /= unit_to_SI[u]**pow
            cgs_to_SI[key] = round(factor,10)
    infile = open(__file__, 'r')
    outfile = open(__file__ + '.tmp', 'w')
    skip = False
    for line in infile:
        line = line.decode('utf-8')
        if line.startswith('#<start'):
            outfile.write(line)
            for key in keys:
                factor = cgs_to_SI.get(key)
                if factor is None:
                    outfile.write("#   '%s': UNDEFINED,\n" % key)
                else:
                    outfile.write("    '%s' : %s,\n" % (key, factor))
            skip = True
        elif line.startswith('#<end'):
            outfile.write(line)
            skip = False
            continue
        elif not skip:
            outfile.write(line.encode('utf-8'))

    # We could replace the input file, but let's not go crazy
    print("Generated %s" % (__file__ + '.tmp'))



"""
Simulation Units
Simulations can be setup using the International System of Units (SI) or the centimetergram-second system (CGS). Although the majority of units are consistent with the specified systems, there are exceptions. The following table provides the SI and CGS units employed by MFIX for various quantities.


Quantity                      MFIX SI unit              MFIX CGS unit

length, position              meter(m)                  centimeter (cm)
mass                          kilogram (kg)             gram (g)
time                          second (s)                second (s)
thermal temperature           Kelvin (K)                Kelvin (K)
energy†                       Joule (J)                 calorie (cal)
amount of substance‡          kilomole (kmol)           mole (mol)
force                         Newton (1 N = 1 kg·m·s-2) dyne (1 dyn = 1 g·cm·s-2)
pressure                      Pascal (1 Pa = 1 N·m-2)   barye (1 Ba = 1 dyn·cm-2)
dynamic viscosity             Pa·s                      poise (1 P = 1 g·cm-1·s-1)
kinematic viscosity           m2·s-1                    Stokes (1 St = 1 cm2·s-1)
gas constant                  J·K-1·kmol-1              erg·K-1·mol-1
enthalpy                      J                         cal
specific heat                 J·kg-1·K-1                cal·g-1·K-1
thermal conductivity          J·m-1·K-1·s-1             cal·cm-1·K-1·s-1

† The CGS unit for energy is the ergon (1 erg = 1 dyne·cm). This is reflected in MFIX
through the gas constant. However, all thermochemical properties related to the energy
equations are based in calories for CGS units. Entries in the Burcat database are
always specified in terms of calories regardless of the simulation units. MFIX converts
the entries to joules after reading the database when SI units are used.
‡ The SI unit for the amount of a substance is the mole (mol). These units are needed
when specifying reaction rates:
• amount per time per volume for Eulerian model reactions
• amount per time for Lagrangian model reactions
"""
