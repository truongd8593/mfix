#!/usr/bin/env python

data = """

Phase_Change { chem_eq = "Solid --> Gas" }


Drying { chem_eq = "Liquid --> Vapor" }


Char_Combustion1 { chem_eq = "2FC1 + O2 --> 2CO"}

Char_Combustion2 { chem_eq = "2FC2 + O2 --> 2CO"}

Char_Combustion3 { chem_eq = '2.2FC1 + O2 --> 2CO + 0.2sSoot'}

CH4_Combustion { chem_eq = "CH4 + 2O2 --> CO2 + 2H2O" }

CH4_Catalytic1 { chem_eq = "CH4 + 2O2 + 0.FC1 --> CO2 + 2H2O" }

CH4_Catalytic2 { chem_eq = "CH4 + 2O2 --> CO2 + 2H2O + 0.FC1" }

CH4_Catalytic3 { chem_eq = "CH4 + 2O2 + 0.FC1 --> CO2 + 2H2O + 0.FC1" }

Char_to_Char { chem_eq = 'FC1 --> FC2' }  ! Hot --> Cold

Ash_to_Ash { chem_eq = 'Ash2 --> FlyAsh' } ! Cold --> Hot

Skippy_RxN { chem_eq = 'None' }


Ozone_Decomp { chem_eq = "O3 --> 1.5O2" }


RX1F { chem_eq = "SiH4 --> SiH2 + H2"}
RX1R { chem_eq = "SiH2 + H2 --> SiH4"}

RX2F { chem_eq = "Si2H6 --> SiH4 + SiH2"}
RX2R { chem_eq = "SiH4 + SiH2 --> Si2H6"}

RX3  { chem_eq = "SiH4 --> Si + 2H2"}

RX4  { chem_eq = "SiH2 --> Si + H2"}


CHAR_COMBUSTION { chem_eq = "2FC1 + O2 --> 2CO"}

CHAR_CO2  { chem_eq = "FC1 + CO2 --> 2CO"}           ! Forward
CHAR_CO2r { chem_eq = "2CO + 0.FC1 --> Soot + CO2"}  ! Reverse

CO_Combustion { chem_eq = "CO + 0.5O2 --> CO2" }


Combustion_s1 { chem_eq = "2FC1 + O2 --> 2CO"}
Combustion_s2 { chem_eq = "2FC2 + O2 --> 2CO"}

Char_CO2_s1 { chem_eq = "FC1 + CO2 --> 2CO"}            ! Forward
Char_CO2_s1r { chem_eq = "2CO + 0.FC1 --> Soot + CO2"}  ! Reverse

Char_CO2_s2 { chem_eq = "FC2 + CO2 --> 2CO"}            ! Forward
Char_CO2_s2r { chem_eq = "2CO + 0.FC2 --> Soot + CO2"}  ! Reverse

CO_Combustion { chem_eq = "CO + 0.5O2 --> CO2" }

Char_to_Char { chem_eq = 'FC2 --> FC1' }

Ash_to_Ash { chem_eq = 'Ash2 --> Ash1' }


Drying { chem_eq = 'Moisture --> H2O' }


Pyrolysis {

chem_eq = 'Biomass --> ' &
'0.9639 * CO + 0.8771 * CO2 + 0.3491 * CH4 + ' &
'1.6276 * H2 + 1.4210 * H2O'

DH = 3.585d3     ! (cal/mol-biomass)  150.0 J/g-biomass
fracDH(1) = 1.0  ! assign to the coal phase
}


Charring {

chem_eq = 'Biomass --> 8.3334 * Char'

DH = 3.585d3     ! (cal/mol-biomass)  150.0 J/g-biomass
fracDH(1) = 1.0  ! assign to the coal phase
}


Tarring {

chem_eq = 'Biomass --> Tar'

DH = 3.585d3     ! (cal/mol-biomass)  150.0 J/g-biomass
fracDH(1) = 1.0  ! assign to the coal phase
}


AtoR { chem_eq = "A --> 3R" }

"""

def EMIT():
    print(rxn_id, d, idx)

def tokenize(data):
    for line in data.split('\n'):
        line = line.strip()
        tok = ''
        quote = None
        for c in line:
            if c == quote: # matched quotes
                quote = None
                if tok:
                    yield tok
                    tok = ''
            elif c in ('"', "'"):
                quote = c
            elif c in ('{}()='):
                if tok:
                    yield tok
                    tok = ''
                yield c
            elif c=='!' and not quote:
                if tok:
                    yield tok
                    tok = ''
                break
            elif quote:
                tok += c
            elif c.isspace():
                if tok:
                    yield tok
                    tok = ''
            else:
                tok += c
        if tok:
            yield tok
            tok = ''

def is_alnum(s):
    return all((c.isalnum() or c=='_') for c in s)

class ParseError(Exception):
    pass

def valid_rxn_id(s):
    return len(s) <= 32 and is_alnum(s)

def err(expected, got):
    raise ParseError('expected %s, got %s state=%s' % (expected, got, state))

state = 0
def push(t):
    global state, key, val, index, rxn_id, d, idx
    if state == 0: # start of reaction
        d = {}
        idx = {}
        key = ''
        index = None
        if not valid_rxn_id(t):
            err('reaction id', t)
        rxn_id = t
        state = 1
    elif state == 1: # start of { block
        if t != '{':
            err('{', t)
        state = 2
    elif state == 2: # start of key
        if not is_alnum(t):
            err('key', t)
        key = t
        d[key] = ''
        state = 3
    elif state == 3: # key may have args
        if t == '(':
            state = 4
        elif t == '=':
            state = 7
        else:
            err('( or =', t)
    elif state == 4: # start of arg
        if not t.isdigit():
            err('integer', t)
        index = int(t)
        idx[key] = index
        state = 5
    elif state == 5: # expect closing paren
        if t != ')':
            err(')', t)
        state = 6
    elif state == 6: # expect =
        if t != '=':
            err('=', t)
        state = 7
    elif state == 7: # start of value or contination of value
        if t in '{}()=':
            err('string value', t)
        if key not in d:
            err('key', t)
        d[key] += t
        state = 8
    elif state == 8:
        if t == '&':
            state = 7
        else:
            if t == '}':
                EMIT()
                state = 0
            elif not is_alnum(t):
                err('key', t)
            else:
                key = t
                d[key] = ''
                state = 3


for t in tokenize(data):
    #print(t)
    push(t)




# key := string | string (  index )
# val := string | string & val
# expr :=  key = val
# reaction := reaction ID  { expr [expr]+ }
