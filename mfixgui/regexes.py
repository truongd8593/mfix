"""Regexes used by mfix parser & user interface"""

import re

# Regular Expressions
# Note: parsing with regular expressions is a problem, since they
#  can't handle things like matching parens and quotes.

RE_KEYVALUE = re.compile(r"""
    (\w+)                            # Alphanumeric, key name
    \s*                              # possible whitespace
    (?:\(([\d,: ]+)\))?              # Indices, possibly with : and , non-capturing group
    \s*                              # Possible whitespace
    =                                # Equal sign
    \s*                              # Possible whitespace
    (.*?)                            # Value
    (?=(!|$|\w+(\([\d,: ]+\))?\s*=)) # comment ?  further keywords ? lookahead
    """, re.VERBOSE|re.IGNORECASE)

RE_FLOAT_EXP = re.compile(r"""
    ^                                # Beginning of expr
    [-+]?                            # possible sign
    \d+                              # digits
    (\.\d*)?                         # optional decimal point and more digits
    [ED]                             # E or D
    ([+-]?\d+)?                      # exponent, with 'd' or 'e' (not required)
    $                                # end
    """, re.VERBOSE|re.IGNORECASE)

RE_FLOAT = re.compile(r"""
    ^                            # beginning of expr
    [-+]?                        # possible sign
    (?:                          # non-capturing group
    (\.\d+)                      # no leading digits before decimal point
    |                            # or,
    (\d+\.\d*)                   # leading digits, decimal point,  possibly trailing digits
    )                            # end group
    $                            # end
    """, re.VERBOSE|re.IGNORECASE)

RE_INT = re.compile(r"""
    ^[+-]?\d+$                                 # optional sign, one or more digits
    """, re.VERBOSE|re.IGNORECASE)

RE_SHORTHAND = re.compile(r"""
    [\d]+                                      # unsigned integer
    \*                                         # literal *
    (?:                                        # non-capturing group
      '\w+'                                    # single-quoted word
    | "\w+"                                    # double-quoted word
    | \.\w+\.                                  # .TOKEN.
    | [-+]?\d+\.?\d*(?:[DE][-+]?\d+)?          # number
    )                                          # end group
    """, re.VERBOSE|re.IGNORECASE)

# detect whether int/float keyword needs to be upgraded to Equation
RE_MATH = re.compile(r"""
    pi             |                           # universal constant, or
    ^e$            |                           # a lone 'e', or
    [^\d\.]e[^\d]  |                           # 'e', but not as part of exp. notation, or
    [*/()]         |                           # arithmetic ops, parens, or
    .[-+]                                      # +/-, but not at start of number!
""", re.VERBOSE|re.IGNORECASE)
