#!/usr/bin/env python

# Information about keyword argument types (phase, species, etc) extracted
# from mfix_user_guide

from __future__ import print_function, absolute_import, unicode_literals, division
import io
import os

thisdir = os.path.abspath(os.path.dirname(__file__))

keyword_args = {}
with io.open(os.path.join(thisdir, 'keyword_args.txt'), encoding='utf-8') as f:
    for line in f:
        line = line.lower().strip()
        key, rest = line.split('(')[:2]
        args = [arg.strip() for arg in rest[:-1].split(',')]
        keyword_args[key] = args

arg_types = set()
for (k,v) in keyword_args.items():
    arg_types.update(v)

keys_by_type = {}
for t in arg_types:
    keys_by_type[t] = [k for (k,v) in keyword_args.items() if t in v]
    keys_by_type[t].sort()
