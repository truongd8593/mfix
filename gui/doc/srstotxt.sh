#!/usr/bin/env bash

# Retrieve PDF from MS-Word by requesting a printout

pdftotext -nopgbrk 'MFIX-UI SRS.pdf'

# Clean out page headers, footers, non-ascii chars, duplicate blank lines

sed "s/ /#### /;s/ /## /;s/^o /# /;s/‘/'/g;s/’/'/g" < 'MFIX-UI SRS.txt' |
    egrep -v '^Software Requirements Specification|^UI$|^Page||' |
    uniq > spec.txt
