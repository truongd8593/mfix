#!/bin/bash -lex

RUN_NAME="TFM05"
rm -f POST_VEL.dat &> /dev/null

for DELP_X in -3.0 -2.0 -1.0 0.0 1.0 2.0 3.0; do
  for JMAX in 8 16 32 64; do
    rm -f ${RUN_NAME}* &> /dev/null
    time -p ./mfix JMAX=${JMAX} DELP_X=${DELP_X}
  done
done





#numdiff \
#    -a 0.000001 -r 0.05 \
#    --exclude=1:5 --exclude=2:5 \
#    AUTOTEST/POST_POS.dat POST_POS.dat
#
#numdiff \
#    -a 0.000001 -r 0.05 \
#    --exclude=1:5-6 --exclude=2:5-6 \
#    AUTOTEST/POST_VEL.dat POST_VEL.dat
