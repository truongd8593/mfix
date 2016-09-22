#!/bin/bash -l

RUN_NAME="DEM03"
rm -f ${RUN_NAME}* &> /dev/null

MFIX=./mfix
if [ -n "$1" ]; then
   MFIX=$1
fi

LEVEL=1
if [ -n "$2" ]; then
   LEVEL=$2
fi

PROCS=$(expr ${LEVEL} \* ${LEVEL})
CELLS=$(expr 4 \* ${LEVEL})
LEN=$(awk "BEGIN {printf \"%.10f\n\", 0.0008*${LEVEL}}")

if [ "${LEVEL}" -eq 1 ]; then
  time -p ${MFIX}
else
  time -p mpirun -np ${PROCS} ${MFIX} \
    XLENGTH=${LEN} IMAX=${CELLS} NODESI=${LEVEL} \
    ZLENGTH=${LEN} KMAX=${CELLS} NODESK=${LEVEL} \
    IC_X_E\(1\)=${LEN} IC_Z_T\(1\)=${LEN} \
    BC_X_E\(1\)=${LEN} BC_Z_T\(1\)=${LEN} \
    BC_X_E\(2\)=${LEN} BC_Z_T\(2\)=${LEN}
fi

