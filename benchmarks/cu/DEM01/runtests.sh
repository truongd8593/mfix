#!/bin/bash -l

RUN_NAME="DEM01"

MFIX=./mfix
if [ -n "$1" ]; then
   MFIX=$1
fi

rm -f POST_* &> /dev/null

for LEVEL in 2; do
  rm -f ${RUN_NAME}* &> /dev/null

  PROCS=$(expr ${LEVEL} \* ${LEVEL} \* ${LEVEL})
  CELLS=$(expr 20 \* ${LEVEL})
  LEN=$(awk "BEGIN {printf \"%.10f\n\", 0.004*${LEVEL}}")

  time -p mpirun -np ${PROCS} ${MFIX} \
    XLENGTH=${LEN} IMAX=${CELLS} NODESI=${LEVEL} \
    YLENGTH=${LEN} JMAX=${CELLS} NODESJ=${LEVEL} \
    ZLENGTH=${LEN} KMAX=${CELLS} NODESK=${LEVEL} \
    DESGRIDSEARCH_IMAX=${CELLS} \
    DESGRIDSEARCH_JMAX=${CELLS} \
    DESGRIDSEARCH_KMAX=${CELLS} \
    IC_X_E\(1\)=${LEN} IC_Y_N\(1\)=${LEN} IC_Z_T\(1\)=${LEN}
done

#post_dats=AUTOTEST/POST*.dat
#
#for test_post_file in ${post_dats}; do
#    numdiff -a 0.000001 -r 0.05 ${test_post_file} $(basename ${test_post_file})
#done
