#!/bin/bash -l

module load intel
DIRS=~/exa_profiling/DEM0[1-4]/[1-4]

for DIR in ${DIRS}; do
    RESULTS=$( find ${DIR}/* -maxdepth 0 -type d )
    for RESULT in ${RESULTS}; do
        echo  > ${RESULT}.report
        echo >> ${RESULT}.report
        echo " DATA FOR DIRECTORY: ${RESULT}" >> ${RESULT}.report
        amplxe-cl -report summary -r ${RESULT} >>  ${RESULT}.report
        amplxe-cl -report hotspots -r ${RESULT} | head -n 10 | cut -c 1-80 >> ${RESULT}.report
    done
done
