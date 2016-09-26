#!/bin/bash -le

module load intel
DIRS=~/exa_profiling/DEM0[1-4]/[1-4]/*

for DIR in ${DIRS}; do
    echo  > ${DIR}.report
    echo >> ${DIR}.report
    echo " DATA FOR DIRECTORY: ${DIR}" >> ${DIR}.report
    amplxe-cl -report summary -r ${DIR} >>  ${DIR}.report
    amplxe-cl -report hotspots -r ${DIR} | head -n 10 | cut -c 1-80 >> ${DIR}.report
done
