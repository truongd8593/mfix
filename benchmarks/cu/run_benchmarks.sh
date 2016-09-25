#!/bin/bash -le

DIRNAME=$(date +%Y-%m-%d_%H%M%S)
VTUNE_RESULTS_BASEDIR=~/exa_profiling
if [ -z "${LEVEL}" ]; then
    LEVEL=1
fi

CORES=$( perl -e "use POSIX; print ceil(${LEVEL}**3/16)*16")

QSUB=
if [ -n "${Q}" ]; then
    QSUB="qsub -v MFIX,VTUNE_CMD -q general -pe mpi ${CORES}"
fi

VTUNE_BASE=
if [ -n "${V}" ]; then
    VTUNE_BASE="/nfs/apps/Compilers/Intel/ParallelStudio/2016.3.067/vtune_amplifier_xe/bin64/amplxe-cl"
    VTUNE_BASE="${VTUNE_BASE} -collect hotspots"
    VTUNE_BASE="${VTUNE_BASE} -target-duration-type=long"
fi
export LEVEL MFIX VTUNE_BASE

for benchmark in DEM01 DEM02 DEM03 DEM04; do
    export MFIX=${PWD}/${benchmark}/mfix
    pushd ${benchmark}
    if [ -n "${V}" ]; then
        export VTUNE_CMD="${VTUNE_BASE} -result-dir ${VTUNE_RESULTS_BASEDIR}/${benchmark}/${LEVEL}/${DIRNAME}"
    fi
    ${QSUB} ./runtests.sh
    popd
done
