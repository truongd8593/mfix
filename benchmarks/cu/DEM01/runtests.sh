#!/bin/csh
## Change into the current working directory
#$ -cwd
##
## The name for the job. It will be displayed this way on qstat
#$ -N VTUNE_TEST
##
#$ -r n
##
module load gnu openmpi/1.10.2_gnu6.1
setenv RUN_NAME "DEM01"
rm -f $RUN_NAME* >& /dev/null
rm -f POST_* >& /dev/null

if ( ! $?VTUNE_CMD ) then
   setenv VTUNE_CMD ""
endif

if ( ! $?MFIX ) then
   setenv MFIX "./mfix"
endif

if ( ! $?LEVEL ) then
   setenv LEVEL "1"
endif

setenv PROCS `perl -e "print $LEVEL*$LEVEL*$LEVEL"`
setenv CELLS `perl -e "print 20*$LEVEL"`
setenv LEN `perl -e "print 0.004*$LEVEL"`

if ( "$LEVEL" == "1" ) then
   setenv MPIRUN ""
else
    setenv MPIRUN "mpirun -np $PROCS"
endif

time $MPIRUN $VTUNE_CMD $MFIX \
     XLENGTH=$LEN IMAX=$CELLS NODESI=$LEVEL \
     YLENGTH=$LEN JMAX=$CELLS NODESJ=$LEVEL \
     ZLENGTH=$LEN KMAX=$CELLS NODESK=$LEVEL \
     DESGRIDSEARCH_IMAX=$CELLS \
     DESGRIDSEARCH_JMAX=$CELLS \
     DESGRIDSEARCH_KMAX=$CELLS \
     IC_X_E\(1\)=$LEN IC_Y_N\(1\)=$LEN IC_Z_T\(1\)=$LEN
