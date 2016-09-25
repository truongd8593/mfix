#!/bin/csh
## Change into the current working directory
#$ -cwd
##
## The name for the job. It will be displayed this way on qstat
#$ -N VTUNE_TEST
##
#$ -r n
##
setenv RUN_NAME "DEM03"
rm -f $RUN_NAME* >& /dev/null

if ( ! $?VTUNE_CMD ) then
    setenv VTUNE_CMD ""
endif

if ( ! $?MFIX ) then
   setenv MFIX "./mfix"
endif

if ( ! $?LEVEL ) then
   setenv LEVEL "1"
endif

setenv PROCS `perl -e "print $LEVEL*$LEVEL"`
setenv CELLS `perl -e "print 4*$LEVEL"`
setenv LEN `perl -e "print 0.0008*$LEVEL"`

if ( "$LEVEL" == "1" ) then
    setenv MPIRUN ""
else
    setenv MPIRUN "mpirun -np $PROCS"
endif

/usr/bin/time -p $MPIRUN $VTUNE_CMD $MFIX \
     XLENGTH=$LEN IMAX=$CELLS NODESI=$LEVEL \
     ZLENGTH=$LEN KMAX=$CELLS NODESK=$LEVEL \
     IC_X_E\(1\)=$LEN IC_Z_T\(1\)=$LEN \
     BC_X_E\(1\)=$LEN BC_Z_T\(1\)=$LEN \
     BC_X_E\(2\)=$LEN BC_Z_T\(2\)=$LEN
