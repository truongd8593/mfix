#!/bin/bash -exl

# set case directory
export CASE_DIR=`pwd`

# load modules
module load gnu/4.6.4 openmpi/1.5.5_gnu4.6

# compile MFIX in ./src/
echo "******** Compiling MFIX..."
cd $CASE_DIR
#../../../model/make_mfix --dmp --opt=O3 --compiler=gcc --exe=mfix.exe -j
../../../model/make_mfix --serial --opt=O3 --compiler=gcc --exe=mfix.exe -j


cd $CASE_DIR

# remove old result files
rm -f error_summary.dat
rm -f solution_tec_block.dat

# Run case
echo "******** Running simulation..."
#mpirun -np 1 $CASE_DIR/mfix.exe nodesi=1 nodesj=1 nodesk=1 
$CASE_DIR/mfix.exe Discretize=9*0
rm -f TFM04.*

$CASE_DIR/mfix.exe Discretize=9*3
rm -f TFM04.*

$CASE_DIR/mfix.exe Discretize=9*2
rm -f TFM04.*

$CASE_DIR/mfix.exe Discretize=9*5
rm -f TFM04.*

#$CASE_DIR/mfix.exe Discretize=9*4
#rm -f TFM04.*

$CASE_DIR/mfix.exe Discretize=9*7
rm -f TFM04.*

$CASE_DIR/mfix.exe Discretize=9*6
rm -f TFM04.*

$CASE_DIR/mfix.exe Discretize=9*8
rm -f TFM04.*

$CASE_DIR/mfix.exe Discretize=9*9
rm -f TFM04.*

rm $CASE_DIR/mfix.exe

echo "******** Done."

# uncomment the following to generate plots:
#echo "******** Generating plots..."
#python plot_results.py &
