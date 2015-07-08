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

# Run case
echo "******** Running simulation..."
#mpirun -np 1 $CASE_DIR/mfix.exe nodesi=1 nodesj=1 nodesk=1 
$CASE_DIR/mfix.exe Discretize=9*0
rm -r TFM04.*
mv solution_tec_block.dat solution_tec_block_0.dat
mv error_summary.dat error_summary_0.dat

$CASE_DIR/mfix.exe Discretize=9*3
rm -r TFM04.*
mv solution_tec_block.dat solution_tec_block_3.dat
mv error_summary.dat error_summary_3.dat

$CASE_DIR/mfix.exe Discretize=9*2
rm -r TFM04.*
mv solution_tec_block.dat solution_tec_block_2.dat
mv error_summary.dat error_summary_2.dat

$CASE_DIR/mfix.exe Discretize=9*5
rm -r TFM04.*
mv solution_tec_block.dat solution_tec_block_5.dat
mv error_summary.dat error_summary_5.dat

#$CASE_DIR/mfix.exe Discretize=9*4
#rm -r TFM04.*
#mv solution_tec_block.dat solution_tec_block_4.dat
#mv error_summary.dat error_summary_4.dat

$CASE_DIR/mfix.exe Discretize=9*7
rm -r TFM04.*
mv solution_tec_block.dat solution_tec_block_7.dat
mv error_summary.dat error_summary_7.dat

$CASE_DIR/mfix.exe Discretize=9*6
rm -r TFM04.*
mv solution_tec_block.dat solution_tec_block_6.dat
mv error_summary.dat error_summary_6.dat

$CASE_DIR/mfix.exe Discretize=9*8
rm -r TFM04.*
mv solution_tec_block.dat solution_tec_block_8.dat
mv error_summary.dat error_summary_8.dat

$CASE_DIR/mfix.exe Discretize=9*9
rm -r TFM04.*
mv solution_tec_block.dat solution_tec_block_9.dat
mv error_summary.dat error_summary_9.dat

#rm $CASE_DIR/mfix.exe

echo "******** Done."

# uncomment the following to generate plots:
#echo "******** Generating plots..."
#python plot_results.py &
