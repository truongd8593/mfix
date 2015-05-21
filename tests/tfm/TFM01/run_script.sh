#!/bin/bash -exl

# set case directory
export CASE_DIR=`pwd`

# load modules
module load gnu/4.6.4 openmpi/1.5.5_gnu4.6

# compile MFIX in ./src/
echo "******** Compiling MFIX..."
cd $CASE_DIR
../../../model/make_mfix --dmp --opt=O3 --compiler=gcc --exe=mfix.exe -j

cd $CASE_DIR

# Run mesh_32 (i.e., 32x32 for 2D, 32x32x32 for 3D)
echo "******** Running mesh_32..."
mpirun -np 8 $CASE_DIR/mfix.exe imax=32 jmax=32 nodesi=4 nodesj=2 nodesk=1 > out.log
rm $CASE_DIR/{TFM01.*,out.log}
rm $CASE_DIR/de_norms.dat

rm $CASE_DIR/mfix.exe


echo "******** Done."


# uncomment the following to generate plots:
#echo "******** Generating plots..."
#python plot_results.py &
