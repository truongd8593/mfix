# set home directory
export MFIX_HOME=~/projects/mfix_development

# set case directory
export CASE_DIR=$MFIX_HOME/tests/mms_test_suite/horizontal_channel_superbee

# load modules
module load gnu/4.9.2 openmpi/1.5.5_gnu4.6

# compile MFIX in ./src/
cd $CASE_DIR/src
$MFIX_HOME/model/make_mfix --dmp --opt=O3 --compiler=gcc --exe=mfix.exe -j
echo "******** MFIX compiled successfully."

cd $CASE_DIR

# remove these files if exists:
if [ -e "de_norms_collected.dat" ]
then
  echo "******** Removing old files..."
  rm de_norms_collected.dat
fi

# Run mesh_8 (i.e., 8x8 for 2D, 8x8x8 for 3D)
echo "******** Running mesh_8..."
$CASE_DIR/src/mfix.exe imax=8 jmax=8 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm $CASE_DIR/{HCS.*,de_norms.dat,out.log}

# Run mesh_16 (i.e., 16x16 for 2D, 16x16x16 for 3D)
echo "******** Running mesh_16..."
$CASE_DIR/src/mfix.exe imax=16 jmax=16 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm $CASE_DIR/{HCS.*,de_norms.dat,out.log}

# Run mesh_32 (i.e., 32x32 for 2D, 32x32x32 for 3D)
echo "******** Running mesh_32..."
mpirun -np 8 $CASE_DIR/src/mfix.exe imax=32 jmax=32 nodesi=4 nodesj=2 nodesk=1 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm $CASE_DIR/{HCS.*,de_norms.dat,out.log}

# Run mesh_64 (i.e., 64x64 for 2D, 64x64x64 for 3D)
echo "******** Running mesh_64..."
mpirun -np 16 $CASE_DIR/src/mfix.exe imax=64 jmax=64 nodesi=8 nodesj=2 nodesk=1 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm $CASE_DIR/{HCS.*,de_norms.dat,out.log}

# Run mesh_128 (i.e., 128x128 for 2D, 128x128x128 for 3D)
echo "******** Running mesh_128..."
mpirun -np 16 $CASE_DIR/src/mfix.exe imax=128 jmax=128 nodesi=8 nodesj=2 nodesk=1 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm $CASE_DIR/{HCS.*,de_norms.dat,out.log}

## Evaluate observed orders
cp $CASE_DIR/../ooa_test.f95 $CASE_DIR
echo "******** Calculating observed orders..."
gfortran -o ooa_test ooa_test.f95
./ooa_test
rm $CASE_DIR/{ooa_test,ooa_test.f95,de_norms_collected.dat}
mv $CASE_DIR/de_l2.dat $CASE_DIR/AUTOTEST/POST_de_l2.dat
mv $CASE_DIR/de_linf.dat $CASE_DIR/AUTOTEST/POST_de_linf.dat
mv $CASE_DIR/ooa_l2.dat $CASE_DIR/AUTOTEST/POST_ooa_l2.dat
mv $CASE_DIR/ooa_linf.dat $CASE_DIR/AUTOTEST/POST_ooa_linf.dat

echo "******** Done."
