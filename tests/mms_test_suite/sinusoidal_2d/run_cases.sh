# set case directory
export CASE_DIR=`pwd`

# load modules
module load gnu/4.6.4 openmpi/1.5.5_gnu4.6
#module load intel/2013.5.192 intelmpi/4.1.1.036

# copy common files to ./src/
cp ../usr_common/usr_mod.f ./src/usr_mod.f
cp ../usr_common/usr3.f ./src/usr3.f

# compile MFIX in ./src/
echo "******** Compiling MFIX..."
cd $CASE_DIR/src
../../../../model/make_mfix --dmp --opt=O3 --compiler=gcc --exe=mfix.exe -j
#../../../../model/make_mfix --dmp --opt=O0 --compiler=intel --exe=mfix.exe -j

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
rm $CASE_DIR/{MMS2D.*,de_norms.dat,out.log}

# Run mesh_16 (i.e., 16x16 for 2D, 16x16x16 for 3D)
echo "******** Running mesh_16..."
$CASE_DIR/src/mfix.exe imax=16 jmax=16 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm $CASE_DIR/{MMS2D.*,de_norms.dat,out.log}

# Run mesh_32 (i.e., 32x32 for 2D, 32x32x32 for 3D)
echo "******** Running mesh_32..."
mpirun -np 2 $CASE_DIR/src/mfix.exe imax=32 jmax=32 nodesi=2 nodesj=1 nodesk=1 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm $CASE_DIR/{MMS2D.*,de_norms.dat,out.log}

# Run mesh_64 (i.e., 64x64 for 2D, 64x64x64 for 3D)
echo "******** Running mesh_64..."
mpirun -np 16 $CASE_DIR/src/mfix.exe imax=64 jmax=64 nodesi=8 nodesj=2 nodesk=1 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm $CASE_DIR/{MMS2D.*,de_norms.dat,out.log}

# Run mesh_128 (i.e., 128x128 for 2D, 128x128x128 for 3D)
echo "******** Running mesh_128..."
mpirun -np 16 $CASE_DIR/src/mfix.exe imax=128 jmax=128 nodesi=8 nodesj=2 nodesk=1 > out.log
cat $CASE_DIR/de_norms.dat >> $CASE_DIR/de_norms_collected.dat
rm $CASE_DIR/{MMS2D.*,de_norms.dat,out.log}

# Evaluate observed orders
cp ../usr_common/ooa_test.f95 $CASE_DIR
echo "******** Calculating observed orders..."
gfortran -o ooa_test ooa_test.f95
./ooa_test
rm $CASE_DIR/{ooa_test,ooa_test.f95,de_norms_collected.dat}
mv $CASE_DIR/de_l2.dat $CASE_DIR/AUTOTEST/de_l2.dat
mv $CASE_DIR/de_linf.dat $CASE_DIR/AUTOTEST/de_linf.dat
mv $CASE_DIR/ooa_l2.dat $CASE_DIR/AUTOTEST/ooa_l2.dat
mv $CASE_DIR/ooa_linf.dat $CASE_DIR/AUTOTEST/ooa_linf.dat

rm $CASE_DIR/src/usr_mod.f
rm $CASE_DIR/src/usr3.f
rm $CASE_DIR/src/mfix.exe

echo "******** Done."
