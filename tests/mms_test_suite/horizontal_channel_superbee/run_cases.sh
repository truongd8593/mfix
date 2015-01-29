## set home directory
export MFIX_HOME=~/projects/mfix_development
export CASE_DIR=$MFIX_HOME/tests/mms_test_suite/horizontal_channel_superbee

## compile MFIX in ./src/
cd $CASE_DIR/src
$MFIX_HOME/model/make_mfix --dmp --opt=O0 --compiler=gcc --exe=mfix.exe -j

## Run mesh_8 (i.e., 8x8 for 2D, 8x8x8 for 3D)
cd $CASE_DIR
$CASE_DIR/src/mfix.exe imax=8 jmax=8 > out.log


## gfortran -o ooa_test ooa_test.f95
## ./ooa_test
## rm ooa_test
