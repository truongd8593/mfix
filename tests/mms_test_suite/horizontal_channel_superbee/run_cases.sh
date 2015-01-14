rm de_norms.dat

cd mesh_8
mkmfix -r > tmp.dat 2>&1
./mfix.exe > result.dat
cat de_norms.dat >> ../de_norms.dat
# cleanup
rm -rf HCS.* scr* solution* mfix.exe de_norms.dat result.dat tmp.dat
cd ..

cd mesh_16
mkmfix -r > tmp.dat 2>&1
./mfix.exe > result.dat
cat de_norms.dat >> ../de_norms.dat
# cleanup
rm -rf HCS.* scr* solution* mfix.exe de_norms.dat result.dat tmp.dat
cd ..

cd mesh_32
mkmfix -r > tmp.dat 2>&1
./mfix.exe > result.dat
cat de_norms.dat >> ../de_norms.dat
# cleanup
rm -rf HCS.* scr* solution* mfix.exe de_norms.dat result.dat tmp.dat
cd ..

#cd mesh_64
#mkmfix -r > tmp.dat 2>&1
#./mfix.exe > result.dat
#cat de_norms.dat >> ../de_norms.dat
## cleanup
#rm -rf HCS.* scr* solution* mfix.exe de_norms.dat result.dat tmp.dat
#cd ..
#
#cd mesh_128
#mkmfix -r > tmp.dat 2>&1
#./mfix.exe > result.dat
#cat de_norms.dat >> ../de_norms.dat
## cleanup
#rm -rf HCS.* scr* solution* mfix.exe de_norms.dat result.dat tmp.dat
#cd ..

gfortran -o ooa_test ooa_test.f95
./ooa_test
rm ooa_test
