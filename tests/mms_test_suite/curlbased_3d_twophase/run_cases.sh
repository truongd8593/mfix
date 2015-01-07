rm de_norms.dat

cd mesh_9
mkmfix -r > tmp.dat 2>&1
mpirun -np 1 mfix.exe > result.dat
cat de_norms.dat >> ../de_norms.dat
# cleanup
rm -rf MMS3D* scr* solution* mfix.exe de_norms.dat result.dat tmp.dat
cd ..

cd mesh_17
mkmfix -r > tmp.dat 2>&1
mpirun -np 1 mfix.exe > result.dat
cat de_norms.dat >> ../de_norms.dat
# cleanup
rm -rf MMS3D* scr* solution* mfix.exe de_norms.dat result.dat tmp.dat
cd ..

cd mesh_33
mkmfix -r > tmp.dat 2>&1
mpirun -np 8 mfix.exe > result.dat
cat de_norms.dat >> ../de_norms.dat
# cleanup
rm -rf MMS3D* scr* solution* mfix.exe de_norms.dat result.dat tmp.dat
cd ..

#cd mesh_65
#mkmfix -r > tmp.dat 2>&1
#mpirun -np 8 mfix.exe > result.dat
#cat de_norms.dat >> ../de_norms.dat
## cleanup
#rm -rf MMS3D* scr* solution* mfix.exe de_norms.dat result.dat tmp.dat
#cd ..
#
#cd mesh_129
#mkmfix -r > tmp.dat 2>&1
#mpirun -np 8 mfix.exe > result.dat
#cat de_norms.dat >> ../de_norms.dat
## cleanup
#rm -rf MMS3D* scr* solution* mfix.exe de_norms.dat result.dat tmp.dat
#cd ..

gfortran -o ooa_test ooa_test.f95
./ooa_test
rm ooa_test
