#!/bin/bash -exl

# Run case
echo "******** Running simulation..."

# Set up N = 128
#mpirun -np 4 ./mfix imax=4 jmax=128 no_k=.T. NODESI=1 NODESJ=4 NODESK=1 >& out.log
./mfix imax=4 jmax=128 no_k=.T. >& out.log
rm -f TFM05.* out.log
mv solution_x_velocity_profile.dat N128_Uprofile.dat
mv de_norms.dat N128_denorm.dat

# Set up N = 256
#mpirun -np 4 ./mfix imax=4 jmax=256 no_k=.T. NODESI=1 NODESJ=4 NODESK=1 >& out.log
./mfix imax=4 jmax=256 no_k=.T. >& out.log
rm -f TFM05.* out.log
mv solution_x_velocity_profile.dat N256_Uprofile.dat
mv de_norms.dat N256_denorm.dat

# Set up N = 512
#mpirun -np 4 ./mfix imax=4 jmax=512 no_k=.T. NODESI=1 NODESJ=4 NODESK=1 >& out.log
./mfix imax=4 jmax=512 no_k=.T. >& out.log
rm -f TFM05.* out.log
mv solution_x_velocity_profile.dat N512_Uprofile.dat
mv de_norms.dat N512_denorm.dat

# Set up N = 1024
#mpirun -np 4 ./mfix imax=4 jmax=1024 no_k=.T. NODESI=1 NODESJ=4 NODESK=1 >& out.log
./mfix imax=4 jmax=1024 no_k=.T. >& out.log
rm -f TFM05.* 
mv solution_x_velocity_profile.dat N1024_Uprofile.dat
mv de_norms.dat N1024_denorm.dat

# uncomment the following to generate plots:
#echo "******** Generating plots..."
#python plot_results.py &
