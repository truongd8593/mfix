#Main script
cd horizontal_channel_superbee
echo "Running: horizontal_channel_superbee"
. run_cases.sh
echo "Completed: horizontal_channel_superbee"
cd ..

cd vertical_channel_superbee
echo "Running: vertical_channel_superbee"
. run_cases.sh
echo "Completed: vertical_channel_superbee"
cd ..

cd sinusoidal_2d_superbee
echo "Running: sinusoidal_2d_superbee"
. run_cases.sh
echo "Completed: sinusoidal_2d_superbee"
cd ..

cd sinusoidal_2d
echo "Running: sinusoidal_2d"
. run_cases.sh
echo "Completed: sinusoidal_2d"
cd ..

cd curlbased_3d
echo "Running: curlbased_3d"
. run_cases.sh
echo "Completed: curlbased_3d"
cd ..

cd curlbased_3d_twophase
echo "Running: curlbased_3d_twophase"
. run_cases.sh
echo "Completed: curlbased_3d_twophase"
cd ..

cd bc_nsw_2d
echo "Running: bc_nsw_2d"
. run_cases.sh
echo "Completed: bc_nsw_2d"
cd ..

cd bc_fsw_2d
echo "Running: bc_fsw_2d"
. run_cases.sh
echo "Completed: bc_fsw_2d"
cd ..

cd bc_po_2d
echo "Running: bc_po_2d"
. run_cases.sh
echo "Completed: bc_po_2d"
cd ..

cd bc_nsw_3d
echo "Running: bc_nsw_3d"
. run_cases.sh
echo "Completed: bc_nsw_3d"
cd ..

cd bc_fsw_3d
echo "Running: bc_fsw_3d"
. run_cases.sh
echo "Completed: bc_fsw_3d"
cd ..

cd bc_po_3d
echo "Running: bc_po_3d"
. run_cases.sh
echo "Completed: bc_po_3d"
cd ..

# collect all results
. collect_results
