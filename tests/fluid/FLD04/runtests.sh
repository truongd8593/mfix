#!/bin/bash -exl

# remove old result files
rm -f error_summary.dat
rm -f solution_tec_block.dat

# set case directory
RUN_NAME="FLD04"

rm -f POST_* &> /dev/null


for SCHEME in 0 1 2 3 5 6 7 8 9; do
   rm -f ${RUN_NAME}* &> /dev/null
   ./mfix Discretize=9*${SCHEME}
done

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    numdiff -a 0.000001 -r 0.05 ${test_post_file} \
      $(basename ${test_post_file})
done
