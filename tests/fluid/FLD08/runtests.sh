#!/bin/bash -lex

RUN_NAME="FLD08"

rm -f POST_*.dat VTU* &> /dev/null

rm -f ${RUN_NAME}* &> /dev/null
time -p ./mfix

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    numdiff -a 0.000001 -r 0.05 ${test_post_file} \
      $(basename ${test_post_file})
done
