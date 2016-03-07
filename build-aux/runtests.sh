#!/bin/bash -lex

if [ -f runtests.sh ]; then
    exec ./runtests.sh
fi

post_script=AUTOTEST/post.script.NEW

./mfix${EXEEXT}
if [ -e ${post_script} ]; then
    OMP_NUM_THREADS=1 ./postmfix${EXEEXT} < ${post_script}
fi

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
	  numdiff -a 0.000001 -r 0.05 ${test_post_file} $(basename ${test_post_file})
done
