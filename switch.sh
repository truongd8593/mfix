#!/bin/bash

find . -name "*.f" | xargs perl -pi -e 's/(DES_POS_NEW)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'
find . -name "*.f" | xargs perl -pi -e 's/(DES_POS_OLD)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'
find . -name "*.f" | xargs perl -pi -e 's/(DES_VEL_NEW)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'
find . -name "*.f" | xargs perl -pi -e 's/(DES_VEL_OLD)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'

find . -name "*.f" | xargs perl -pi -e 's/(OMEGA_NEW)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'
find . -name "*.f" | xargs perl -pi -e 's/(OMEGA_OLD)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'
find . -name "*.f" | xargs perl -pi -e 's/(DES_ACC_OLD)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'
find . -name "*.f" | xargs perl -pi -e 's/(DES_ROT_OLD)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'
find . -name "*.f" | xargs perl -pi -e 's/(PPOS)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'
find . -name "*.f" | xargs perl -pi -e 's/(FC)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'
find . -name "*.f" | xargs perl -pi -e 's/(TOW)\(([^\)]*),([^\)]*)\)/\1\(\3,\2\)/gi'
