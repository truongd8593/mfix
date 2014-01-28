
mach_file="LINUX.F"

FORTRAN_CMD=pgf90
LINK_CMD=pgf90

FORTRAN_EXT=f
OBJ_EXT=o
MODULE_CODE=1

compile=" -c -O3 -Mnosave -Mfreeform -Mextend -byteswapio "
link="  "

compile_d=" -c -g -Mbounds -Mchkptr -Mchkfpstk -Mchkstk -Mfreeform -byteswapio "
link_d=" "
 
