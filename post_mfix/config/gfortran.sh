
mach_file="LINUX.F"

FORTRAN_CMD=gfortran
LINK_CMD=gfortran
 
FORTRAN_EXT=f
OBJ_EXT=o
MODULE_CODE=1
                                                                                
compile=" -c -O -I. -fconvert='big-endian' -ffree-form -ffree-line-length-0 "     
link="  "
                                                                   
compile_d=" -c -g -I. -fconvert='big-endian' -ffree-form -ffree-line-length-0 "
link_d="  "
 