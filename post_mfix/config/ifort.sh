
mach_file="LINUX.F"

FORTRAN_CMD=ifort
LINK_CMD=ifort

FORTRAN_EXT=f
OBJ_EXT=o
MODULE_CODE=1

                                                                                 
compile="-c -O3 -axP -I. -w -w95 -i_dynamic -ip -convert big_endian -assume byterecl -FR"
link="-I.  "

compile_d="-c -g -w -w95 -convert big_endian -assume byterecl -I. -FR -O0 -fpe0  -CB "
link_d="-I.  "

