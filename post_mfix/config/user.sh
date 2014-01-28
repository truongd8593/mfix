
# name of the machine dependent file in the model folder to use
mach_file="LINUX.F"

# fortran and linker commands
FORTRAN_CMD=gfortran
LINK_CMD=gfortran

# This should remain as is 
FORTRAN_EXT=f

# when compiling Fortran code, this is the extension
#      of the files that are created
OBJ_EXT=o

# if Fortran module names end in .mod (lowercase) , use 1
# if they end in .MOD , use 0
#    After compiling MFIX, look in the main MFIX folder
#    for a sub-folder that starts with : OBJECT_FILES_
#    Check whether bc.mod or bc.MOD was created.

MODULE_CODE=1

# non-debug compiler and linker flags                                                                                
compile=" -c -O -I. -fconvert='big-endian' -ffree-form -ffree-line-length-0 "     
link="  "
     
# debug compiler and linker flags                                                              
compile_d=" -c -g -I. -fconvert='big-endian' -ffree-form -ffree-line-length-0 "
link_d="  "
