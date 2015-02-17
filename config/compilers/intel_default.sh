omp=
incs=

mpi_libs=
misc_libs=

echo "Intel Fortarn Compiler"

MODDIRPREFIX="-module "

# Add some additinal flags to the object directory
DPO=${DPO_BASE}/${DPO}_INTEL/
if test ! -d ${DPO}; then mkdir ${DPO}; fi

# Set OpenMP flags.
if test ${USE_SMP} = 1; then omp="-openmp"; fi

# Set MPI flags.
if test ${USE_DMP} = 1; then
  FORTRAN_CMD=mpiifort
  LINK_CMD=mpiifort
else
  FORTRAN_CMD=ifort
  LINK_CMD=ifort
fi

# Set the MPI include path.  This must occur after
# the Fortran command is defined.
SET_MPI_INCLUDE
if test ${USE_DMP} = 1; then
  incs=${incs}" -I${MPI_INCLUDE_PATH}"
fi


# --> Verify compiler is in $PATH <-- #


# Set library paths:
ode="${DPO}odepack.a"
blas="${DPO}blas90.a"
dgtsv="${DPO}dgtsv90.a"

mkl_libs="${blas} ${dgtsv}"
LIB_FLAGS="${ode} ${mkl_libs} ${misc_libs} ${mpi_libs} "

# Debug flags for Intel Fortran
dbg=
if test ${USE_DEBUG} = 1; then dbg="-g"; fi

# Base flags for GNU Fortran
common="-c -I. -convert big_endian -assume byterecl -cpp"

case $OPT in
  0)echo "Setting flags for debugging."
    dbg="-traceback -check all -fpe:0 -fp-model precise"
    FORT_FLAGS="${omp} ${incs} ${common} -FR ${dbg} -O0 -g "
    FORT_FLAGS3="${common} ${dbg} ${incs} -O0 -g "
    LINK_FLAGS="${omp} -g";;

  1)echo "Setting flags for low optimization."
    FORT_FLAGS="${omp} ${incs} ${common} -FR -O1 ${dbg} "
    FORT_FLAGS3="${common} ${incs} -O1 ${dbg} "
    LINK_FLAGS="${omp} ${dbg}";;

  2)echo "Setting flags for medium optimization."
    FORT_FLAGS="${omp} ${common} ${incs} -FR -O2 ${dbg}"
    FORT_FLAGS3="${common} -O1 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  3)echo "Setting flags for high optimization."
    FORT_FLAGS="${omp} ${incs} ${common} -FR -O3"
    FORT_FLAGS3="${common} -O2"
    LINK_FLAGS="${omp}";;

  *)echo "Unsupported optimization level."
    exit;;
esac
