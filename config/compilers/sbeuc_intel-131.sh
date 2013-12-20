# Initialize some variables.
omp=
mpi=
mkl=

mkl_libs=
mpi_libs=
misc_libs=

echo "SBEUC@NETL :: Intel Fortran Compiler"

# Set the Intel module code.
MODULE_CODE=0

# Add some additinal flags to the object directory
DPO=${DPO_BASE}/${DPO}_INTEL_SBEUC/
if test ! -d $DPO; then  mkdir $DPO; fi

# Set OpenMP flags.
if test $USE_SMP = 1; then omp="-openmp"; fi


# Set compiler commands for DMP or serial.
if test $USE_DMP = 1; then
  FORTRAN_CMD=mpiifort
  LINK_CMD=mpiifort
else
  FORTRAN_CMD=ifort
  LINK_CMD=ifort
fi


# --> Verify compiler is in $PATH <-- #


SET_MPI_INCLUDE
if test $USE_DMP = 1; then
  mpi="-I$MPI_INCLUDE_PATH"
  mpi_libs=
fi

# Set generic library information:
ode="${DPO}odepack.a"
blas="${DPO}blas90.a"
dgtsv="${DPO}dgtsv90.a"


# 12-05-2013 // Intel Math Kernel Library Link Line Advisor
# http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#
# Product: Intel MKL 11.1, OS: Linux, Coprocessor: None
# Compiler: Intel Fortran, Architecture: Intel 64
# Static Linking, Interface Layer:  ILP64, Sequential
# with Fortran 95 interfaces for BLAS95 and LAPACK95
#
# Compile line options.
# These compile line options cause MPI runtime issues.
#mkl="-i8 -I${MKL}/include/intel64/ilp64 -I${MKL}/include"
#
# Link line
##kl_libs="-L ${MKLPATH} ${MKLPATH}/libmkl_blas95_ilp64.a "
#mkl_libs=${mkl_libs}" ${MKLPATH}/libmkl_lapack95_ilp64.a"
#mkl_libs=${mkl_libs}" -Wl,--start-group  ${MKLPATH}/libmkl_intel_ilp64.a "
#mkl_libs=${mkl_libs}" ${MKLPATH}/libmkl_core.a ${MKLPATH}/libmkl_sequential.a"
#mkl_libs=${mkl_libs}" -Wl,--end-group -lpthread -lm"

# The Intel MKL is disabled .
mkl_libs="${blas} ${dgtsv}"


LIB_FLAGS="${ode} ${mkl_libs} ${mpi_libs} ${misc_libs}"

# Setup inline object lists.
inline_objs="${DPO}compare.o ${DPO}eosg.o ${DPO}discretize.o"
inline_files="compare.f eosg.f discretize.f"

# Common compile flags.
common="-c -I. -convert big_endian -assume byterecl"
common=${common}" -diag-disable remark -mcmodel medium -shared-intel"

case $OPT in
  0)echo " Setting flags for debugging."
    dbg="-traceback -check all -fpe:0 -fp-model precise -O0"
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR ${dbg} -g"
    FORT_FLAGS3="${common} ${mkl} -O0 -g"
    LINK_FLAGS="${omp} -g";;

  1)echo " Setting flags for low optimization."
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR -O1 -g"
    FORT_FLAGS3="${common} ${mkl} -O1 -g"
    LINK_FLAGS="${omp} -g ";;

  2)echo " Setting flags for medium optimization."
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR -O2"
    FORT_FLAGS3="${common} ${mkl} -O1"
    LINK_FLAGS="${omp}";;

  3)echo " Setting flags for high optimization."
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR -O3"
    FORT_FLAGS3="${common} ${mkl} -O1"
    LINK_FLAGS="${omp}";;

  *)echo "Unsupported optimization level."
    exit;;
esac
