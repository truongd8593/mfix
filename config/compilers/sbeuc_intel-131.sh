# Initialize some variables.
omp=
mpi=
mkl=

mkl_libs=
mpi_libs=
misc_libs=

echo "SBEUC@NETL :: Intel Fortran Compiler"

MODDIRPREFIX="-module "

# Add some additinal flags to the object directory
if [[ -n $USE_MIC ]]; then
    DPO=${DPO_BASE}/${DPO}_INTEL_MIC_SBEUC/;
else
    DPO=${DPO_BASE}/${DPO}_INTEL_SBEUC/;
fi
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

# Debug flags for Intel Fortran
dbg=
if test "${USE_DEBUG}" = "1"; then dbg="-g"; fi

# Common compile flags.
common="-c -I. -convert big_endian -assume byterecl  -diag-disable remark"
if [[ -n $USE_MIC ]]; then
    common=${common}" -mmic"
else
    common=${common}" -xHost"
fi

if test "${USE_CODECOV}" = "1"; then common=${common}" -prof-gen=srcpos"; fi

case $OPT in
  0)echo " Setting flags for debugging."
    dbg="-traceback -check all -fp-model source -O0"
    if [[ -z $USE_MIC ]]; then dbg="$dbg -fpe:0" ; fi
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR ${dbg} -g"
    FORT_FLAGS3="${common} ${mkl} -O0 -g"
    LINK_FLAGS="${omp} -g";;

  1)echo " Setting flags for low optimization."
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR -O1 ${dbg}"
    FORT_FLAGS3="${common} ${mkl} -O1 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  2)echo " Setting flags for medium optimization."
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR -O2 ${dbg}"
    FORT_FLAGS3="${common} ${mkl} -O2 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  3)echo " Setting flags for high optimization."
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR -O3 -no-prec-div -static ${dbg}"
    FORT_FLAGS3="${common} ${mkl} -O3 -no-prec-div -static ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  *)echo "Unsupported optimization level."
    exit;;
esac

if [[ -n $USE_MIC ]]; then LINK_FLAGS=${LINK_FLAGS}" -mmic"; fi
