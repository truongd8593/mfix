# Initialize some variables.
omp=
mpi=
mkl=

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

LIB_FLAGS="${ode} ${mpi_libs} ${misc_libs}"

if [[ -n $USE_MKL ]]; then
    mkl='-mkl'
else
    LIB_FLAGS="${LIB_FLAGS} ${blas} ${dgtsv}"
fi

# Setup inline object lists.
inline_objs="${DPO}compare.o ${DPO}eosg.o ${DPO}discretize.o"
inline_files="compare.f eosg.f discretize.f"

# Debug flags for Intel Fortran
dbg=
if test "${USE_DEBUG}" = "1"; then dbg="-g"; fi

# Common compile flags.
common="-c -I. -grecord-gcc-switches -convert big_endian -assume byterecl"

# To display the flags mfix.exe is compiled with, run:
# >
# > readelf -p .debug_str mfix.exe | grep switches

if [[ -n $USE_MIC ]]; then
    common=${common}" -mmic"
else
    common=${common}" -xHost"
fi

if test "${USE_CODECOV}" = "1"; then common=${common}" -prof-gen=srcpos"; fi

case $OPT in
  0)echo "Setting flags for debugging."
    dbg="-traceback -check all -fp-model source -O0"
    if [[ -z $USE_MIC ]]; then dbg="$dbg -fpe:0" ; fi
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR ${dbg} -g"
    FORT_FLAGS3="${common} ${mkl} -O0 -g"
    LINK_FLAGS="${omp} ${mkl} -g";;

  1)echo "Setting flags for low optimization."
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR -O1 ${dbg}"
    FORT_FLAGS3="${common} ${mkl} -O1 ${dbg}"
    LINK_FLAGS="${omp} ${mkl} ${dbg}";;

  2)echo "Setting flags for medium optimization."
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR -O2 ${dbg}"
    FORT_FLAGS3="${common} ${mkl} -O2 ${dbg}"
    LINK_FLAGS="${omp} ${mkl} ${dbg}";;

  3)echo "Setting flags for high optimization."
    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR -O3 -no-prec-div -static ${dbg}"
    FORT_FLAGS3="${common} ${mkl} -O3 -no-prec-div -static ${dbg}"
    LINK_FLAGS="${omp} ${mkl} ${dbg}";;

  *)echo "Unsupported optimization level."
    exit;;
esac

if [[ -n $USE_MIC ]]; then LINK_FLAGS=${LINK_FLAGS}" -mmic"; fi
