# Functions for MPI
. $MFIX_CONFIG/mpi_fun.sh

echo " "
echo "Gfortran on 64 bit machine"

MODULE_CODE=2

# Add some additinal flags to the object directory
DPO=${DPO_BASE}/${DPO}_GNU/
if test ! -d $DPO; then  mkdir $DPO; fi


# Set OpenMP flags.
if test $USE_SMP = 1; then
  omp="-fopenmp"
else
  omp=""
fi


SET_MPI_INCLUDE

# Set MPI flags.
if test $USE_DMP = 1; then
  FORTRAN_CMD=mpif90
  LINK_CMD=mpif90
  mpi="-I$mpi_include"
else
  FORTRAN_CMD=gfortran
  LINK_CMD=gfortran
  mpi=""
fi


# Set library paths:
ode="${DPO}odepack.a"
blas="${DPO}blas90.a"
dgtsv="${DPO}dgtsv90.a"
mpi_libs=
misc_libs=

LIB_FLAGS="${ode} ${blas} ${dgtsv} ${misc_libs} ${mpi_libs} "


# Setup inline object lists.
inline_objs="${DPO}compare.o ${DPO}eosg.o ${DPO}discretize.o "
inline_files="compare.f eosg.f discretize.f "


# --> Verify compiler is in $PATH <-- #


# Base flags for GNU Fortran
common="-c -I. -fconvert='big-endian' "

case $OPT in
  0)echo " Setting flags for debugging."
    dbg="-fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow "
    FORT_FLAGS="$omp $common -ffree-form -ffree-line-length-0 $dbg $mpi -I$DPO -g "
    FORT_FLAGS3="$omp $common $dbg $mpi -I$DPO -g "
    LINK_FLAGS="$omp -g";;

  1)echo " Setting flags for low optimization."
    FORT_FLAGS="$omp $common -ffree-form -ffree-line-length-0 $mpi -I$DPO -O1 "
    FORT_FLAGS3="$omp $common $mpi -I$DPO -O1 "
    LINK_FLAGS="$omp -g";;

  2)echo " Setting flags for medium optimization."
    FORT_FLAGS="$omp $common -ffree-form -ffree-line-length-0 $mpi -I$DPO -O2 $inline "
    FORT_FLAGS3="$omp $common $mpi -I$DPO -O1 "
    LINK_FLAGS="$omp";;

  3)echo " Setting flags for high optimization."
    FORT_FLAGS="$omp $common -ffree-form -ffree-line-length-0 $mpi -I$DPO -O3 $inline "
    FORT_FLAGS3="$omp $common $mpi -I$DPO -O1 "
    LINK_FLAGS="$omp";;

  *)echo "Unsupported optimization level."
    exit;;
esac
