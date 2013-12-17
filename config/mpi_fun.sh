#######################################################################
# Function: set_MPI_INCLUDE                                           #
#                                                                     #
# Try to determine the location of the mpif.h file.                   #
#                                                                     #
#######################################################################
SET_MPI_INCLUDE(){

  if test ${USE_DMP} = 1; then
    echo "Searching for mpif.h file"
# User supplied file include path.
    if test ! -z ${MPI_INCLUDE_PATH}; then
      echo "User defined" > /dev/null 2>&1;
# User supplied mpi path.
    elif test ! -z ${MPI_PATH}; then
      inc=${MPI_PATH}/include
      if test -f ${inc}; then
        MPI_INCLUDE_PATH=${inc}
      fi

# We -might- be able to get it with showme: (Open MPI)
    elif (eval ${FORTRAN_CMD} -showme:incdirs) > /dev/null 2>&1; then
      MPI_INCLUDE_PATH=$(eval ${FORTRAN_CMD} -showme:incdirs)

# We -might- be able to get it with show (Intel MPI)
    elif (eval ${FORTRAN_CMD} -show) > /dev/null 2>&1; then
      out=$(eval ${FORTRAN_CMD} -show)
      for arg in ${out}; do
        dir=$(echo ${arg} | grep include | sed -e 's/-I\(.*\)\/include/\1/')
        if test ! -z ${dir}; then
          inc=${dir}/include
          if test -f ${inc}/mpif.h; then
            MPI_INCLUDE_PATH=${inc}
          fi
        fi
      done
    fi

# Serial runs uses the fake mpi include file.
  else
    MPI_INCLUDE_PATH="${MFIX_SRC}/dmp_modules/mpi_donothing"
  fi

  if test -z ${MPI_INCLUDE_PATH}; then
    echo "  Unable to locate mpif.h"
    echo -n "  Please provide the location of the mpif.h file: "
    read MPI_INCLUDE_PATH
  fi

  if test -f ${MPI_INCLUDE_PATH}/mpif.h; then
    ln -sf ${MPI_INCLUDE_PATH}/mpif.h ${MFIX_SRC}/mpif.h
  else
    echo "  Fatal Error: Unable to locate the mpif.h header file!"
    echo "  Aborting."
    exit
  fi

} # END set_MPI_INCLUDE
