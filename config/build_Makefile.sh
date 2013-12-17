echo "Checking Makefile for updates."

# Change the name to the full path.
MAKEFILE=${MFIX_SRC}/${MAKEFILE}
tmpMFILE=${MFIX_SRC}/tmp.make


# Remove the previous tmp file if it exists.
if test -f ${tmpMFILE}; then rm ${tmpMFILE}; fi


# Include any TAU definitions.
if test ${USE_TAU} = 1; then
  echo "include ${TAU_MAKEFILE}" >> ${tmpMFILE}
  echo "TAU_OPTIONS=${TAU_OPTS}" >> ${tmpMFILE}
fi


# Include any NetCDF definitions.
if test ${USE_NETCDF} = 1; then
  LIB_FLAGS="${LIB_FLAGS} ${NETCDF_LIBS}"
  FORT_FLAGS="${FORT_FLAGS} ${NETCDF_INCS}"
fi


# Include any Trilinos definitions.
#if test ${USE_TRILINOS} = 1; then
#  echo "using Trilinos ?"
#fi


# Include the base definitions:
echo "DPO=${DPO}" >> ${tmpMFILE}
echo "OBJ_EXT=${OBJ_EXT}" >> ${tmpMFILE}
echo "FORTRAN_EXT=${FORTRAN_EXT}" >> ${tmpMFILE}
echo "FORT_FLAGS=${FORT_FLAGS}" >> ${tmpMFILE}
echo "FORT_FLAGS3=${FORT_FLAGS3}" >> ${tmpMFILE}
echo "FORTRAN_CMD=${FORTRAN_CMD}" >> ${tmpMFILE}
echo "LINK_FLAGS=${LINK_FLAGS}" >> ${tmpMFILE}
echo "LINK_CMD=${LINK_CMD}" >> ${tmpMFILE}
echo "LIB_FLAGS=${LIB_FLAGS}" >> ${tmpMFILE}
echo "EXEC_FILE=${EXEC_FILE}" >> ${tmpMFILE}


# Copy of the correct base makefile.
case ${MODULE_CODE} in
  0) cat ${MFIX_SRC}/mfix_l.make >> ${tmpMFILE};;
  1) cat ${MFIX_SRC}/mfix_u.make >> ${tmpMFILE};;
  2) cat ${MFIX_SRC}/mfix_l.make_gfor >> ${tmpMFILE};;
  *) echo "Fatal Error.  Unknown module code: ${MODULE_CODE}"
     exit;;
esac


# Swap the module folder references.
if test ${USE_DMP} = 0; then

  ex ${tmpMFILE} < ${MFIX_SRC}/dmp_modules/ex.commands
  ex ${tmpMFILE} < ${MFIX_SRC}/des/ex.commands

  src=${MFIX_SRC}/dmp_modules
  cmp -s ${src}/compar_mod.f ${src}/mpi_donothing/compar_mod.f
  if test ! $? = 0; then
    /bin/cp -f ${src}/compar_mod.f ${src}/mpi_donothing/
  fi

  cmp -s ${src}/gridmap_mod.f ${src}/mpi_donothing/gridmap_mod.f
  if test ! $? = 0; then
    /bin/cp -f ${src}/gridmap_mod.f ${src}/mpi_donothing/
  fi
fi


# Only copy over the tmp make file if it differs from
# the existing make file.
if test -f ${MAKEFILE}; then
  cmp -s ${tmpMFILE} ${MAKEFILE}
  if test ! $? = 0; then
    /bin/cp -f ${tmpMFILE} ${MAKEFILE}
    /bin/rm ${tmpMFILE}
    echo "Makefile was updated."
  else
    echo "Makefile is up to date."
  fi
else
  /bin/cp -f ${tmpMFILE} ${MAKEFILE}
  /bin/rm ${tmpMFILE}
  echo "Makefile was created."
fi

# Ensure that the object directory has write permission.
if test ! -w ${DPO}; then
  echo "Write permission added to object directory"
  chmod -R u+w ${DPO}
else
  echo "Object directory has write permission"
fi
