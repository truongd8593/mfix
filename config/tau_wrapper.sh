# Load some MFIX specific TAU functions.
. $MFIX_CONFIG/tau_fun.sh

# Query the user to select a TAU Makefile.
if test -z $TAU_MAKEFILE; then
  set_TAUMAKE
  export TAU_MAKEFILE
else
  echo " Using TAU_MAKEFILE from environment"
fi


# Build the TAU instrumentation file.
#build_tau_instrumentation

#-tau:show

#if test -z $PDTROOT; then
#  echo "PDToolkit is unavailable"
#  exit
#fi



# Specify PDT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PDT_OPTS=
#MFIX_PDB="mfix_merged.pdb"
# Verbose mode
#PDT_OPTS=$PDT_OPTS"-v "
# Specify free-form Fortran
#PDT_OPTS=$PDT_OPTS"-R free "
# Specify the PDB file.
#PDT_OPTS=$PDT_OPTS"-o"${MFIX_PDB}


# Specify TAU option
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TAU_OPTS=
# Verbose mode
TAU_OPTS=${TAU_OPTS}" -optVerbose"
# Do not remove intermediate .pdb and .inst.* files.
TAU_OPTS=${TAU_OPTS}" -optKeepFiles"
# Incorporate the PDT options.
#TAU_OPTS=${TAU_OPTS}" -optPdtF95Opts=\"${PDT_OPTS}\""
# Specify selective instrumentation file
#TAU_OPTS=${TAU_OPTS}" -optTauSelectFile="${TAU_SEL_FILE}

FORT_FLAGS="\$(TAU_OPTIONS) "${FORT_FLAGS}
FORT_FLAGS3="\$(TAU_OPTIONS) "${FORT_FLAGS3}
FORTRAN_CMD="\$(TAU_COMPILER) "${FORTRAN_CMD}
LINK_CMD="\$(TAU_COMPILER) "${LINK_CMD}











#echo "  TAU_OPTS= ${TAU_OPTS}" >> ${MFIX_TAU_SETTINGS}
#echo "  FORTRAN_CMD=tau_f90.sh \$(TAU_OPTS)" $FORTRAN_CMD >> ${MFIX_TAU_SETTING$
#echo "  FORTRAN77_CMD=\$(TAU_COMPILER) \$(MYOPTIONS)" $FORTRAN_CMD >> ${MFIX_TA$
#echo "  LINK_CMD=\$(TAU_COMPILER) \$(MYOPTIONS)" $FORTRAN_CMD >> ${MFIX_TAU_SET$
#echo "  LINK_FLAGS=" >> ${MFIX_TAU_SETTINGS}

#export FORTRAN_CMD



#if test -z $PDTROOT; then
#  PDTHOME=
#  echo "Checking for fortran parser"
#  command -v f95parse >${PDTHOME} 2>&1
#  if test ! $? = 0; then
# || { echo >&2 "I require foo but it's not installed.  Aborting."; exit 1; }
#{ (eval echo "$as_me:$LINENO: \"$ac_compiler --version </dev/null >&5\"") >&5

#MFIX_TAU_SETTINGS="mfix_tau.make"${MAKE_EXT}
#echo "ifdef USE_TAU" > ${MFIX_TAU_SETTINGS}
#echo "  include $TAU_MAKEFILE" >> ${MFIX_TAU_SETTINGS}
#echo "  TAU_OPTS= ${TAU_OPTS}" >> ${MFIX_TAU_SETTINGS}
#echo "  FORTRAN_CMD=tau_f90.sh \$(TAU_OPTS)" $FORTRAN_CMD >> ${MFIX_TAU_SETTINGS}
#echo "  FORTRAN77_CMD=\$(TAU_COMPILER) \$(MYOPTIONS)" $FORTRAN_CMD >> ${MFIX_TAU_SETTINGS}
#echo "  LINK_CMD=\$(TAU_COMPILER) \$(MYOPTIONS)" $FORTRAN_CMD >> ${MFIX_TAU_SETTINGS}
#echo "  LINK_FLAGS=" >> ${MFIX_TAU_SETTINGS}
#echo "endif" >> ${MFIX_TAU_SETTINGS}
#echo " " >> ${MFIX_TAU_SETTINGS}
