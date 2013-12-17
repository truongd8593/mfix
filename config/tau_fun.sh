#######################################################################
# Function: set_TAUMAKE                                               #
#                                                                     #
# Provide a list of TAU Makefiles if any are listed in $TAUROOT/lib/  #
# This routine is invoked when the environment variable TAU_MAKEFILE  #
# is not set.                                                         #
#                                                                     #
#######################################################################
set_TAUMAKE(){

  mlist=
  for lmake in ${TAUROOT}/lib/Makefile.*; do
    mlist=${mlist}" "${lmake}
  done

  if test -z "${mlist}"; then
    echo "Fatal Error. No Makefiles found!"
    exit
  fi

  echo
  echo "=============================================================="
  echo "Available TAU Makefiles:"
  echo "=============================================================="
  opt=1
  for lmake in $mlist; do
    short=$(basename ${lmake})
    echo "  [${opt}] ${short}"
    opt=$((${opt}+1))
  done
  echo
  echo -n "Select a Makefile: [1] "
  read input
  if test -z $input; then input=1; fi

  opt=1
  for lmake in $mlist; do
    if test "${opt}" -eq "${input}"; then
      TAU_MAKEFILE=${lmake}
      break
    else
      opt=$((${opt}+1))
    fi
  done

  if test -z $TAU_MAKEFILE; then
    echo "Fatal Error. TAU_MAKEFILE failed to set!"
    exit
  else
    export TAU_MAKEFILE
  fi
} # END set_TAUMAKE


#######################################################################
# Notes concerning selecting instrumentation:                         #
#                                                                     #
# > TAU is provided either include or exclude lists, never both.      #
# > Lists specify either functions/subroutines or files               #
# > Wildcards can be used in lists. (* and ?)                         #
#######################################################################
build_tau_instrumentation() {

  _tfile=$MFIX_SRC/select.tau

# Use any existing instrumentation file.
  if test -r ${_tfile}; then
    echo " Found previous TAU instrumentation list."

# Begin selective instrumentation:
  else
    echo " Building TAU instrumentation list."
# List functions to omit from profiling.
    echo "# " > ${_tfile}
    echo "# Functions excluded from profiling." >> ${_tfile}
    echo "BEGIN_EXCLUDE_LIST" >> ${_tfile}
    echo "SET_WALL_BC1" >> ${_tfile}
    echo "PARALLEL_MPI::MPI_CHECK" >> ${_tfile}
    echo "G_0" >> ${_tfile}
    echo "EQUAL" >> ${_tfile}
    echo "DROODP_G" >> ${_tfile}
    echo "DG_0DNU" >> ${_tfile}
    echo "BOUNDFUNIJK::BOUND_FUNIJK" >> ${_tfile}
    echo "G_0CS" >> ${_tfile}
    echo "COMPARE" >> ${_tfile}
    echo "BOUNDFUNIJK3::BOUND_FUNIJK3" >> ${_tfile}
    echo "DEBUG::ASSERT_I" >> ${_tfile}
    echo "EOSG" >> ${_tfile}
    echo "CALC_XSI" >> ${_tfile}
    echo "END_EXCLUDE_LIST" >> ${_tfile}
    echo "# " >> ${_tfile}

# List any files to completely exclude from profiling.
    echo "#" >> ${_tfile}
    echo "BEGIN_FILE_EXCLUDE_LIST" >> ${_tfile}
    echo "END_FILE_EXCLUDE_LIST" >> ${_tfile}
    echo "#" >> ${_tfile}

# If an include list is provided, only the funtions listed are profiled
# This list is intentionally left commented.
    echo "# " >> ${_tfile}
    echo "#BEGIN_INCLUDE_LIST" >> ${_tfile}
    echo "#END_INCLUDE_LIST" >> ${_tfile}
    echo "# " >> ${_tfile}

# List any files who's contents should be omitte from profiling.
    echo "# " >> ${_tfile}
    echo "#BEGIN_FILE_INCLUDE_LIST" >> ${_tfile}
    echo "#END_FILE_INCLUDE_LIST" >> ${_tfile}
    echo "#" >> ${_tfile}
  fi

  echo "TAU instrumentation set."
  echo "--> pause"; read input
} # END build_tau_instrumentation


skip_this(){
# Check if MFIX_PDB exists and user still wants to run the parser
# Inform user if overwriting an existing *.pdb else print nothing
    if test -r $MFIX_SRC/$MFIX_PDB; then
      echo " * Overwriting the previously generated $MFIX_PDB"
    fi

# Check if any TAU instrumented file (*.inst.f) exists in the source directory
# before attempting to parse with PDToolkit parser
#    for i in `ls *.f`; do
#      outfile=`echo $i | sed 's/\./\.inst./'`
#      if test -r $MFIX_SRC/$outfile; then
#        echo " "
#        echo "   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!"
#        echo "    TAU instrumented file found in source directory: $MFIX_SRC/$outfile"
#        echo "   Temporarily moving all instrumented files (*.inst.f) to "
#        echo "   $MFIX_SRC/inst_files_temporary before running PDT parser"
#        echo "   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!"	   
#        echo " "


# Create a temporary location to store the *.inst.f files during parsing
#        if test -d $MFIX_SRC/inst_files_temporary; then
#          touch $MFIX_SRC/inst_files_temporary
#        else
#          mkdir $MFIX_SRC/inst_files_temporary
#        fi

# Move the previously instrumented files (*.inst.f) to temporary location
#        mv $MFIX_SRC/*.inst.f $MFIX_SRC/inst_files_temporary
#        mv $MFIX_SRC/*.inst.F $MFIX_SRC/inst_files_temporary
#        inst_moved=1
#      fi
#      break
#    done


#    if test $miss_dir = 0; then
#      if test $opsys = "Linux_64_pgi_xt_psc" -o $opsys = "Linux_64_pgi_xt_nccs"; then
#        PDTPARSER=$PDTROOTDIR/xt3/bin/f95parse
#      else
#        PDTPARSER=$PDTROOTDIR/bin/f95parse
#      fi


  PDTPARSER=f95parse

  if test $USE_DMP = 1; then
    echo " Running parser for parallel compilation"
    $PDTPARSER -I$MFIX_SRC -I$MFIX_SRC/des \
      $MFIX_SRC/*.f $MFIX_SRC/dmp_modules/*.f $MFIX_SRC/des/*.f \
      $MFIX_SRC/dqmom/*.f $MFIX_SRC/thermochemical/*.f  \
      $MFIX_SRC/chem/*.f \
      $MFIX_SRC/cohesion/*.f -R free \
      -o$MFIX_SRC/$MFIX_PDB

  else

    $PDTPARSER -I$MFIX_SRC -I$MFIX_SRC/des \
      $MFIX_SRC/*.f \
      $MFIX_SRC/dmp_modules/mpi_donothing/*.f \
      $MFIX_SRC/des/*.f \
      $MFIX_SRC/dqmom/*.f $MFIX_SRC/thermochemical/*.f  \
      $MFIX_SRC/chem/*.f \
      $MFIX_SRC/cohesion/*.f -R free \
      -o$MFIX_SRC/$MFIX_PDB
  fi



#instrumented files (*.inst.f) from the temporary location
# and remove the temporary directory created for storing these files during parsing
#----
#if test $inst_moved = 1; then
#  mv $MFIX_SRC/inst_files_temporary/*.inst.f $MFIX_SRC
#  mv $MFIX_SRC/inst_files_temporary/*.inst.F $MFIX_SRC
#  rm -f -r $MFIX_SRC/inst_files_temporary
#fi



#----
# Patch the mfix_?.make makefile with the necessary directives for TAU automatic instrumentation
#----
#  case $MODULE_CODE in
#  0)
#    NEWMakfile=$MFIX_SRC/mfix_u.make
#    #
#    # Check if the existing mfix_?.make already includes TAU directives or not
#    if grep '^ifdef USE_TAU' $NEWMakfile >/dev/null ; then
#      echo " "
#      echo " " 
#      echo " *******"
#      echo " * Makefile $NEWMakfile already includes directives for TAU"
#    else
#      echo " "
#      echo " " 
#      echo " *******"
#      echo " * Updating the makefile $NEWMakfile to include directives for TAU"
#      cp -f -r $NEWMakfile $MFIX_SRC/mfix_u.make_ORG
#      cat $OUTfile > $NEWMakfile
#      cat $MFIX_SRC/mfix_u.make_ORG >> $NEWMakfile
#      rm -f -r $MFIX_SRC/$OUTfile
#    fi
#   ;;
#  1)
#    NEWMakfile=$MFIX_SRC/mfix_l.make
#    if grep '^ifdef USE_TAU' $NEWMakfile >/dev/null ; then 
#      echo " "
#      echo " " 
#      echo " *******"
#      echo " * Makefile $NEWMakfile already includes directives for TAU"
#    else
#      echo " "
#      echo " " 
#      echo " *******"
#      echo " * Updating the makefile $NEWMakfile to include directives for TAU"
#      cp -f -r $NEWMakfile $MFIX_SRC/mfix_l.make_ORG
#      cat $OUTfile > $NEWMakfile
#      cat $MFIX_SRC/mfix_l.make_ORG >> $NEWMakfile
#      rm -f -r $MFIX_SRC/$OUTfile
#    fi
#   ;;
#  esac
#  fi
#  echo " *******"
#  echo " "
}
