#!/bin/bsh -f

show_version() {
echo
echo "=============================================================="
echo ""
echo "       MFIX: Multiphase Flow with Intephase eXchanges"
echo ""
echo "                      mfix.netl.doe.gov"
echo ""
echo "                        Version ${MFIX_VERSION}"
echo ""
echo "=============================================================="
echo
echo "MFIX source: ${MFIX_SRC}"
echo "Operating System: ${opsys}"
echo "Processor type: ${proctyp}"
echo

}


show_usage() {
echo ""
echo " ======================================================================"
echo ""
echo " Usage: mfix/model/make_mfix [OPTIONS]... [VAR=VALUE]..."
echo ""
echo " Invoke make_mfix without arguments for a guided build."
echo ""
echo " Configuration:"
echo ""
echo " -h, --help           Display this help and exit"
echo " -V, --version        Display version information and exit"
echo " -l, --long           Display all build options"
echo " -c, --clean          Remove previous build directory"
echo " -d, --default        Enable default build (serial, GCC compiler)"
echo ""
echo " Optional Features:"
echo ""
echo " --force-recompile    Forced recompile of modified files"
echo " --serial             Disable shared and distributed memory parallel"
echo " --smp                Enable shared memory parallel (OpenMP)"
echo " --dmp                Enable distributed memory paralle (MPI)"
echo " --debug              Compile with debugging information (e.g., -g flag)"
echo ""
echo " --opt=level - Specify compiler optimization level"
echo "     O0     -   No optimization with debugging flags"
echo "     O1     -   O1 optimization - Low"
echo "     O2     -   O2 optimization - Moderate"
echo "     O3     -   O3 optimization - Aggressive"
echo ""
echo " --compiler=option - Specify compiler flag file. There are three"
echo "                     generic options that cover most user needs:"
echo "    gcc            - GCC compiler (gfortran)"
echo "    intel          - Intel compiler (ifort)"
echo "    portland       - Portland Groupe (pgf90)"
echo ""
echo " Specific hardware configuration files: This script will first check for"
echo " the specified file in the run directory, then within mfix/config directory."
echo ""
echo "    gcc_default.sh         - same as gcc"
echo "    intel_default.sh       - same as intel"
echo "    portland_default.sh    - same as portland"
echo ""
echo "    SBEUC system files:"
echo "      sbeuc_gcc-46.sh      - GCC 4.6 and Open MPI 1.5.5"
echo "      sbeuc_intel-131.sh   - Intel 13.1 and Intel MPI"
#echo ""
#echo "    OLCF system files:"
#echo "      olcf_xt4_pg.sh       - Cray XT4 - Portland Group (pgf90)"
#echo ""
#echo "    ALCF system files:"
#echo "      alcf_bgp_ibm.sh      - IBM BG/P (xlf90)"
#echo "      alcf_bgq_ibm.sh      - IBM BG/Q (xlf90)"
echo ""
echo " Advanced Features:"
echo ""
echo " --mpi=PATH           MPI instalation directory"
echo " --mpi_include=PATH   MPI include directory"
echo " --mpi_lib=PATH       MPI library directory"
echo ""
echo " --enable-tau         Enable TAU profiling"
echo ""
echo ""
echo " Some influential environment variables:"
echo "  FORTRAN_CMD"
echo "  FORT_FLAGS"
echo ""
echo " ======================================================================"
echo ""
echo ""
exit
}

# Display version.
#-------------------------------------------------------------------------->>
show_version

# If version information is explictly requested, exit.
for arg in $input; do
  case ${arg} in
    "--version"|"-V") exit;;
    *) echo "" > /dev/null 2>&1;
  esac
done



# Loop the the argument list.
for arg in $input; do

  case ${arg} in

# Display help.
#-------------------------------------------------------------------------->>
    "--help"|"-h")
      show_usage;;


# Clean out the last build.
#-------------------------------------------------------------------------->>
    "--clean" | "-c")
      log=${MFIX_SRC}/${MAKEFILE}
      if ! test -f ${log}; then
        echo "Unable to locate previous build log."
        echo "Cannot remove previous build directory"
        exit
      fi
      build_dir=$(grep "DPO=" ${log} | cut -d "=" -f2)
      echo "Cleaning previous build: ${build_dir}"

      if test ! -d ${build_dir}; then
        echo "Nothing to clean"
      else
        /bin/rm -rf ${build_dir}
        if test $? = 0; then
          echo "Previous build directory removed."
        else
          echo "Error removing previous build directory."
          exit
        fi
      fi;;

# Specify default compile options.
#-------------------------------------------------------------------------->>
    "--default"|"-d" )

      EXPERT=0
      FORCE_COMPILE=0
      COMP_FILE="gnu_default.sh"

      if test -z ${USE_SMP}; then USE_SMP=0; fi
      if test -z ${USE_DMP}; then USE_DMP=0; fi
      REQ_MODE=0

      OPT=3
      REQ_OPT=0;;


# Show all available options.
#-------------------------------------------------------------------------->>
    "--long" | "-l" )
      echo "All compilation options will be shown"
      EXPERT=1;;


# A recompile command is executed. Invoke the make file
# that was previously generated.
#-------------------------------------------------------------------------->>
    "--repeat" | "-r" )
      if test ! -f "${MFIX_SRC}/${MAKEFILE}"; then
        echo "Cannot locate Makefile in model directory."
        echo "Unable to repeat last compile."
        exit
      fi
      mfile=${MFIX_SRC}/${MAKEFILE}
      echo "Using last compile settings."
      DPO=$(grep "DPO=" ${mfile} | cut -d "=" -f2)
      AUTOCOMPILE=1;;

# Option for multiple make jobs
#-------------------------------------------------------------------------->>
    "-j" ) MAKE_ARGS="-j";;

# Enable specify debug flags.
#-------------------------------------------------------------------------->>
    "--debug" ) USE_DEBUG="1";;

# Code coverage flag
#-------------------------------------------------------------------------->>
    "--codecov" ) USE_CODECOV="1";;

# Specify optimization level.
#-------------------------------------------------------------------------->>
    "--opt="* )
      OPT=$(echo ${arg} | cut -d "=" -f2)
      case ${OPT} in
        O0) OPT=0; USE_DEBUG=1;;
        O1) OPT=1;;
        O2) OPT=2;;
        O3) OPT=3;;
        *) echo "Error unknown optimization level: ${arg}"
           echo "First character is upper case letter O, not zero"
           echo "Aborting."
           exit;;
      esac
      echo "Specified optimization level: ${OPT}"
      REQ_OPT=0;;


# Enable SMP model
#-------------------------------------------------------------------------->>
    "--smp" )
      USE_SMP=1
      echo "SMP build mode enabled."
      if test -z ${USE_DMP}; then USE_DMP=0; fi
      REQ_MODE=0;;


# Enable DMP model
#-------------------------------------------------------------------------->>
    "--dmp" )
      USE_DMP=1
      echo "DMP build mode enabled."
      if test -z ${USE_SMP}; then USE_SMP=0; fi
      REQ_MODE=0;;


# Enable Serial model
#-------------------------------------------------------------------------->>
    "--serial" )
      USE_SMP=0
      USE_DMP=0
      REQ_MODE=0;;

# Enable MIC model
#-------------------------------------------------------------------------->>
    "--mic" )
      USE_MIC=1
      REQ_MODE=0;;

# Bypass chemical reaction preprocessing.
#-------------------------------------------------------------------------->>
    "--skip-rxn" )
      if test ${REQ_RXNS} = 1; then
        echo "WARNING: Skipping chemical reaction preprocessing."
      fi
      REQ_RXNS=0;;


# Compiler selection
#-------------------------------------------------------------------------->>
    "--compiler="*)
      compiler=$(echo ${arg} | cut -d "=" -f2)
      case ${compiler} in
        gcc) COMP_FILE=${MFIX_CONFIG}/compilers/gcc_default.sh;;
        intel) COMP_FILE=${MFIX_CONFIG}/compilers/intel_default.sh;;
        portland) COMP_FILE=${MFIX_CONFIG}/compilers/portland_default.sh;;
        *) comp=$(echo ${arg} | cut -d "=" -f2)
          if test -f ${RUN_DIR}/${comp}; then
            COMP_FILE=${RUN_DIR}/${comp}
          elif test -f ${MFIX_CONFIG}/compilers/${comp}; then
            COMP_FILE=${MFIX_CONFIG}/compilers/${comp}
          elif test -f ${MFIX_CONFIG}/${comp}; then
            COMP_FILE=${MFIX_CONFIG}/${comp}
          elif test -f ${comp}; then
            COMP_FILE=${comp}
          else
            echo "  Error: Unable to locate compiler file!"
            echo "   >>> ${arg}"
            exit
          fi ;;
      esac
      REQ_COMP=0;;


# Enable force source file recompiling.
#-------------------------------------------------------------------------->>
    "--force-recompile"* )
      FORCE_COMPILE=1;;


# Use specified build directory.
#-------------------------------------------------------------------------->>
    "--build="*)
      dir=$(echo ${arg} | cut -d "=" -f2)
      if test ! -d ${dir}; then
        echo "  Specified build directory not found!"
        echo "   >>> ${dir}"
        exit
      fi
      cd ${dir}
      set `pwd` ; DPO_BASE=$1
      cd ${MFIX_SRC}
      echo "Build directory: ${DPO_BASE}";;


# Use specified executable name. (mfix.exe is the default)
#-------------------------------------------------------------------------->>
    "--exe="*)
      EXEC_FILE=$(echo ${arg} | cut -d "=" -f2)
      if test -z ${EXEC_FILE}; then
        echo "Specified executable name is empty!"
        echo "Aborting make_mfix."
        exit
      fi
      echo "User specified executable name: ${EXEC_FILE}";;


# MPI installation directory. This folder should contain the include
# and library directories.
#-------------------------------------------------------------------------->>
    "--mpi="*)
      dir=$(echo ${arg} | cut -d "=" -f2)
      if test ! -d ${dir}; then
        echo "   Specified MPI path not found!"
        echo "   >>> ${dir}"
        exit
      fi
      cd ${dir}
      set `pwd` ; MPI_PATH=$1
      cd ${MFIX_SRC}
      echo "MPI path: ${MPI_PATH}";;


# User specified path to MPI include directory.
#-------------------------------------------------------------------------->>
    "--mpi_include="*)
      dir=$(echo ${arg} | cut -d "=" -f2)
      if test ! -d ${dir}; then
        echo "   Specified MPI include path not found!"
        echo "   >>> ${dir}"
        exit
      fi
      cd ${dir}
      set `pwd` ; MPI_INCLUDE_PATH=$1
      if test ! -f "mpif.h"; then
        echo "  Specified mpi_include does not contain mpif.h"
        echo "   >>> $PWD"
        exit
      fi
      cd ${MFIX_SRC}
      echo "MPI include path: ${MPI_INCLUDE_PATH}";;


# MPI Library directory
#-------------------------------------------------------------------------->>
    "--mpi_lib="*)
      dir=$(echo ${arg} | cut -d "=" -f2)
      if test ! -d ${dir}; then
        echo "  Specified mpi library path not found!"
        echo "   >>> ${dir}"
        exit
      fi
      cd ${dir}
      set `pwd` ; MPI_LIB_PATH=$1
      cd ${MFIX_SRC}
      echo "MPI library path: ${MPI_LIB_PATH}";;


# Enable TAU profiling
#-------------------------------------------------------------------------->>
    "--enable-tau" )
      if test -z ${TAUROOT}; then
        echo "  Fatal Error: TAUROOT not set!"
        exit
      fi
      USE_TAU=1
      DPO=${DPO}_TAU;;


# Enable NetCDF output
#-------------------------------------------------------------------------->>
    "--enable-netcdf" )
      if test -z ${NETCDF_HOME}; then
        echo "  Fatal Error: NETCDFROOT not set!"
        exit
      fi
      USE_NETCDF=1
      DPO=${DPO}_NCDF;;


# An unknown build command was specified. Print the usage
# information and exit.
#-------------------------------------------------------------------------->>
    *)echo "Unknown flag: ${arg}"
      show_usage;;

  esac
done




#echo "Forcing a stop in parse_args.sh"
#exit
