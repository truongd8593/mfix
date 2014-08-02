#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT(mfix, 2014.1, develop@mfix.netl.doe.gov)

AC_MSG_NOTICE([==============================================================])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([       MFIX: Multiphase Flow with Intephase eXchanges])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([                      mfix.netl.doe.gov])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([                        Version ${MFIX_VERSION}])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([==============================================================])

AC_LANG(Fortran)

# Checks for programs.
#AC_PROG_AWK
AC_PROG_F77
AC_PROG_FC([gfortran ifort])
AC_PROG_RANLIB
AC_PROG_GREP
AC_PROG_SED
#AC_PROG_LN_S
#AC_PROG_MAKE_SET

AC_OPENMP()

# Checks for source files.

AC_CHECK_FILES(usr0.f usr1.f usr2.f usr3.f usr_mod.f usr_write_out0.f usr_rates.f mfix.dat)
AM_CONDITIONAL([USR0], [test "$ac_cv_file_usr0_f" = yes])
AM_CONDITIONAL([USR1], [test "$ac_cv_file_usr1_f" = yes])
AM_CONDITIONAL([USR2], [test "$ac_cv_file_usr2_f" = yes])
AM_CONDITIONAL([USR3], [test "$ac_cv_file_usr3_f" = yes])
AM_CONDITIONAL([USR_MOD], [test "$ac_cv_file_usr_mod_f" = yes])
AM_CONDITIONAL([USR_WRITE_OUT0], [test "$ac_cv_file_usr_write_out0_f" = yes])
AM_CONDITIONAL([USR_RATES], [test "$ac_cv_file_usr_rates_f" = yes])

if test "$ac_cv_file_usr_rates_f" = yes; then
   if test "$ac_cv_file_mfix_dat" = no; then
	echo "If usr_rate.f is present then mfix.dat must be present"
	exit -1
	fi
fi

AC_CHECK_FILES(des/usr0_des.f des/usr1_des.f des/usr2_des.f des/usr3_des.f des/usr4_des.f des/usr_rates_des.f dqmom/usr_dqmom.f)
AM_CONDITIONAL([DES_USR0_DES], [test "$ac_cv_file_des_usr0_des_f" = yes])
AM_CONDITIONAL([DES_USR1_DES], [test "$ac_cv_file_des_usr1_des_f" = yes])
AM_CONDITIONAL([DES_USR2_DES], [test "$ac_cv_file_des_usr2_des_f" = yes])
AM_CONDITIONAL([DES_USR3_DES], [test "$ac_cv_file_des_usr3_des_f" = yes])
AM_CONDITIONAL([DES_USR4_DES], [test "$ac_cv_file_des_usr4_des_f" = yes])
AM_CONDITIONAL([DES_USR_RATES_DES], [test "$ac_cv_file_des_usr_rates_des_f" = yes])
AM_CONDITIONAL([DQMOM_USR_DQMOM], [test "$ac_cv_file_dqmom_usr_dqmom_f" = yes])

# Checks for typedefs, structures, and compiler characteristics.
# AC_HEADER_STDBOOL
# AC_TYPE_SIZE_T

# Checks for library functions.
#AC_CHECK_FUNCS([memset])

AC_CONFIG_SRCDIR([model/mfix.f])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([subdir-objects foreign -Wall -Wno-portability])

AC_ARG_ENABLE([smp],
[  --enable-smp    Turn on OpenMP support],
[case "${enableval}" in
  yes) smp=true ;;
  no)  smp=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-smp]) ;;
esac],[smp=false])
AM_CONDITIONAL([OPENMP], [test x$smp = xtrue])

AC_ARG_ENABLE([dmp],
[  --enable-dmp    Turn on MPI support],
[case "${enableval}" in
  yes) dmp=true ;;
  no)  dmp=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-dmp]) ;;
esac],[dmp=false])
AM_CONDITIONAL([MPI], [test x$dmp = xtrue])

AC_FC_SRCEXT(f)
AC_FC_FREEFORM()
AC_FC_LINE_LENGTH([unlimited])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT

AC_MSG_NOTICE([==============================================================])
AC_MSG_NOTICE([Now run "make" to build mfix])
AC_MSG_NOTICE([==============================================================])
