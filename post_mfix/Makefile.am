################################################################
#                                                              #
#      Automake input file                                     #
#                                                              #
################################################################

F77 = $(FC)
F77COMPILE = $(FC) $(AM_FCFLAGS) $(FCFLAGS)
FCCOMPILE  = $(FC) $(AM_FCFLAGS) $(FCFLAGS)

# Fortran preprocessor is used (-cpp)
AM_FCFLAGS = -cpp
AM_LDFLAGS =
FC_MODEXT=mod
CLEANFILES=*.$(FC_MODEXT)
MODDIR=$(top_builddir)/model
AM_FCFLAGS += $(FC_MODINC)$(top_srcdir)/model/include
AM_FCFLAGS += $(FC_MODINC)$(top_srcdir)/post_mfix/include
AM_FCFLAGS += $(FC_MODINC)$(MODDIR)

.DEFAULT_GOAL := postmfix

################################################################
#                                                              #
#      main exectuable:  postmfix                              #
#                                                              #
################################################################
bin_PROGRAMS = postmfix

#  SOURCES are fully defined further below
postmfix_SOURCES = main_f.f mpi_utility_mod.f
postmfix_LDADD = ../model/ODEPACK.$(OBJEXT) $(top_builddir)/libmfix.a

# build main directory first
main_f.f: $(top_builddir)/libmfix.a

mpi_utility_mod.f: ../model/dmp_modules/mpi_utility_mod.f
	cp $^ $@

$(top_builddir)/libmfix.a:
	$(MAKE) -C $(top_srcdir) libmfix.a

################################################################
#                                                              #
#  Use script "fortran-depcomp" to generate makefile rules     #
#                                                              #
################################################################

postmodule.deps: $(postmfix_SOURCES)
	SRCDIR="$(srcdir)" OBJEXT="$(OBJEXT)" MODEXT="$(FC_MODEXT)" MODDIR="$(MODDIR)" \
	SED="$(SED)" GREP="$(GREP)" $(SHELL) $(top_srcdir)/build-aux/fortran-depcomp $^ > $@
	SRCDIR="$(srcdir)" OBJEXT="$(OBJEXT)" MODEXT="$(FC_MODEXT)" MODDIR="$(MODDIR)" \
	SED="$(SED)" GREP="$(GREP)" $(SHELL) $(top_srcdir)/build-aux/update-usrdeps $@
	@for main_mod in $(${GREP} -o -E "\S+\.${MODEXT}" $@); do \
	    echo "${main_mod}: main" >> $@; \
	done;

# Use generated module dependency rules
$(eval -include postmodule.deps)
CLEANFILES += postmodule.deps

################################################################
#                                                              #
#      OpenMP and MPI support                                  #
#                                                              #
#      --enable-smp and --enable-dmp options in configure.ac   #
#                                                              #
################################################################

if OPENMP
AM_FCFLAGS += $(OPENMP_FCFLAGS)
AM_LDFLAGS += $(OPENMP_FCFLAGS)
endif

################################################################
#                                                              #
#      postmfix                                                #
#                                                              #
################################################################

postmfix_SOURCES += \
	any_more_data.f \
	calc_cell2.f \
	calc_corr_01.f \
	calc_corr_type_1.f \
	calc_distance.f \
	calc_dpody.f \
	calc_ep_g.f \
	calc_quantities.f \
	calc_ro_g.f \
	calc_vol.f \
	correl_mod.f \
	deallocate_arrays.f \
	examine_data.f \
	f_init_data.f \
	file_handle.f \
	finit.f \
	flow_gx.f \
	flow_gy.f \
	flow_gz.f \
	flow_sx.f \
	flow_sy.f \
	flow_sz.f \
	gas_flux.f \
	get_file_name.f \
	get_file_status.f \
	get_index.f \
	get_location.f \
	get_mu_s.f \
	get_same_time.f \
	get_selection.f \
	get_substr.f \
	granular_qty.f \
	header_main.f \
	ik_avg.f \
	ik_avg_out.f \
	include/xforms.inc \
	interp_res.f \
	ornl_corr.f \
	ornl_corr_c.f \
	ornl_filt.f \
	ornl_filt_c.f \
	ornl_ft.f \
	ornl_ft_c.f \
	ornl_header.f \
	ornl_pca.f \
	ornl_stats.f \
	ornl_stats_c.f \
	ornl_sym.f \
	ornl_util.f \
	ornl_zone.f \
	out_from_res.f \
	out_from_spx.f \
	out_spec_time.f \
	out_time.f \
	paralleldata_mod.f \
	post3d_mod.f \
	post_precision_mod.f \
	print_out.f \
	read_spx0.f \
	read_spx1.f \
	res_from_spx.f \
	seek_time.f \
	select_spx_rec.f \
	set_dollar.f \
	set_read_spx.f \
	sol_flux.f \
	strcmp.f \
	streqs.f \
	time_avg.f \
	usr_init_namelist.f \
	usr_input.f \
	usr_input_mod.f \
	usr_post.f \
	usr_write_out1.f
