##
## Programs and options...
##
include /home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/make-variables
TOPDIR=/home/bdotson/applications/visit/visit/1.5.2/linux-intel


##
## Libraries and includes...
##
VTK_INCLUDE= \
  -I/home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/vtk \
  -I/home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/vtk/Common \
  -I/home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/vtk/Filtering \
  -I/home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/vtk/Graphics \
  -I/home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/vtk/Hybrid \
  -I/home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/vtk/IO \
  -I/home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/vtk/Imaging \
  -I/home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/vtk/Rendering
MOC=/home/bdotson/applications/visit/visit/1.5.2/linux-intel/bin/moc
CXXFLAGS=$(CXXFLAGSORIG)
CPPFLAGS=$(CPPFLAGSORIG) $(VTK_INCLUDE) -I. -I/home/bdotson/applications/visit/visit/1.5.2/linux-intel/include -I/home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/visit
LDFLAGS=$(LDFLAGSORIG)  -L/home/bdotson/.visit/linux-intel/plugins/databases

##
## Files...
##
PLUGINDIR=databases
PLUGINNAME=MFIXDatabase

WIDGETS=
ISRC=MFIXPluginInfo.C
COMMONSRC=MFIXPluginInfo.C MFIXCommonPluginInfo.C
MSRC=MFIXMDServerPluginInfo.C avtMFIXFileFormat.C
ESRC=MFIXEnginePluginInfo.C avtMFIXFileFormat.C


SRC=$(ISRC) $(COMMONSRC) $(MSRC) $(ESRC)

SHLIB_FORCED=

ELIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS)
MLIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS)
ILIBS=
MLIBS=-lpipeline_ser -ldatabase_ser -lmir_ser -lavtmath_ser -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk  -L/home/bdotson/.visit/linux-intel/plugins/databases -lI$(PLUGINNAME) $(MLIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
ESERLIBS=-lpipeline_ser -ldatabase_ser -lmir_ser -lavtmath_ser -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk  -L/home/bdotson/.visit/linux-intel/plugins/databases -lI$(PLUGINNAME) $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
EPARLIBS=-lpipeline_par -ldatabase_par -lmir_par -lavtmath_par -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk  -L/home/bdotson/.visit/linux-intel/plugins/databases -lI$(PLUGINNAME) $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS) $(SHLIB_MPI_LIBS)

IDSO=/home/bdotson/.visit/linux-intel/plugins/databases/libIMFIXDatabase.so
MDSO=/home/bdotson/.visit/linux-intel/plugins/databases/libMMFIXDatabase.so
ESERDSO=/home/bdotson/.visit/linux-intel/plugins/databases/libEMFIXDatabase_ser.so
EPARDSO=/home/bdotson/.visit/linux-intel/plugins/databases/libEMFIXDatabase_par.so


DISTRIB=

##
## Derived objects
##
IOBJ=$(ISRC:.C=.o)
GOBJ=$(COMMONSRC:.C=.o) $(GSRC:.C=.o)
SOBJ=$(COMMONSRC:.C=.o) $(SSRC:.C=.o)
VOBJ=$(COMMONSRC:.C=.o) $(VSRC:.C=.o)
MOBJ=$(COMMONSRC:.C=.o) $(MSRC:.C=.o) $(MSPECIFICSRC:.C=_mds.o)
ESEROBJ=$(COMMONSRC:.C=.o) $(ESRC:.C=.o)
EPAROBJ=$(COMMONSRC:.C=.o) $(ESRC:.C=_par.o)

MOCSRC = $(WIDGETS:.h=_moc.C)
MOCOBJ = $(MOCSRC:.C=.o)

##
## Standard targets...
##
all: message $(IDSO) $(MDSERVERLIB) $(ENGINELIBSER) $(ENGINELIBPAR)

clean:
	$(RM) $(IOBJ) $(COMMONOBJ)
	$(RM) $(GOBJ) $(SOBJ) $(VOBJ) $(MOBJ) $(ESEROBJ) $(EPAROBJ)
	$(RM) $(MOCSRC) $(MOCOBJ)
	$(RM) $(IDSO) $(GDSO) $(SDSO) $(VDSO) $(MDSO) $(ESERDSO) $(EPARDSO)
	$(RM) $(JAVAOBJ)

##
## Other targets...
##
message:
	@echo
	@echo "****************************************************************************"
	@echo "*** Building MFIX Database Plugin"
	@echo "****************************************************************************"

##
## moc
##
$(MOCSRC) or_no_widgets: $(WIDGETS)
	@rm -f $@
	$(MOC) $(@:_moc.C=.h) > $@


##
## Automatic dependency stuff
##
include /home/bdotson/applications/visit/visit/1.5.2/linux-intel/include/make-targets
