################################################################################
#
# DO NOT FORGET TO change module for gcc on jaguar
# module swap PrgEnv-pgi PrgEnv-gnu
#
#
####################                                        ####################
####################   Compile MFIX to Ensight Translator   ####################
####################                                        ####################
################################################################################
#
# Purpose
#    Compile the C code that translates MFIX files to Ensight 6 files.
# Original
#    Author: Kent Eschenberg eschenbe@psc.edu
#    Organization: Pittsburgh Supercomputing Center, CMU
#    Date: 2002/1/29 10:41
# Current
#    $Source$
#    $Author$
#    $Revision$
#    $Date$

MAINH = m2e.h
MAINC = m2e.c
MAINO = m2e.o
MAIN = m2e

####
####  Pieces
####

MRH = MfixReader.h
MRC = MfixReader.c
MRO = MfixReader.o

CALCC = calc.c
CALCO = calc.o

GH = cell2geom.h
GC = cell2geom.c
GO = cell2geom.o

KH = kLib.h
KC = kLib.c
KO = kLib.o

BH = bio.h
BC = libbio.c
BO = libbio.o

####
####  Compile and link options
####
#COMP = /opt/gcc/4.2.0.quadcore/bin/gcc
COMP = gcc
DEBUG = -g
#DEBUG = -O -Wuninitialized
WARN1 = -W -Wcomment -Wformat -Wimplicit -Wparentheses -Wreturn-type -Wswitch
WARN2 = -Wunused -Wsign-compare -Waggregate-return
WARN3 = -Wstrict-prototypes -Wmissing-prototypes -Wredundant-decls
WARN4 = -Wnested-externs -Wlong-long
INCDIR = -I.
CCFLAGS = -c $(DEBUG) $(WARN1) $(WARN2) $(WARN3) $(WARN4) $(INCDIR)
ITCHES = $(MAINH) $(MRH) $(GH) $(KH) $(BH)

####
####  Dependencies
####

$(MAIN): $(MAINO) $(CALCO) $(MRO) $(GO) $(KO) $(BO)
	$(COMP) -o $(MAIN) $(CALCO) $(MAINO) $(MRO) $(GO) $(KO) $(BO) -lm

$(MAINO): $(MAINC) $(ITCHES)
	$(COMP) $(CCFLAGS) $(MAINC) 

$(CALCO): $(CALCC) $(ITCHES)
	$(COMP) $(CCFLAGS) -o $(CALCO) $(CALCC)

$(MRO): $(MRC) $(MRH)
	$(COMP) $(CCFLAGS) -o $(MRO) $(MRC)

$(GO): $(GC) $(GH)
	$(COMP) $(CCFLAGS) -o $(GO) $(GC)

$(KO): $(KC) $(KH)
	$(COMP) $(CCFLAGS) -o $(KO) $(KC)

$(BO): $(BC) $(BH)
	$(COMP) $(CCFLAGS) -o $(BO) $(BC)

clean:
	-/bin/rm -f *.o *~ $(MAIN) 
#############################  End of makefile  ################################
