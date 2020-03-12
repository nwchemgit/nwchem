#
# $Id$
#

include ../../../config/makefile.h
include ../../DEFS
NOCASENS = CYGNUS MACX64 MACX CYGWIN CYGWIN64
ifneq ($(peigs_TARGET),$(findstring $(peigs_TARGET), $(NOCASENS)))
	REMOVESF=TRUE
endif

ifeq ($(peigs_CODEOBJ),DBLE)
FOBJ = inverm.o \
	damax.o \
       dlaebz2.o \
       dstebz3.o \
       heapsort.o \
       neblw2.o\
       dspevx2.o \
	dgetavec.o \
	dgetavec2.o \
	dgetavec3.o \
	types1.o 
endif

ifeq ($(peigs_CODEOBJ),SINGLE)
FOBJ = samax.o \
       slaebz2.o \
       sstebz3.o \
       sheapsort.o \
       sneblw2.o \
       dspevx2.o \
	dgetavec.o \
	dgetavec2.o \
	dgetavec3.o \
	types1.o 
endif


EOBJ	= xerbl2.o  menode.o pairup.o xstop.o

all:		$(FOBJ) $(EOBJ)

#.EXPORT_ALL_VARIABLES:

.PHONY :	clean


		-rm -f F*.f

realclean:
		clean

#
# Note that under Cygwin the rm commands are case-insensitive, so the
# commands below intended to remove intermediate .f files actually remove
# the original .F files!  Since intermediate files are not generated for
# Cygwin, this is avoided by commenting out the rm statements.
#
clean:
		-rm -f *.o
		-rm -f *.s
ifdef REMOVESF
		-rm -f dgetavec.f dlaebz2.f  \
		dspevx2.f dstebz3.f \
		neblw2.f sneblw2.f sstebz3.f
endif
		-rm -f F*.f

 CNFDIR := $(NWCHEM_TOP)/src/config
.PHONY:	 64_to_32 32_to_64
64_to_32:
	$(CNFDIR)/64_to_32 types.f inverm.F \
	 dgetavec.F dgetavec2.f dgetavec3.f dgetavec4.F \
	dptsvx2.f dspevx2.F

32_to_64:
	$(CNFDIR)/32_to_64 types.f inverm.F \
	 dgetavec.F dgetavec2.f dgetavec3.f dgetavec4.F \
	dptsvx2.f dspevx2.F

.f.o:; $(peigs_FC) -c $<
.F.o:; $(peigs_FC) -c $<
.c.o:; $(peigs_CC) -c $<
