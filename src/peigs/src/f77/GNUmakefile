#
# $Id: GNUmakefile,v 1.21 2004-10-20 17:52:29 edo Exp $
#

include ../../../config/makefile.h
include ../../DEFS
NOCASENS = CYGNUS MACX CYGWIN INTERIX
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
	dlarnv.o \
	dlasrt.o \
	dlasq3.o \
	dlasq4.o \
	dlasrt.o \
	dlaebz.o \
	dlascl.o \
	types1.o \
	dlasq3a.o \
	dlasq4a.o \
	dlasq5a.o \
	dlasq5.o 
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
	dlarnv.o \
	dlasrt.o \
	dlasq3.o \
	dlasq4.o \
	dlasrt.o \
	dlaebz.o \
	dlascl.o \
	types1.o \
	dlasq3a.o \
	dlasq4a.o \
	dlasq5a.o 
endif


EOBJ	= xerbl2.o ilaenv1.o menode.o pairup.o xstop.o

all:		$(FOBJ) $(EOBJ)

.EXPORT_ALL_VARIABLES:

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
		-rm -f dgetavec.f dlaebz.f dlaebz2.f dlagtf1.f \
		dlas2.f  \
		dlasq3.f dlasq4.f dlasr.f \
		dlasrt.f dlassq.f dlasv2.f dlaswp.f \
		dlasy2.f dlasyf.f dspevx2.f dstebz3.f \
		neblw2.f sneblw2.f sstebz3.f
endif
		-rm -f F*.f

.f.o:; $(peigs_FC) -c $<
.F.o:; $(peigs_FC) -c $<
.c.o:; $(peigs_CC) -c $<
