#
# $Id: GNUmakefile,v 1.17 1999-10-28 17:12:01 d3g270 Exp $
#

include ../../../config/makefile.h
include ../../DEFS

ifeq ($(peigs_CODEOBJ),DBLE)
FOBJ = damax.o \
       dlaebz2.o \
       dstebz3.o \
       heapsort.o \
       neblw2.o\
       dspevx2.o \
	dgetavec.o \
	dgetavec2.o \
	dgetavec3.o \
	dlarnv.o \
	dlasq1.o \
	dlasrt.o \
	dlasq2.o \
	dlasq3.o \
	dlasq4.o \
	dlasrt.o \
	dlaebz.o \
	dlascl.o \
	types1.o \
	dlasq2a.o \
	dlasq3a.o \
	dlasq4a.o \
	dlasq5a.o 
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
	dlasq1.o \
	dlasrt.o \
	dlasq2.o \
	dlasq3.o \
	dlasq4.o \
	dlasrt.o \
	dlaebz.o \
	dlascl.o \
	types1.o \
	dlasq2a.o \
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

clean:
		-rm -f *.o
		-rm -f *.s
		-rm -f dgetavec.f dlaebz.f dlaebz2.f dlagtf1.f \
		dlas2.f dlasq1.f dlasq2.f \
		dlasq3.f dlasq4.f dlasr.f \
		dlasrt.f dlassq.f dlasv2.f dlaswp.f \
		dlasy2.f dlasyf.f dspevx2.f dstebz3.f \
		neblw2.f sneblw2.f sstebz3.f
		-rm -f F*.f


dlasq1.o:	dlasq1.F
		$(peigs_FC) -c -g $<

dlasq2.o:	dlasq2.F
		$(peigs_FC) -c -g $<

.f.o:; $(peigs_FC) -c $<
.F.o:; $(peigs_FC) -c $<
.c.o:; $(peigs_CC) -c $<
