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
	dlascl.o
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
	dlascl.o
endif


EOBJ	= xerbl2.o ilaenv1.o menode.o pairup.o xstop.o


all:		$(FOBJ) $(EOBJ)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean



		-rm -f F*.f

realclean:
		clean

dgetavec.o:	dgetavec.F
		$(peigs_CPP) < dgetavec.F >  dgetavec.f
		$(peigs_FC) -c $(F77NFLAGS) dgetavec.f

dlaebz.o:	dlaebz.F
		$(peigs_CPP) < dlaebz.F >  dlaebz.f
		$(peigs_FC) -c $(F77NFLAGS) dlaebz.f

dlaebz2.o:	dlaebz2.F
		$(peigs_CPP) < dlaebz2.F >  dlaebz2.f
		$(peigs_FC) -c $(F77NFLAGS) dlaebz2.f

dlagtf1.o:	dlagtf1.F
		$(peigs_CPP) < dlagtf1.F >  dlagtf1.f
		$(peigs_FC) -c $(F77NFLAGS) dlagtf1.f

dlas2.o:	dlas2.F
		$(peigs_CPP) < dlas2.F >  dlas2.f
		$(peigs_FC) -c $(F77NFLAGS) dlas2.f

dlasq1.o:	dlasq1.F
		$(peigs_CPP) < dlasq1.F >  dlasq1.f
		$(peigs_FC) -c $(F77NFLAGS) dlasq1.f

dlasq2.o:	dlasq2.F
		$(peigs_CPP) < dlasq2.F >  dlasq2.f
		$(peigs_FC) -c $(F77NFLAGS) dlasq2.f

dlasq3.o:	dlasq3.F
		$(peigs_CPP) < dlasq3.F >  dlasq3.f
		$(peigs_FC) -c $(F77NFLAGS) dlasq3.f

dlasq4.o:	dlasq4.F
		$(peigs_CPP) < dlasq4.F >  dlasq4.f
		$(peigs_FC) -c $(F77NFLAGS) dlasq4.f

dlasr.o:	dlasr.F
		$(peigs_CPP) < dlasr.F >  dlasr.f
		$(peigs_FC) -c $(F77NFLAGS) dlasr.f

dlasrt.o:	dlasrt.F
		$(peigs_CPP) < dlasrt.F >  dlasrt.f
		$(peigs_FC) -c $(F77NFLAGS) dlasrt.f

dlassq.o:	dlassq.F
		$(peigs_CPP) < dlassq.F >  dlassq.f
		$(peigs_FC) -c $(F77NFLAGS) dlassq.f

dlasv2.o:	dlasv2.F
		$(peigs_CPP) < dlasv2.F >  dlasv2.f
		$(peigs_FC) -c $(F77NFLAGS) dlasv2.f

dlaswp.o:	dlaswp.F
		$(peigs_CPP) < dlaswp.F >  dlaswp.f
		$(peigs_FC) -c $(F77NFLAGS) dlaswp.f

dlasy2.o:	dlasy2.F
		$(peigs_CPP) < dlasy2.F >  dlasy2.f
		$(peigs_FC) -c $(F77NFLAGS) dlasy2.f


dlasyf.o:	dlasyf.F
		$(peigs_CPP) < dlasyf.F >  dlasyf.f
		$(peigs_FC) -c $(F77NFLAGS) dlasyf.f

dspevx2.o:	dspevx2.F
		$(peigs_CPP) < dspevx2.F >  dspevx2.f
		$(peigs_FC) -c $(F77NFLAGS) dspevx2.f

dstebz3.o:	dstebz3.F
		$(peigs_CPP) < dstebz3.F >  dstebz3.f
		$(peigs_FC) -c $(F77NFLAGS) dstebz3.f

neblw2.o:	neblw2.F
		$(peigs_CPP) < neblw2.F >  neblw2.f
		$(peigs_FC) -c $(F77NFLAGS) neblw2.f

sneblw2.o:	sneblw2.F
		$(peigs_CPP) < sneblw2.F >  sneblw2.f
		$(peigs_FC) -c $(F77NFLAGS) sneblw2.f

sstebz3.o:	sstebz3.F
		$(peigs_CPP) < sstebz3.F >  sstebz3.f
		$(peigs_FC) -c $(F77NFLAGS) sstebz3.f


.c.o:;	$(peigs_CC) $(OPTC) -c $<
.f.o:;	$(peigs_FC) $(OPT) -c $<
.F.o:;	$(peigs_FC) $(OPT) -c $<


clean:
		-rm -f *.o
		-rm -f *.s
		-rm -f dgetavec.f dlaebz.f dlaebz2.f dlagtf1.f \
		dlas2.f dlasq1.f dlasq2.f \
		dlasq3.f dlasq4.f dlasr.f \
		dlasrt.f dlassq.f dlasv2.f dlaswp.f \
		dlasy2.f dlasyf.f dspevx2.f dstebz3.f \
		neblw2.f sneblw2.f sstebz3.f
