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
	dlaebz.o

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
	dlaebz.o

endif


EOBJ	= xerbl2.o ilaenv1.o menode.o pairup.o xstop.o


all:		$(FOBJ) $(EOBJ)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:
		-rm -f *.o
		-rm -f *.s
		-rm -f dgetavec.f dlasq1.f dspevx2.f dstebz3.f \
			 neblw2.f sneblw2.f sstebz3.f
		-rm -f F*.f

realclean:
		clean

dgetavec.o:	dgetavec.F
		$(peigs_CPP) < dgetavec.F >  dgetavec.f
		$(peigs_FC) -c $(F77NFLAGS) dgetavec.f

dlasrt.o:	dlasrt.F
		$(peigs_CPP) < dlasrt.F >  dlasrt.f
		$(peigs_FC) -c $(F77NFLAGS) dlasrt.f


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

