include ../../../config/makefile.h
include ../../DEFS

ifeq ($(peigs_CODEOBJ),DBLE)
FOBJ = damax.o \
       dlaebz2.o \
       dstebz3.o \
       dspgv2.o \
       dspevx2.o \
       heapsort.o \
       neblw2.o
endif

ifeq ($(peigs_CODEOBJ),SINGLE)
FOBJ = samax.o \
       slaebz2.o \
       sstebz3.o \
       sheapsort.o \
       sneblw2.o
endif


EOBJ	= xerbl2.o ilaenv1.o menode.o pairup.o xstop.o


all:		$(FOBJ) $(EOBJ)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:
		-rm -f *.o

realclean:
		clean

dstebz1.o:	dstebz1.f
		$(peigs_FC) -c $(peigs_FOPT) dstebz1.f  

pairup.o:	pairup.f
		$(peigs_FC) -c $(peigs_FOPT) pairup.f  

neblw1.o:	neblw1.f
		$(peigs_FC) -c $(peigs_FOPT) neblw1.f  

dlaebz1.o:	dlaebz1.f
		$(peigs_FC) -c $(peigs_FOPT) dlaebz1.f  

ilaenv1.o:	ilaenv1.f
		$(peigs_FC) -c ilaenv1.f  

heapsort.o:	heapsort.f
		$(peigs_FC) -c heapsort.f  

menode.o:	menode.f
		$(peigs_FC) -c menode.f  

xerbl2.o:	xerbl2.f
		$(peigs_FC) -c $(peigs_FOPT) xerbl2.f  

unif.o:		unif.f
		$(peigs_FC) -c $(peigs_FOPT) unif.f  


.c.o:;	$(peigs_CC) $(OPTC) -c $<
.f.o:;	$(peigs_FC) $(FOPT) -c $<


