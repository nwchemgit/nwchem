include ../../DEFS

ifeq ($(CODEOBJ),DBLE)
FOBJ = damax.o \
       dlaebz2.o \
       dstebz3.o \
       dspgv2.o \
       dspevx2.o \
       heapsort.o \
       neblw2.o
endif

ifeq ($(CODEOBJ),SINGLE)
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
		$(FC) -c $(FCNFLAGS) dstebz1.f  

pairup.o:	pairup.f
		$(FC) -c $(FCNFLAGS) pairup.f  

neblw1.o:	neblw1.f
		$(FC) -c $(FCNFLAGS) neblw1.f  

dlaebz1.o:	dlaebz1.f
		$(FC) -c $(FCNFLAGS) dlaebz1.f  

ilaenv1.o:	ilaenv1.f
		$(FC) -c ilaenv1.f  

heapsort.o:	heapsort.f
		$(FC) -c heapsort.f  

menode.o:	menode.f
		$(FC) -c menode.f  

xerbl2.o:	xerbl2.f
		$(FC) -c $(FCNFLAGS) xerbl2.f  

unif.o:		unif.f
		$(FC) -c $(FCNFLAGS) unif.f  


.c.o:;	$(CC) $(OPTC) -c $<

.f.o:;	$(FC) $(FOPT) -c $<

