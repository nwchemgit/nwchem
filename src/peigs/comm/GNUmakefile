include ../DEFS

MXOBJ	=	mxbrod.o \
			mxcombv1.o \
			mxsubs.o \
			sum.o \
			maxx.o

all:		$(MXOBJ)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:		
		-rm -f *.o
		-rm -f mxbrod.f
		-rm -f mxcombv1.f
		-rm -f mxsubs.f
		-rm -f sum.f

realclean:
		clean

.SUFFIXES:	.o .s .F .f .c

.c.o:;	$(peigs_CC) $(peigs_COPT) -c $<
.f.o:;	$(peigs_FC) $(peigs_FOPT) -c $<

#ifdef EXPLICITF
.F.o:	
	@echo Converting $*.F '->' $*.f
	@$(peigs_FCONVERT)
	$(FC) -c $(FFLAGS) $*.f
	@$(RM) $*.f

.F.f:
	@echo Converting $*.F '->' $*.f
	@$(peigs_FCONVERT)

#else
.F.o:;  $(peigs_FC) $(peigs_FOPT) -c $<
#endif






