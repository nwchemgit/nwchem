include ../../config/makefile.h
include ../DEFS

MXOBJ	=	mxbrod.o \
		mxcombv1.o \
		mxsubs.o \
		sum.o \
		maxx.o \
		peigs_dgop.o

all:		$(MXOBJ)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:		
		-rm -f *.o
		-rm -f mxbrod.f
		-rm -f mxcombv1.f
		-rm -f mxsubs.f
		-rm -f sum.f peigs_dgop.f

realclean:
		clean

.SUFFIXES:	.o .s .F .f .c

.c.o:;	$(peigs_CC) $(peigs_COPT) -c $<
.f.o:;	$(peigs_FC) $(peigs_FOPT) -c $<


mxbrod.o:	mxbrod.cpp
		-rm -f mxbrod.f
		$(peigs_CPP) < mxbrod.cpp | sed '/^\#/D' > mxbrod.f
		$(peigs_FC) $(peigs_FOPT) -c mxbrod.f

mxsubs.o:	mxsubs.cpp
		-rm -f mxsubs.f
		$(peigs_CPP) < mxsubs.cpp | sed '/^\#/D' > mxsubs.f
		$(peigs_FC) $(peigs_FOPT) -c mxsubs.f

mxcombv1.o:	mxcombv1.cpp
		-rm -f mxcombv1.f
		$(peigs_CPP) < mxcombv1.cpp | sed '/^\#/D' > mxcombv1.f
		$(peigs_FC) $(peigs_FOPT) -c mxcombv1.f

sum.o:		sum.cpp
		-rm -f sum.f
		$(peigs_CPP) < sum.cpp | sed '/^\#/D' > sum.f
		$(peigs_FC) $(peigs_FOPT) -c sum.f

peigs_dgop.o:	peigs_dgop.cpp
		-rm -f peigs_dgop.f
		$(peigs_CPP) < peigs_dgop.cpp | sed '/^\#/D' > peigs_dgop.f
		$(peigs_FC) $(peigs_FOPT) -c peigs_dgop.f




