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

realclean:
		clean

mxsubs.o : 	mxsubs.cpp
	       	rm -f mxsubs.f
		$(peigs_CPP) < mxsubs.cpp > mxsubs.f
		$(peigs_FC) -c mxsubs.f

mxbrod.o :	mxbrod.cpp
	       	rm -f mxbrod.f
		$(peigs_CPP) < mxbrod.cpp > mxbrod.f
		$(peigs_FC) -c mxbrod.f 

mxcomb.o :	mxcomb.cpp
	       	rm -f mxcomb.f
		$(peigs_CPP) < mxcomb.cpp > mxcomb.f
		$(peigs_FC) -g -c mxcomb.f 

clock.o :	clock.c
		$(peigs_CC) -c clock.c 

mxcombv1.o :	mxcombv1.cpp
	       	rm -f mxcombv1.f
		$(peigs_CPP) < mxcombv1.cpp > mxcombv1.f
		$(peigs_FC) -c mxcombv1.f 

mxcombv.o :	mxcombv.cpp
	       	rm -f mxcombv.f
		$(peigs_CPP) < mxcombv.cpp > mxcombv.f
		$(peigs_FC) -g  -c mxcombv.f 

sum.o :		sum.cpp
		rm -f sum.f
		$(peigs_CPP) < sum.cpp > sum.f
		$(peigs_FC) -c sum.f 

.c.o:;	$(peigs_CC) $(peigs_COPT) -c $<

.f.o:;	$(peigs_FC) $(peigs_FOPT) -c $<

