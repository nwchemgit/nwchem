include ../DEFS

MXOBJ	=	mxbrod.o \
		mxcombv1.o \
		mxsubs.o \
		sum.o

all:		$(MXOBJ)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:		
		-rm -f *.o
		-rm -f *.f

realclean:
		clean

mxsubs.o : 	mxsubs.cpp
	       	rm -f mxsubs.f
		$(CPP) -I$(HDIR) < mxsubs.cpp > mxsubs.f
		$(FC) -c -o mxsubs.o mxsubs.f

mxbrod.o :	mxbrod.cpp
	       	rm -f mxbrod.f
		$(CPP)  -I$(HDIR) < mxbrod.cpp > mxbrod.f
		$(FC) -c mxbrod.f -o mxbrod.o

mxcomb.o :	mxcomb.cpp
	       	rm -f mxcomb.f
		$(CPP) -I$(HDIR) < mxcomb.cpp > mxcomb.f
		$(FC) -g -c mxcomb.f  -o mxcomb.o

clock.o :	clock.c
		$(CC) $(OPT) -c clock.c  -o clock.o

mxcombv1.o :	mxcombv1.cpp
	       	rm -f mxcombv1.f
		$(CPP) -I$(HDIR) < mxcombv1.cpp > mxcombv1.f
		$(FC) -c mxcombv1.f -o mxcombv1.o

mxcombv.o :	mxcombv.cpp
	       	rm -f mxcombv.f
		$(CPP) -I$(HDIR) < mxcombv.cpp > mxcombv.f
		$(FC) -g  -c mxcombv.f -o mxcombv.o

sum.o :		sum.cpp
		rm -f sum.f
		$(CPP) -I$(HDIR) < sum.cpp > sum.f
		$(FC) -c sum.f  -o sum.o
