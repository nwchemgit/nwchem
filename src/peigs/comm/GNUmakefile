#
# $Id: GNUmakefile,v 1.14 1999-07-28 00:39:02 d3e129 Exp $
#
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
.F.o:;	$(peigs_FC) $(peigs_FOPT) -c $<


