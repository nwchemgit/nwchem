#
# $Id$
#

include ../../config/makefile.h
include ../DEFS

peigs_LIBS = -lpeigs -llapack $(peigs_PEIGSLIB) $(peigs_LAPACKLIB) $(peigs_BLASLIB) $(peigs_MATHLIB) $(peigs_COMMLIB) -lblas
#
# Do not use -lkmath on delta, it occasionally causes 
# erroneous results (e.g. with tri_test.c with robert_tri as input
# on many different processor counts, including (8,8).
#

OBJECT = teslamch.o slamch.o dlamch.o lsame.o xerbla.o tim_com.o \
	tsterr.o gendrv.o geneig2.o wilk2.o tritest.o geneig_la.o

all:		geneig tritest wilk

teslamch:	teslamch.o slamch.o dlamch.o lsame.o xerbla.o
		$(peigs_LINK) $(peigs_FOPT) -o teslamch teslamch.o slamch.o dlamch.o lsame.o  xerbla.o

geneig_la:	gendrv.o geneig_la.o tim_com.o xerbla.o
		$(peigs_LINK) $(peigs_FOPT) -o geneig_la \
		gendrv.o geneig_la.o tim_com.o  xerbla.o $(peigs_LIBS)

tsterr:		gendrv.o tsterr.o xerbla.o
		$(peigs_LINK) $(peigs_FOPT) -o tsterr gendrv.o  xerbla.o \
		tsterr.o $(peigs_LIBS)

geneig:		gendrv.o geneig2.o tim_com.o xerbla.o
		$(peigs_LINK) $(peigs_FOPT) -o geneig gendrv.o  xerbla.o \
		geneig2.o tim_com.o $(peigs_LIBS)

wilk:		gendrv.o wilk2.o tim_com.o xerbla.o
		$(peigs_LINK) $(peigs_FOPT) -o wilk gendrv.o wilk2.o  xerbla.o \
		 tim_com.o  $(peigs_LIBS)

tritest:	gendrv.o tri_test.o tim_com.o xerbla.o
		$(peigs_LINK) $(peigs_FOPT) -o tritest gendrv.o  xerbla.o tri_test.o tim_com.o $(peigs_LIBS)

#.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:
		-rm -f *.o
		-rm -f geneig tritest wilk teslamch
		-rm -f events* core

realclean:
		clean


eig.o:		eig.f
		$(peigs_FC) -c eig.f

xerbla.o:	xerbla.f
		$(peigs_FC) -c xerbla.f

gendrv.o:	gendrv.F
		$(peigs_CPP) < gendrv.F > gendrvf.f
		$(peigs_FC) $(peigs_FOPT) -c gendrvf.f
		mv gendrvf.o gendrv.o
		rm gendrvf.f


.f.o: ;	 $(peigs_FC) -g $(peigs_FOPT) -c $<
.c.o: ;  $(peigs_CC) -g $(peigs_COPT)  -c $<


