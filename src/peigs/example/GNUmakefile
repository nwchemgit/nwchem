include ../../config/makefile.h
include ../DEFS

peigs_LIBS = $(peigs_PEIGSLIB) $(peigs_LAPACKLIB) $(peigs_BLASLIB) $(peigs_MATHLIB) $(peigs_COMMLIB)
#
# Do not use -lkmath on delta, it occasionally causes 
# erroneous results (e.g. with tri_test.c with robert_tri as input
# on many different processor counts, including (8,8).
#

all:		geneig tritest wilk

teslamch:	teslamch.o slamch.o dlamch.o lsame.o
		$(peigs_LINK) $(peigs_FOPT) -o teslamch teslamch.o slamch.o dlamch.o lsame.o

geneig_la:	gendrv.o geneig_la.o tim_com.o
		$(peigs_LINK) $(peigs_FOPT) -o geneig_la \
		gendrv.o geneig_la.o tim_com.o $(peigs_LIBS)

tsterr:		gendrv.o tsterr.o
		$(peigs_LINK) $(peigs_FOPT) -o tsterr gendrv.o \
		tsterr.o $(peigs_LIBS)

geneig:		gendrv.o geneig2.o tim_com.o
		$(peigs_LINK) $(peigs_FOPT) -o geneig gendrv.o \
		geneig2.o tim_com.o $(peigs_LIBS) \
		-lpeigs -llapack -lblas -ltcgmsg -lsci

wilk:		gendrv.o wilk2.o tim_com.o
		$(peigs_LINK) $(peigs_FOPT) -o wilk gendrv.o wilk2.o \
		 tim_com.o  \
		-lpeigs -llapack -lblas -ltcgmsg -lsci

tritest:	gendrv.o tri_test.o tim_com.o
		$(peigs_LINK) $(peigs_FOPT) -o tritest gendrv.o tri_test.o tim_com.o $(peigs_LIBS)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:
		-rm -f *.o
		-rm -f geneig tritest wilk teslamch
		-rm -f events* core

realclean:
		clean


eig.o:		eig.f
		$(peigs_FC) -c eig.f

gendrv.o:	gendrv.F
		$(peigs_CPP) < gendrv.F > gendrv.f
		$(peigs_FC) $(peigs_FOPT) -c gendrv.f
		rm gendrv.f


.f.o: ;	 $(peigs_FC) $(peigs_FOPT) -c $<
.c.o: ;  $(peigs_CC) $(peigs_COPT)  -c $<
