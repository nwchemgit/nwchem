include ../DEFS

LIBS 	=	$(PEIGSLIB) $(LAPACKLIB) $(BLASLIB) $(MATHLIB) $(COMMLIB)
#
# Do not use -lkmath on delta, it occasionally causes 
# erroneous results (e.g. with tri_test.c with robert_tri as input
# on many different processor counts, including (8,8).
#

all:		geneig tritest wilk

teslamch:	teslamch.o slamch.o dlamch.o lsame.o
		$(LINK) $(OPT) -o teslamch teslamch.o slamch.o dlamch.o lsame.o

geneig_la:	gendrv.o geneig_la.o tim_com.o
		$(LINK) $(OPT) -o geneig_la gendrv.o geneig_la.o tim_com.o $(LIBS)

tsterr:		gendrv.o tsterr.o
		$(LINK) $(OPT) -o tsterr gendrv.o tsterr.o $(LIBS)

geneig:		gendrv.o geneig2.o tim_com.o
		$(LINK) $(OPT) -o geneig gendrv.o geneig2.o tim_com.o $(LIBS)

wilk:		gendrv.o wilk2.o tim_com.o
		$(LINK) $(OPT) -o wilk gendrv.o wilk2.o tim_com.o $(LIBS)

tritest:	gendrv.o tri_test.o tim_com.o
		$(LINK) $(OPT) -o tritest gendrv.o tri_test.o tim_com.o $(LIBS)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:
		-rm -f *.o
		-rm -f geneig tritest wilk teslamch
		-rm -f events* core

realclean:
		clean


eig.o:		eig.f
		$(FC) -c eig.f

gendrv.o:	gendrv.F
		$(CPP) < gendrv.F > gendrv.f
		$(FC) $(OPT) -c $(FCNFLAGS) gendrv.f
		rm gendrv.f


.f.o: ;	 $(FC) $(OPT) -c $(FCNFLAGS) $<
.c.o: ;  $(CC) $(OPTC)  -c $(CFLAGS) -I$(HDIR) $<
