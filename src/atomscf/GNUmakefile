# $Id: GNUmakefile,v 1.7 1995-10-17 08:52:41 d3g681 Exp $

# optimization breaks on the SGITFP ... there's no point
# optimizing this stuff anyway.

          OBJ = orderd.o \
                atcond.o \
                denat.o \
		hamild.o \
		outpud.o \
		teigd.o \
		atomd.o \
		denmad.o \
		jacod.o \
		pdfded.o \
		tracd.o \
		cmergd.o \
		densid.o \
		oeigd.o \
		shalfd.o \
		trafsd.o \
		creded.o \
		square.o \
		tramad.o \
		datoms.o \
		guess_mem.o  \
		starcd.o \
		writel.o \
		printbasis.o

  USES_BLAS = atomd.F


  LIBRARY = libguess.a


include ../config/makefile.h
include ../config/makelib.h


$(OBJ):	cguess.fh
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
atomd.o :	 $(INCDIR)/inp.fh cguess.fh 
creded.o :	 $(INCDIR)/global.fh cguess.fh 
datoms.o :	 $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cguess.fh $(INCDIR)/bas.fh 
denat.o :	 $(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cguess.fh $(INCDIR)/bas.fh 
guess_mem.o :	 cguess.fh 
hamild.o :	 cguess.fh 
oeigd.o :	 cguess.fh 
osatod.o :	 $(INCDIR)/geom.fh cguess.fh $(INCDIR)/bas.fh 
outpud.o :	 cguess.fh 
printbasis.o :	 cguess.fh 
teigd.o :	 cguess.fh 
