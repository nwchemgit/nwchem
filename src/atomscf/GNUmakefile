# $Id: GNUmakefile,v 1.6 1995-03-31 01:32:19 d3g681 Exp $

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
