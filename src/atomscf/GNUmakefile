# $Id: GNUmakefile,v 1.5 1995-03-16 23:50:24 rg240 Exp $

          OBJ = orderd.o 
#
# KSR optimize barfs on orderd.F !?
#

 OBJ_OPTIMIZE = atcond.o \
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
