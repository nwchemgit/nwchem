# $Id$

# well actually in large basis sets the atomic guess can take a long time

 OBJ_OPTIMIZE = orderd.o \
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
		densid_z.o \
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
		printbasis.o \
                atscf_ecp.o \
                atscf_atscf.o

  USES_BLAS = atomd.F atscf_atscf.F atscf_ecp.F

  LIBRARY = libguess.a

include ../config/makefile.h
include ../config/makelib.h
