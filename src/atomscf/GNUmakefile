# $Id: GNUmakefile,v 1.24 2000-08-08 23:49:31 bjohnson Exp $

# optimization breaks on the SGITFP ... there's no point
# optimizing this stuff anyway.

# well actually in large basis sets the atomic guess can take
# a long time ... optimize except on the SGI for now.

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

# BEGIN_WIN32_IGNORE
ifeq ($(TARGET),SGITFP)
  OBJ := $(OBJ_OPTIMIZE)
  OBJ_OPTIMIZE :=
endif
# END_WIN32_IGNORE

  USES_BLAS = atomd.F atscf_atscf.F atscf_ecp.F

  LIBRARY = libguess.a

include ../config/makefile.h
include ../config/makelib.h
