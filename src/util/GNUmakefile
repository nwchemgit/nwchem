# $Id: GNUmakefile,v 1.7 1994-08-17 06:23:13 d3e129 Exp $

include ../config/makefile.h

ifeq ($(TARGET),SUN)
  EXTRA_OBJ = ieeetrap.o
else
  EXTRA_OBJ = 
endif

       BLAS = dfill.o ifill.o rsg.o

   LIB_TARGETS = testsolve

        OBJ = output.o errquit.o ffflush.o \
              print_center.o util_flush.o \
              ga_iter_orthog.o ga_iter_lsolve.o ga_iter_diag.o \
              ga_maxelt.o ga_pcg_minimize.o line_search.o $(BLAS) $(EXTRA_OBJ) 
    LIBRARY = libutil.a
    HEADERS = util.fh itri.fh

    USES_BLAS = util.fh ga_iter_lsolve.F ga_maxelt.F

include ../config/makelib.h


testsolve:      testsolve.o $(LIBRARY)
	$(FC) $(FLDOPT) -o $@ testsolve.o $(LIBS)

