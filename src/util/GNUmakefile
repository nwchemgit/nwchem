# $Id: GNUmakefile,v 1.14 1994-09-21 23:44:32 d3g681 Exp $

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
              util_hostname.o util_date.o \
	      nxtseg.o \
              ga_iter_orthog.o ga_iter_lsolve.o ga_iter_diag.o \
              ga_maxelt.o ga_pcg_minimize.o line_search.o \
              ga_orthog_vec.o ga_ran_fill.o ga_mix.o \
              ga_iter_project.o ga_screen.o \
              $(BLAS) $(EXTRA_OBJ) 
    LIBRARY = libutil.a
    HEADERS = util.fh itri.fh msgids.fh

    USES_BLAS = util.fh ga_iter_lsolve.F ga_maxelt.F ga_mix.F ga_iter_diag.F

include ../config/makelib.h


testsolve:      testsolve.o $(LIBDIR)/$(LIBRARY)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

ga_iter_project.o:	cproject.fh
