# $Id: GNUmakefile,v 1.22 1994-11-18 18:08:45 d3g681 Exp $

include ../config/makefile.h

ifeq ($(TARGET),SUN)
  EXTRA_OBJ = ieeetrap.o
endif

ifeq ($(TARGET),IBM)
  EXTRA_OBJ = ibm_cputime.o
endif

       BLAS = dfill.o ifill.o rsg.o

   LIB_TARGETS = testsolve

OBJ_OPTIMIZE= ga_orthog.o 
        OBJ = output.o errquit.o ffflush.o \
              print_center.o util_flush.o \
              util_hostname.o util_date.o \
	      nxtseg.o \
              ga_iter_orthog.o ga_iter_lsolve.o ga_iter_diag.o \
              ga_maxelt.o ga_pcg_minimize.o line_search.o \
              ga_orthog_vec.o ga_ran_fill.o ga_mix.o \
              ga_iter_project.o ga_screen.o ga_get_diag.o \
<<<<<<< makefile
	      fortchar.o sequential_output.o ga_mat2col.o $(BLAS) $(EXTRA_OBJ) 
=======
	      fortchar.o sequential_output.o ga_mat2col.o \
	      banner.o $(BLAS) $(EXTRA_OBJ) 
OBJ_OPTIMIZE= ga_orthog.o 
>>>>>>> 1.21
    LIBRARY = libutil.a
    HEADERS = util.fh itri.fh msgids.fh

    USES_BLAS = util.fh ga_iter_lsolve.F ga_maxelt.F ga_mix.F ga_iter_diag.F \
                ga_orthog.F

include ../config/makelib.h


testsolve:      testsolve.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

ga_iter_project.o:	cproject.fh
