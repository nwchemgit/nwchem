# $Id: GNUmakefile,v 1.49 1995-11-27 02:19:56 d3g681 Exp $

include ../config/makefile.h

ifeq ($(TARGET),SUN)
  EXTRA_OBJ = ieeetrap.o
endif

ifeq ($(TARGET),IBM)
  EXTRA_OBJ = ibm_cputime.o
endif

ifeq ($(TARGET),LINUX)
  EXTRA_OBJ = linux_cputime.o
endif

ifeq ($(TARGET),KSR)
  EXTRA_OBJ = ksr_cputime.o
endif

       BLAS = dfill.o ifill.o rsg.o

   LIB_TARGETS = testsolve

 OBJ_OPTIMIZE =  ga_matpow.o util_pack.o dabssum.o dabsmax.o $(BLAS)


        OBJ = output.o errquit.o ffflush.o \
              print_center.o util_flush.o \
              util_hostname.o util_date.o \
	      ga_orthog.o file_prefix.o \
              ga_iter_orthog.o ga_iter_lsolve.o ga_iter_diag.o \
              ga_maxelt.o ga_pcg_minimize.o line_search.o \
              ga_orthog_vec.o ga_ran_fill.o ga_mix.o \
              ga_iter_project.o ga_screen.o ga_get_diag.o \
	      fortchar.o sequential_output.o ga_mat2col.o \
              util_char_ga_brdcst.o two_index_transf.o \
	      util_pname.o sread.o swrite.o \
	      banner.o util_print.o util_version.o \
              mk_fit_xf.o \
              int_2c_ga.o \
              util_cpusec.o util_wallsec.o \
	      gather.o scatter.o \
	      ga_trace_diag.o \
              icopy.o \
              $(EXTRA_OBJ) 

    LIBRARY = libutil.a
    HEADERS = util.fh itri.fh msgids.fh  numerical_constants.fh stdio.fh \
              printlevels.fh bitops.fh bitops_decls.fh bitops_funcs.fh

    USES_BLAS = util.fh ga_iter_lsolve.F ga_maxelt.F ga_mix.F ga_iter_diag.F \
                ga_orthog.F


include ../config/makelib.h

util_version.F:	util_version_stub.F
	cp $^ $@

version:	../config/nwchem_config.h
	catsrc nwchem.F config $(NWSUBDIRS) | \
		awk -f ids.awk > util_version.F

testsolve:      testsolve.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)



# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(error.o):	 $(INCDIR)/errorP.h 
$(LIBRARY_PATH)(file_prefix.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/inp.fh ./cfileprefix.fh 
$(LIBRARY_PATH)(fortchar.o):	 $(INCDIR)/sndrcv.h 
$(LIBRARY_PATH)(ga_get_diag.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_iter_diag.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_iter_lsolve.o):	 ./util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_iter_orthog.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(ga_iter_project.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cproject.fh 
$(LIBRARY_PATH)(ga_mat2col.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_matpow.o):	 ./stdio.fh ./numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_maxelt.o):	 ./msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_mix.o):	 ./util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_orthog.o):	 ./msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_orthog_vec.o):	 $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_pcg_minimize.o):	 ./util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_ran_fill.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_screen.o):	 ./msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_trace_diag.o):	 ./msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(int_2c_ga.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(line_search.o):	 ./util.fh 
$(LIBRARY_PATH)(mk_fit_xf.o):	 ./stdio.fh ./numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(print_center.o):	 $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(sequential_output.o):	 $(INCDIR)/tcgmsg.fh ./msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(testsolve.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(two_index_transf.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(util_char_ga_brdcst.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(util_pack.o):	 $(INCDIR)/mafdecls.fh ./bitops.fh 
$(LIBRARY_PATH)(util_pname.o):	 $(INCDIR)/inp.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(util_print.o):	 $(INCDIR)/rtdb.fh ./printlevels.fh ./printcommon.fh $(INCDIR)/inp.fh 
