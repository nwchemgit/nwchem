# $Id: GNUmakefile,v 1.113 2000-05-11 19:44:23 edo Exp $

include ../config/makefile.h

   LIB_DEFINES = -DSCRATCH_DEF_DIR=$(SCRATCH_DEF_DIR) \
                 -DPERM_DEF_DIR=$(PERM_DEF_DIR)
ifdef JOBTIME_PATH
   LIB_DEFINES += -DJOBTIMEPATH=\"$(JOBTIME_PATH)\"
else
   LIB_DEFINES += -DJOBTIMEPATH=\"$(BINDIR)\"
endif

ifeq ($(TARGET),SOLARIS)
  EXTRA_OBJ = ieeetrap.o
#  SPEECH = YES
  LIB_DEFINES += -DHAVE_FORTRAN_INTEGER2 -DHAVE_FORTRAN_BYTE
endif

ifeq ($(TARGET),SUN)
  EXTRA_OBJ = ieeetrap.o
#  SPEECH = YES
  LIB_DEFINES += -DHAVE_FORTRAN_INTEGER2 -DHAVE_FORTRAN_BYTE
endif

ifeq ($(TARGET),SP)
  EXTRA_OBJ = ibm_cputime.o
  LIB_DEFINES += -DHAVE_FORTRAN_INTEGER2 -DHAVE_FORTRAN_BYTE
  LIB_ALSO_BUILD = $(BINDIR)/jobtime $(BINDIR)/jobtime.pl
endif

ifeq ($(TARGET),LAPI)
  EXTRA_OBJ = ibm_cputime.o
  LIB_DEFINES += -DHAVE_FORTRAN_INTEGER2 -DHAVE_FORTRAN_BYTE
  LIB_ALSO_BUILD = $(BINDIR)/jobtime $(BINDIR)/jobtime.pl
endif

ifeq ($(TARGET),IBM)
  EXTRA_OBJ = ibm_cputime.o
  LIB_DEFINES += -DHAVE_FORTRAN_INTEGER2 -DHAVE_FORTRAN_BYTE
endif
ifeq ($(TARGET),IBM64)
  EXTRA_OBJ = ibm_cputime.o
endif
ifeq ($(TARGET),HPUX64)
  EXTRA_OBJ = util_erfc.o
endif

ifeq ($(TARGET),$(findstring $(TARGET),LINUX CYGNUS))
  EXTRA_OBJ = linux_cpu.o linux_shift.o linux_random.o
#
# If you want to turn traps on for LINUX uncomment this line. That is it.
#  EXTRA_OBJ += linux_gnu_trap.o
#
endif

ifeq ($(TARGET),PGLINUX)
  EXTRA_OBJ = linux_cpu.o linux_shift.o linux_random.o
endif

ifeq ($(TARGET),HPUX)
  EXTRA_OBJ = util_hp_extra.o
endif

ifeq ($(TARGET),KSR)
  EXTRA_OBJ = ksr_cputime.o
endif

ifeq ($(TARGET),DECOSF)
  EXTRA_OBJ = dec_fpe.o
endif

ifeq ($(TARGET),FUJITSU_VPP)
   LIB_DEFINES = -Wp,-DSCRATCH_DEF_DIR=$(SCRATCH_DEF_DIR) \
                 -Wp,-DPERM_DEF_DIR=$(PERM_DEF_DIR) \
                 -Wp,-DBINDIR=\"$(BINDIR)\" \
                 -Wp,-DHAVE_FORTRAN_INTEGER2,-DHAVE_FORTRAN_BYTE
endif

ifeq ($(TARGET),SGITFP)
       BLAS = dfill.o ifill.o mabyte_fill.o
else
       BLAS = dfill.o ifill.o rsg.o mabyte_fill.o
endif
   LIB_TARGETS = testsolve testecce

 OBJ_OPTIMIZE =  ga_matpow.o util_pack.o dabssum.o dabsmax.o $(BLAS) \
                 ga_it_lsolve.o ga_it_orth.o ga_orthog.o idamin.o \
                 util_jacobi.o

  SUBDIRS = md5

        OBJ = output.o errquit.o ffflush.o \
              print_center.o util_flush.o \
              util_host.o util_date.o \
	      input_echo.o util_transpose.o \
              ga_iter_diag.o \
              ga_maxelt.o ga_pcg_min.o line_search.o \
              ga_orth_vec.o ga_ran_fill.o ga_mix.o ga_list.o\
              ga_it_proj.o ga_screen.o ga_get_diag.o \
	      fortchar.o seq_output.o ga_mat2col.o \
              util_ch_brd.o two_ind_trn.o \
	      util_pname.o sread.o swrite.o \
	      banner.o util_print.o util_version.o \
              mk_fit_xf.o \
              int_2c_ga.o \
              ga_local_mdot.o \
              util_cpusec.o util_wallsec.o \
	      gather.o scatter.o \
	      ga_trace_dg.o \
              icopy.o lcopy.o util_legal.o \
              util_file_name.o \
              util_speak.o util_rtdb_speak.o util_file_copy.o \
              util_file_unlink.o util_system.o util_sleep.o \
              util_rtdb_state.o ecce_print.o \
              util_random.o util_job.o util_getenv.o util_getarg.o\
              util_nwchemrc.o util_md.o util_md_sockets.o\
	      dgewr.o atoi.o indint.o util_wall_remain.o \
		ga_normf.o corr_mk_ref.o \
              nw_inp_from_file.o \
              bgj.o movecs_ecce.o\
              get_density.o moeig_read.o\
              util_debug.o \
              $(EXTRA_OBJ)

ifeq ($(TARGET),SGITFP)
	OBJ += rsg.o
endif
#nxtval_ga.o 

    LIBRARY = libutil.a
    HEADERS = util.fh itri.fh msgids.fh  numerical_constants.fh stdio.fh \
              printlevels.fh bitops.fh bitops_decls.fh bitops_funcs.fh \
              bgj.fh bgj_common.fh nwc_const.fh

    USES_BLAS = util.fh ga_it_lsolve.F ga_maxelt.F ga_mix.F ga_iter_diag.F \
                ga_orthog.F dabsmax.F ga_normf.F corr_mk_ref.F

ifdef SPEECH
   LIB_DEFINES += -DSPEECH
   LIB_TARGETS += speechserver speechclient
   OBJ += udp.o
endif

ifeq ($(TARGET),SGI)
   OBJ += sgi_flush6.o
endif
ifeq ($(TARGET),SGI_N32)
   OBJ += sgi_flush6.o
endif
ifeq ($(TARGET),SGITFP)
   OBJ += sgi_flush6.o
endif

include ../config/makelib.h

util_version.F:	util_v_stub.F
	cp $^ $@

version:	largeversion
	@ls -l util_version.F
largeversion:	../config/nwchem_config.h
	catsrc nwchem.F config $(NWSUBDIRS) | \
		awk -f ids.awk > util_version.F
smallversion:	../config/nwchem_config.h
	@echo " "
	@echo " perl must be in your path"
	@echo " "
	catsrc nwchem.F config $(NWSUBDIRS) | \
		awk -f ids.awk > util_version.F
	smallversion.pl

speechserver:	speechserver.c udp.c
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

speechclient:	speechclient.c udp.c
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

testsolve:      testsolve.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

testecce:      testecce.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)


ifeq ($(TARGET),SP)
$(BINDIR)/jobtime:	jobtime
	/bin/cp $^ $@

$(BINDIR)/jobtime.pl:	jobtime.pl
	/bin/cp $^ $@
endif
ifeq ($(TARGET),LAPI)
$(BINDIR)/jobtime:	jobtime
	/bin/cp $^ $@

$(BINDIR)/jobtime.pl:	jobtime.pl
	/bin/cp $^ $@
endif

