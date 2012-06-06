# $Id$

include ../config/makefile.h

ifdef JOBTIME_PATH
   LIB_DEFINES += -DJOBTIMEPATH=\"$(JOBTIME_PATH)\"
else
   LIB_DEFINES += -DJOBTIMEPATH=\"$(BINDIR)\"
endif

# The packing routines have been modified so that they work
# for all combinations of 32-bit, 64-bit, big-endian, little-endian

# The old HAVE_FORTRAN_BYTE/INTEGER2 macros have been replaced with 
# BYTE_SHORT_PACKING to indicate that packing should use integer*1/*2 
# rather than the shift/and/or loops which are used if this macro
# is not defined.  This is because all machines have the short
# integers and the packing routines should now work everywhere
# (though some machines may be faster doing the shifts).

ifeq ($(TARGET),SOLARIS64)
# This assumes f77 ... use -DINTEGER_1='integer*1' for f90
  EXTRA_OBJ = ieeetrap.o
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
endif
ifeq ($(TARGET),SOLARIS)
# This assumes f77 ... use -DINTEGER_1='integer*1' for f90
  EXTRA_OBJ = ieeetrap.o
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
endif

ifeq ($(TARGET),SUN)
  EXTRA_OBJ = ieeetrap.o
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
endif

ifeq ($(TARGET),LAPI)
  EXTRA_OBJ = ibm_cputime.o
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
  LIB_ALSO_BUILD = $(BINDIR)/jobtime $(BINDIR)/jobtime.pl
endif

ifeq ($(TARGET),IBM)
  EXTRA_OBJ = ibm_cputime.o
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
endif

ifeq ($(TARGET),IBM64)
  EXTRA_OBJ = ibm_cputime.o 
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
endif

ifeq ($(TARGET),LAPI64)
  EXTRA_OBJ = ibm_cputime.o 
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
  LIB_ALSO_BUILD = $(BINDIR)/jobtime $(BINDIR)/jobtime.pl
endif

ifeq ($(TARGET),HPUX64)
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
  EXTRA_OBJ = ibm_cputime.o 
endif
ifeq ($(TARGET),MACX)
#  LIB_DEFINES += -DBYTE_SHORT_PACKING  
  EXTRA_OBJ += macx_trapfpe.o
endif

ifeq ($(TARGET),$(findstring $(TARGET),BGL BGP LINUX CYGWIN CYGNUS INTERIX))
  EXTRA_OBJ = linux_cpu.o linux_shift.o linux_random.o  
  ifeq ($(LINUXCPU),x86)
    ifneq ($(TARGET),INTERIX)
      EXTRA_OBJ += linux_setfpucw.o
    endif
  endif
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
  ifeq ($(FC),ifc)
    EXTRA_OBJ += erfc.o
  endif
#
# If you want to turn traps on for LINUX uncomment this line. That is it.
#  EXTRA_OBJ += linux_gnu_trap.o
#
endif


ifeq ($(TARGET),HPUX)
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
#  EXTRA_OBJ = util_hp_extra.o
endif

ifeq ($(TARGET),DECOSF)
  LIB_DEFINES += -DBYTE_SHORT_PACKING  -D_POSIX_PII_SOCKET
  EXTRA_OBJ = dec_fpe.o
endif
ifeq ($(NWCHEM_TARGET),$(findstring $(NWCHEM_TARGET),LINUX64 CATAMOUNT))
ifeq ($(_CPU),alpha)
ifeq ($(FC),fort)
  EXTRA_OBJ += dec_fpe.o
endif
endif
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
ifeq ($(_CPU),ia64)
  EXTRA_OBJ = win32_cpu.o linux_shift.o linux_random.o  linux_setfpucw.o erfc.o
endif
ifeq ($(_CPU),x86_64)
  EXTRA_OBJ = linux_cpu.o linux_shift.o linux_random.o   erfc.o  linux_setfpucw.o
  ifeq ($(FC),pathf90)
    EXTRA_OBJ += util_cflush.o
  endif
ifeq ($(FC),ftn)
  EXTRA_OBJ += util_allocga.o
endif
endif
ifeq ($(_CPU),ppc64)
  EXTRA_OBJ = linux_cpu.o linux_random.o # linux_setfpucw.o
endif
endif                                                                           
ifeq ($(NWCHEM_TARGET),CATAMOUNT)
   LIB_DEFINES += -I/opt/xt-catamount/default/catamount/linux/include/ -I/opt/xt-catamount/default/catamount/linux/include/amd64/
  EXTRA_OBJ += util_allocga.o xt3_dclock.o  util_fullbuff.o util_sigprof.o
endif
ifeq ($(ARMCI_NETWORK),PORTALS)
   LIB_DEFINES += -DGANXTVAL
endif
ifeq ($(TARGET),FUJITSU_VPP)
   LIB_DEFINES = \
                 -Wp,-DBINDIR=\"$(BINDIR)\" \
                 -Wp,-DBYTE_SHORT_PACKING \
                 -Wp,-DINTEGER_1=byte
endif

       BLAS = dfill.o ifill.o  mabyte_fill.o 

   LIB_TARGETS = testsolve testecce

 OBJ_OPTIMIZE =  ga_matpow.o util_pack.o dabssum.o dabsmax.o $(BLAS) \
                 ga_it_lsolve.o ga_it_orth.o ga_orthog.o idamin.o \
                 util_jacobi.o stpr_sjacobi.o util_memcpy.o ga_accback.o \
                 ga_asymmetr.o util_gnxtval.o nxtask.o util_mirror.o  \
                 util_sgroup.o

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
	      banner.o util_print.o \
              util_version.o util_nwchem_version.o util_ga_version.o \
              util_nwchem_paper.o\
              mk_fit_xf.o \
              int_2c_ga.o \
              ga_local_mdot.o \
              util_cpusec.o util_wallsec.o \
	      gather.o scatter.o \
	      ga_trace_dg.o \
              lcopy.o util_legal.o \
              util_file_name.o \
              util_io_unit.o \
              util_speak.o util_rtdb_speak.o util_file_copy.o \
              util_file_unlink.o util_system.o util_sleep.o \
              util_rtdb_state.o ecce_print.o \
              util_random.o util_job.o util_getenv.o util_getarg.o\
              util_nwchemrc.o util_md.o util_md_c.o util_md_sockets.o\
	      dgewr.o atoi.o indint.o util_wall_remain.o \
		ga_normf.o corr_mk_ref.o \
              nw_inp_from_file.o \
              bgj.o movecs_ecce.o\
              get_density.o moeig_read.o\
              util_file_copy.o \
              util_debug.o \
              util_erf.o  \
              util_length_scale.o  \
              ga_it2.o \
              ma_print.o freeze_input.o ga_extra.o \
              util_test.o util_ga_test.o util.o util_patch_test.o \
              util_ndim_test.o util_perf_test.o util_test_lu.o util_test_eig.o \
              util_dra_test.o util_eaf_test.o util_sf_test.o \
              ga_lkain_2cpl3.o util_io.o util_xyz.o util_ma.o util_mpinap.o \
              $(EXTRA_OBJ)
ifndef USE_MLIB
              OBJ_OPTIMIZE += icopy.o dsum.o dgefa.o
endif
#nxtval_ga.o 

    LIBRARY = libnwcutil.a
    HEADERS = util.fh itri.fh msgids.fh  numerical_constants.fh stdio.fh \
              printlevels.fh bitops.fh bitops_decls.fh bitops_funcs.fh \
              bgj.fh bgj_common.fh nwc_const.fh errquit.fh util_sgroup.fh \
	      util_params.fh

    USES_BLAS = util.fh ga_it_lsolve.F ga_maxelt.F ga_mix.F ga_iter_diag.F \
                ga_orthog.F dabsmax.F ga_normf.F corr_mk_ref.F ga_it2.F util_file_name.F dgefa.f util_patch_test.F stpr_sjacobi.F

ifdef SPEECH
   LIB_DEFINES += -DSPEECH
   LIB_TARGETS += speechserver speechclient
   OBJ += udp.o
endif

ifeq ($(TARGET),INTERIX)
DEFINES += -DNOSOCKLEN
endif
ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif
ifeq ($(TARGET),LINUX64)
DEFINES += -DNEED_LOC
endif
ifeq ($(TARGET),BGL)
DEFINES += -DNEED_LOC 
endif
ifeq ($(TARGET),BGP)
DEFINES += -DNEED_LOC
endif

ifeq ($(TARGET),SGI)
   OBJ += sgi_flush6.o
DEFINES += -DNOSOCKLEN
endif
ifeq ($(TARGET),SGI_N32)
   OBJ += sgi_flush6.o
DEFINES += -DNOSOCKLEN
endif
ifeq ($(TARGET),SGITFP)
   OBJ += sgi_flush6.o
DEFINES += -DNOSOCKLEN
endif
ifdef SLURM
  ifdef SLURMOPT
    DEFINES += -DSLURM
  endif
endif
ifdef LSF
  WHICHBJOBS  += $(shell which bjobs)
  DEFINES += -DLSF
  LIB_ALSO_BUILD = $(BINDIR)/jobtime_lsf $(BINDIR)/jobtime_lsf.pl
endif
ifdef PBS
  WHICHQSTAT  += $(shell which qstat)
  DEFINES += -DPBS
  LIB_ALSO_BUILD = $(BINDIR)/jobtime_pbs $(BINDIR)/jobtime_pbs.pl
endif
ifdef SLURM
  DEFINES += -DSLURM
endif

ifdef USE_NOFSCHECK
  DEFINES += -DNOFSCHECK
endif

include ../config/makelib.h


util_version.F:	util_v_stub.F
	cp $^ $@

version:	largeversion
	@ls -l util_version.F
largeversion:	../config/nwchem_config.h
	./catsrc nwchem.F config $(NWSUBDIRS) | \
		awk -f ids.awk > util_version.F
smallversion:	../config/nwchem_config.h
	@echo " "
	@echo " perl must be in your path"
	@echo " "
	./catsrc nwchem.F config $(NWSUBDIRS) | \
		awk -f ids.awk > util_version.F
	smallversion.pl
shortversion:	../config/nwchem_config.h
	-$(RM) util_version.F

speechserver:	speechserver.c udp.c
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

speechclient:	speechclient.c udp.c
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

testsolve:      testsolve.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

testecce:      testecce.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)


USE_LL = LAPI LAPI64

ifeq ($(TARGET),$(findstring $(TARGET),$(USE_LL)))
 LIB_DEFINES += -DUSE_LL
$(BINDIR)/jobtime:	jobtime
	/bin/cp $^ $@

$(BINDIR)/jobtime.pl:	jobtime.pl
	/bin/cp $^ $@
endif

ifdef LSF
$(BINDIR)/jobtime_lsf:	jobtime_lsf
	/bin/cp $^ $@

jobtime_lsf.pl:	
	/bin/rm -f jobtime_lsf.pl
	cat toplsf.pl > jobtime_lsf.pl
	echo "$$"bjobs=\"$(WHICHBJOBS)\"\;  >> jobtime_lsf.pl
	cat botlsf.pl >> jobtime_lsf.pl

$(BINDIR)/jobtime_lsf.pl:	jobtime_lsf.pl
	/bin/cp $^ $@
endif
ifdef PBS
$(BINDIR)/jobtime_pbs:	jobtime_pbs
	/bin/cp $^ $@

jobtime_pbs.pl:	
	/bin/rm -f jobtime_pbs.pl
	cat toplsf.pl > jobtime_pbs.pl
	echo "$$"qstat=\"$(WHICHQSTAT)\"\;  >> jobtime_pbs.pl
	cat botpbs.pl >> jobtime_pbs.pl

$(BINDIR)/jobtime_pbs.pl:	jobtime_pbs.pl
	/bin/cp $^ $@
endif

clean:
	-$(RM) util_version.F
