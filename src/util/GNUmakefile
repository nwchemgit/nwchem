# $Id$

include ../config/makefile.h

ifdef JOBTIME_PATH
   LIB_DEFINES += -DJOBTIMEPATH=\"$(JOBTIME_PATH)\"
else
ifdef XLFMAC
   LIB_DEFINES += -DJOBTIMEPATH="$(BINDIR)"
else
   LIB_DEFINES += -DJOBTIMEPATH=\"$(BINDIR)\"
endif
endif

ifeq ($(TARGET),MACX)
  LIB_DEFINES += -DOLD_GA
endif

   MXLENPATH = 200
   LIB_DEFINES += -DNWCHEM_LONG_PATHS
   LONGNWTOP = $(shell if [ `echo $(SRCDIR) |wc -m` -gt $(MXLENPATH) ]; then echo "Y"; fi  )
ifeq ($(LONGNWTOP),Y)
errorlongpwd:
	@echo " "
	@echo "The directory name chosen for NWCHEM_TOP is longer than "
	@echo "the maximum allowed value of" $(MXLENPATH) "characters"
	@echo "      current NWCHEM_TOP/src="$(SRCDIR)" equal to" $(shell echo $(SRCDIR) |wc -m) "characters"
	@echo "please chose a directory with a shorter name"
	@echo "E.g., setenv NWCHEM_TOP /home/nwchem_user/nwchem"
	@exit 1
endif
   LIB_DEFINES += -DNWCHEM_SRCDIR="\"$(SRCDIR)\""



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
  LIB_DEFINES += -DNO_MPIF
  EXTRA_OBJ += macx_trapfpe.o
endif
ifeq ($(TARGET),MACX64)
#  LIB_DEFINES += -DBYTE_SHORT_PACKING  
  LIB_DEFINES += -DNO_MPIF
  EXTRA_OBJ += macx_trapfpe.o
endif

ifdef NO_MPIF
  LIB_DEFINES += -DNO_MPIF
endif

ifeq ($(TARGET),$(findstring $(TARGET),BGL BGP BGQ LINUX CYGWIN CYGWIN64 CYGNUS LINUX64 MACX MACX64))
  EXTRA_OBJ = linux_cpu.o linux_shift.o linux_random.o  erfc.o linux_printaff.o
  ifeq ($(LINUXCPU),x86)
    EXTRA_OBJ += linux_setfpucw.o
  endif
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
  ifeq ($(FC),ifc)
    EXTRA_OBJ += erfc.o
  endif
endif

ifeq ($(NWCHEM_TARGET),$(findstring $(NWCHEM_TARGET),LINUX64 CATAMOUNT))
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
  ifeq ($(_CPU),ia64)
    EXTRA_OBJ += win32_cpu.o linux_setfpucw.o 
  endif
  ifeq ($(GOTMINGW64),1)
    EXTRA_OBJ += util_mingw.o
  endif
  ifeq ($(_CPU),x86_64)
  ifneq ($(GOTMINGW64),1)
    EXTRA_OBJ  += linux_setfpucw.o
  endif
  endif
  ifeq ($(FC),pathf90)
    EXTRA_OBJ += util_cflush.o
  endif
endif
ifeq ($(GOTMINGW32),1)
    LIB_DEFINES += -DNO_MPIF
    EXTRA_OBJ += util_mingw.o
endif

ifeq ($(TARGET),HPUX)
  LIB_DEFINES += -DBYTE_SHORT_PACKING  
endif

ifeq ($(NWCHEM_TARGET),CATAMOUNT)
   LIB_DEFINES += -I/opt/xt-catamount/default/catamount/linux/include/ -I/opt/xt-catamount/default/catamount/linux/include/amd64/
  EXTRA_OBJ += util_allocga.o xt3_dclock.o  util_fullbuff.o util_sigprof.o
endif
ifeq ($(ARMCI_NETWORK),PORTALS)
   LIB_DEFINES += -DGANXTVAL
endif

ifneq ($(OLD_GA),)
  EXTRA_OBJ += ga_zgop.o
endif
#LIB_DEFINES += -I$(MPI_INCLUDE)
ifneq ($(MPI_INCLUDE),)
   LIB_INCLUDES += $(patsubst -I-I%,-I%,-I$(MPI_INCLUDE))
endif

ifdef USE_VTUNE
  EXTRA_OBJ += vtune_bind.o vtune.o
endif

ifeq ($(shell echo $(BLASOPT) |awk '/Accelerate/ {print "Y"; exit}'),Y)
  EXTRA_OBJ += veclib_threads.o
endif
_PDSYEVR = $(shell [ -e ${NWCHEM_TOP}/src/tools/install/lib/libga.a ] &&  nm ${NWCHEM_TOP}/src/tools/install/lib/libga.a|egrep T |egrep  ga_pdsyevr |awk '/ga_pdsyevr/ {print "Y"}')
ifeq ($(_PDSYEVR),Y) 
 LIB_DEFINES += -DUSE_PDSYEVR
endif
_USE_ELPA = $(shell ${GA_PATH}/bin/ga-config  --use_elpa| awk ' /1/ {print "Y"}')
ifeq ($(_USE_ELPA),Y)
 LIB_DEFINES += -DUSE_ELPA
endif

       BLAS = dfill.o ifill.o  mabyte_fill.o zfill.o

   LIB_TARGETS = testsolve testecce

 OBJ_OPTIMIZE =  ga_matpow.o util_pack.o dabssum.o dabsmax.o $(BLAS) \
                 ga_it_lsolve.o ga_it_orth.o ga_orthog.o idamin.o \
                 util_jacobi.o util_rsg.o stpr_sjacobi.o util_memcpy.o ga_accback.o \
                 ga_asymmetr.o util_gnxtval.o nxtask.o util_mirror.o  \
                 ga_copy_dz.o dgami.o util_eri_cholesky.o \
                 util_test_cholesky.o util_sillysort.o util_talker.o

  SUBDIRS = md5

        OBJ = output.o errquit.o ffflush.o \
              print_center.o util_flush.o util_sgroup.o \
              util_host.o util_date.o \
	      input_echo.o util_transpose.o \
              ga_iter_diag.o \
              ga_maxelt.o nga_idamax.o ga_pcg_min.o line_search.o \
              ga_orth_vec.o ga_ran_fill.o ga_mix.o ga_list.o\
              ga_it_proj.o ga_screen.o ga_get_diag.o \
	      fortchar.o seq_output.o ga_mat2col.o \
              util_ch_brd.o two_ind_trn.o \
	      util_pname.o sread.o swrite.o \
	      banner.o util_print.o util_module_avail.o \
              util_version.o util_nwchem_version.o util_ga_version.o \
              util_nwchem_paper.o util_nwchem_srcdir.o utilc_nwchem_srcdir.o \
              mk_fit_xf.o \
              int_2c_ga.o \
              ma_solve.o \
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
              util_cflush.o  \
              util_ascii_plot.o \
              util_length_scale.o  \
              ga_it2.o ga_lkain_ext.o\
              ma_print.o freeze_input.o ga_extra.o \
              util_test.o util_ga_test.o util.o util_patch_test.o \
              util_ndim_test.o util_perf_test.o util_test_lu.o util_test_eig.o \
              util_dra_test.o util_eaf_test.o util_sf_test.o \
              ga_lkain_2cpl3.o util_io.o util_xyz.o util_ma.o util_mpinap.o \
	      ga_lkain_2cpl3_ext.o  util_dgeev.o util_nwsort.o \
	      util_mygabcast.o util_getppn.o util_units.o util_align.o\
              util_intrsc2008.o util_fadvise.o  util_scalapack_info.o\
              ga_conjg.o nwmolden.o ao_1prdm_read.o ao_1prdm_write.o\
              util_strips.o util_poltensor.o  util_mygop.o util_isnan.o \
              util_blasthreads.o util_gpu_affinity.o util_scalapack.o \
              util_norm2.o util_dpofa.o util_dgedi.o \
              $(EXTRA_OBJ)
ifndef USE_MLIB
              OBJ_OPTIMIZE += icopy.o dsum.o dgefa.o zsum.o
endif
#nxtval_ga.o 

    LIBRARY = libnwcutil.a
    HEADERS = util.fh itri.fh msgids.fh  numerical_constants.fh stdio.fh \
              printlevels.fh bitops.fh bitops_decls.fh bitops_funcs.fh \
              bgj.fh bgj_common.fh nwc_const.fh errquit.fh util_sgroup.fh \
	      util_params.fh utilc_nwchem_srcdir.fh util_maxlength.fh util_intsize.fh
#if USE_VTUNE
    HEADERS += vtune.fh
#endif

    USES_BLAS = util.fh ga_it_lsolve.F ga_maxelt.F ga_mix.F ga_iter_diag.F \
                ga_orthog.F dabsmax.F ga_normf.F corr_mk_ref.F ga_it2.F ga_lkain_ext.F util_file_name.F dgefa.f util_patch_test.F stpr_sjacobi.F util_dgeev.F \
		util_test_cholesky.F dfill.F ga_lkain_2cpl3_ext.F ga_it2.F util_poltensor.F util_dpofa.F util_dgedi.F

ifdef SPEECH
   LIB_DEFINES += -DSPEECH
   LIB_TARGETS += speechserver speechclient
   OBJ += udp.o
endif

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif
ifeq ($(TARGET),LINUX64)
DEFINES += -DNEED_LOC
endif
ifeq ($(TARGET),$(findstring $(TARGET),BGL BGP BGQ))
DEFINES += -DNEED_LOC
LIB_DEFINES += -DNO_UTIL_TESTS
LIB_DEFINES += -I/bgsys/drivers/ppcfloor/firmware/include -I/bgsys/drivers/ppcfloor/spi/include/kernel
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


ifeq ($(ARMCI_NETWORK),DMAPP)
  DEFINES += -DDISABLE_GAMIRROR
endif
ifdef DISABLE_GAMIRROR
  DEFINES += -DDISABLE_GAMIRROR
endif

ifdef USE_NOIO
  DEFINES += -DNOFSCHECK
endif
ifdef USE_NOFSCHECK
  DEFINES += -DNOFSCHECK
endif
ifdef TCE_CUDA
  LIB_DEFINES += $(CUDA_INCLUDE)
  DEFINES += -DTCE_CUDA
  OBJ += util_cuda_support.o
endif

ifdef TCE_HIP
  LIB_DEFINES += $(HIP_INCLUDE)
  OBJ += util_cuda_support.o
endif

ifdef USE_OPENACC_TRPDRV
  FOPTIONS += -cuda -DUSE_CUDA_AFFINITY
endif

ifdef USE_FORTRAN2003
  LIB_DEFINES += -DUSE_FORTRAN2003
endif
ifdef USE_NONBLOCKING_BCAST
  LIB_DEFINES += -DUSE_IBCAST
endif

include ../config/makelib.h


util_scalapack_info.F: 
	./util_scalapack_info.bash
ifdef EXTERNAL_GA_PATH
util_ga_version.F: 
	./util_ga_version.bash
endif


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
	./shortsrc

speechserver:	speechserver.c udp.c
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

speechclient:	speechclient.c udp.c
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

testsolve:      testsolve.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

testecce:      testecce.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

test_intrsc2008: test_intrsc2008.o $(LIBRARY_PATH)
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
