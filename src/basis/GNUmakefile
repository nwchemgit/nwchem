#$Id: GNUmakefile,v 1.41 2000-06-08 21:53:54 bert Exp $
include ../config/makefile.h

       LIBRARY = libutil.a
       HEADERS = bas.fh basP.fh bas_staticP.fh basdeclsP.fh geobasmapP.fh  \
                 bas_exndcf.fh bas_exndcf_dec.fh bas_exndcf_sfn.fh \
                 bas_ibs.fh bas_ibs_dec.fh bas_ibs_sfn.fh 
           OBJ = bas_input.o bas_contrib.o bas_checksum.o basis.o basisP.o bas_blas.o bas_blasP.o bas_vec_info.o

ifeq ($(TARGET),CRAY-T3D)
   OBJ_OPTIMIZE = basis.o 
   OBJ =  basisP.o bas_input.o bas_contrib.o bas_checksum.o bas_blas.o bas_blasP.o bas_vec_info.o
endif
ifeq ($(TARGET),CRAY-T3E)
   OBJ_OPTIMIZE = basis.o 
   OBJ =  basisP.o bas_input.o bas_contrib.o bas_checksum.o bas_blas.o bas_blasP.o bas_vec_info.o
endif

   LIB_TARGETS = testbasis testbasis.o libcheck

ifeq ($(TARGET),DELTA)
   LIB_DEFINES = -DBASIS_LIBRARY="'library'"
else
ifdef NWCHEM_BASIS_LIBRARY_PATH
   LIB_DEFINES = -DBASIS_LIBRARY="'$(NWCHEM_BASIS_LIBRARY_PATH)'"
else
   LIB_DEFINES = -DBASIS_LIBRARY="'$(SRCDIR)/basis/library'"
endif
endif

ifeq ($(TARGET),FUJITSU_VPP)
ifdef NWCHEM_BASIS_LIBRARY_PATH
   LIB_DEFINES = -Wp,-DBASIS_LIBRARY=$(NWCHEM_BASIS_LIBRARY_PATH)
else
   LIB_DEFINES = -Wp,-DBASIS_LIBRARY="'$(SRCDIR)/basis/library'"
endif
endif

     USES_BLAS = bas_blas.F bas_blasP.F

     TESTLIBS = -lutil

include ../config/makelib.h

compare:	compare.o
	$(FC) $(FFLAGS) -o $@ compare.o  $(LDFLAGS) $(LIBS) $(TESTLIBS)
libcheck:	libcheck.o
	$(FC) $(FFLAGS) -o $@ libcheck.o  $(LDFLAGS) $(LIBS) $(TESTLIBS)
testbasis:	testbasis.o
	$(FC) $(FFLAGS) -o $@ testbasis.o  $(LDFLAGS) $(LIBS) $(TESTLIBS)

basP.fh:	bas_staticP.fh
	@touch basP.fh

