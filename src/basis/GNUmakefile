#$Id: GNUmakefile,v 1.44 2001-12-14 19:24:56 bert Exp $
TARGET = $(NWCHEM_TARGET)
include ../tools/armci/config/makecoms.h
include ../tools/armci/config/makemp.h
include ../config/makefile.h

       LIBRARY = libnwcutil.a
       HEADERS = bas.fh basP.fh bas_staticP.fh basdeclsP.fh geobasmapP.fh  \
                 bas_exndcf.fh bas_exndcf_dec.fh bas_exndcf_sfn.fh \
                 bas_ibs.fh bas_ibs_dec.fh bas_ibs_sfn.fh ../NWints/rel/rel_nwc.fh
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
   LIB_DEFINES = -DBASIS_LIBRARY="'libraries/'"
else
ifdef NWCHEM_BASIS_LIBRARY_PATH
   LIB_DEFINES = -DBASIS_LIBRARY="'$(NWCHEM_BASIS_LIBRARY_PATH)'"
else
   LIB_DEFINES = -DBASIS_LIBRARY="'$(SRCDIR)/basis/libraries/'"
endif
endif

ifeq ($(TARGET),FUJITSU_VPP)
ifdef NWCHEM_BASIS_LIBRARY_PATH
   LIB_DEFINES = -Wp,-DBASIS_LIBRARY=$(NWCHEM_BASIS_LIBRARY_PATH)
else
   LIB_DEFINES = -Wp,-DBASIS_LIBRARY="'$(SRCDIR)/basis/libraries/'"
endif
endif

     USES_BLAS = bas_blas.F bas_blasP.F

include ../config/makelib.h

compare:	compare.o
	$(LINK.f) -o $@ compare.o  $(LIBS) 
libcheck:	libcheck.o
	$(LINK.f) -o $@ libcheck.o $(LIBS) 
testbasis:	testbasis.o
	$(LINK.f) -o $@ testbasis.o  $(LIBS)

basP.fh:	bas_staticP.fh
	@touch basP.fh

