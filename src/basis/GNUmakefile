#$Id: GNUmakefile,v 1.32 1997-02-17 20:51:28 d3g681 Exp $
include ../config/makefile.h

       LIBRARY = libutil.a
       HEADERS = bas.fh basP.fh bas_staticP.fh basdeclsP.fh geobasmapP.fh  \
                 bas_exndcf.fh bas_exndcf_dec.fh bas_exndcf_sfn.fh \
                 bas_ibs.fh bas_ibs_dec.fh bas_ibs_sfn.fh 
#           OBJ = bas_input.o bas_contrib.o bas_checksum.o
#  OBJ_OPTIMIZE = basis.o basisP.o
           OBJ = bas_input.o bas_contrib.o bas_checksum.o basis.o basisP.o

ifeq ($(TARGET),CRAY-T3D)
   OBJ_OPTIMIZE = basis.o 
   OBJ =  basisP.o bas_input.o bas_contrib.o bas_checksum.o
endif

   LIB_TARGETS = testbasis testbasis.o

ifeq ($(TARGET),DELTA)
   LIB_DEFINES = -DBASIS_LIBRARY="'library'"
else
   LIB_DEFINES = -DBASIS_LIBRARY="'$(SRCDIR)/basis/library'"
endif

     USES_BLAS = basis.F basisP.F

     TESTLIBS = -lutil

include ../config/makelib.h

testbasis:	testbasis.o
	$(FC) $(FFLAGS) -o $@ testbasis.o  $(LDFLAGS) $(LIBS) $(TESTLIBS)

basP.fh:	bas_staticP.fh
	@touch basP.fh

