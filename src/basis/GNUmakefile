#$Id$

#######################################################
# This is only necessary to detect Cray T3D and T3E
# for which the only difference is to optimize basis.o,
# which is pointless since these machines are dead.
#
#TARGET = $(NWCHEM_TARGET)
#  export ARMCI_DIR_NAME := armci
#ifeq ($(ARMCI_NETWORK),PORTALS)
#  export ARMCI_DIR_NAME := armci-portals
#endif
#include ../tools/$(ARMCI_DIR_NAME)/config/makecoms.h
#include ../tools/$(ARMCI_DIR_NAME)/config/makemp.h
########################################################

include ../config/makefile.h

       LIBRARY = libnwcutil.a
       HEADERS = bas.fh basP.fh bas_staticP.fh basdeclsP.fh geobasmapP.fh  \
                 bas_exndcf.fh bas_exndcf_dec.fh bas_exndcf_sfn.fh \
                 bas_ibs.fh bas_ibs_dec.fh bas_ibs_sfn.fh ../NWints/rel/rel_nwc.fh
           OBJ  = basis.o bas_input.o bas_contrib.o bas_checksum.o basisP.o bas_blas.o bas_blasP.o bas_vec_info.o

########################################################
#ifeq ($(TARGET),CRAY-T3D)
#   OBJ_OPTIMIZE = basis.o 
#           OBJ  = bas_input.o bas_contrib.o bas_checksum.o basisP.o bas_blas.o bas_blasP.o bas_vec_info.o
#endif
#ifeq ($(TARGET),CRAY-T3E)
#   OBJ_OPTIMIZE = basis.o 
#           OBJ  = bas_input.o bas_contrib.o bas_checksum.o basisP.o bas_blas.o bas_blasP.o bas_vec_info.o
#endif
########################################################

   LIB_TARGETS = testbasis testbasis.o libcheck
   LONGNWTOP = $(shell if [ `echo $(SRCDIR) |wc -m` -gt 64 ]; then echo "Y"; fi  )
ifeq ($(LONGNWTOP),Y)
errorlongpwd:
	@echo " "
	@echo "The directory name chosen for NWCHEM_TOP is longer than "
	@echo "the maximum allowed value of 65 characters"
	@echo "      current NWCHEM_TOP="$(SRCDIR)" equal to" $(shell echo $(BSRCDIR) |wc -m) "characters"
	@echo "please chose a directory with a shorter name"
	@echo "E.g., setenv NWCHEM_TOP /home/nwchem_user/nwchem"
	@exit 1
endif
   LIB_DEFINES = -DBASIS_LIBRARY="'$(SRCDIR)'"

ifeq ($(TARGET),FUJITSU_VPP)
   LIB_DEFINES = -Wp,-DBASIS_LIBRARY="'$(SRCDIR)'"
endif

     USES_BLAS = bas_blas.F bas_blasP.F basis.F

include ../config/makelib.h

compare:	compare.o
	$(LINK.f) -o $@ compare.o  $(LIBS) 
libcheck:	libcheck.o
	$(LINK.f) -o $@ libcheck.o $(LIBS) 
testbasis:	testbasis.o
	$(LINK.f) -o $@ testbasis.o  $(LIBS)

basP.fh:	bas_staticP.fh
	@touch basP.fh

