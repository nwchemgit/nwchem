# DFTD4 interface for NWChem.
# Developed with assistance from OpenAI Codex; integration and numerical
# validation were reviewed by the author.
include ../../config/makefile.h

DFTD4_VERSION ?= 4.2.0
DFTD4_PREFIX = $(CURDIR)/install
DFTD4_LIBDIR = $(DFTD4_PREFIX)/lib
DFTD4_INCDIR = $(DFTD4_PREFIX)/include
DFTD4_MODULEDIR = $(CURDIR)/dftd4/build/include
DFTD4_BRIDGE = $(LIBDIR)/libnwc_dftd4.a
DFTD4_STAMP = $(DFTD4_PREFIX)/.dftd4-$(DFTD4_VERSION)-blas$(BLAS_SIZE)

all: $(DFTD4_BRIDGE)

$(DFTD4_STAMP): build_dftd4.sh
	DFTD4_BLAS_SIZE="$(BLAS_SIZE)" \
	DFTD4_BLAS_LIBRARIES="$(BLASOPT)" \
	DFTD4_LAPACK_LIBRARIES="$(LAPACK_LIB)" \
	./build_dftd4.sh $(DFTD4_VERSION) $(DFTD4_PREFIX)
	@touch $@

$(DFTD4_BRIDGE): nw_dftd4.c nw_dftd4_params.F90 $(DFTD4_STAMP)
	$(CC) $(COPTIONS) -I$(DFTD4_INCDIR) -c nw_dftd4.c -o nw_dftd4.o
	$(FC) $(FOPTIONS) -I$(DFTD4_MODULEDIR) -c nw_dftd4_params.F90 -o nw_dftd4_params.o
	$(RM) $@
	$(AR) $(ARFLAGS) $@ nw_dftd4.o nw_dftd4_params.o
	$(RANLIB) $@

clean:
	$(RM) nw_dftd4.o nw_dftd4_params.o $(DFTD4_BRIDGE)

realclean: clean
	$(RM) -r dftd4-* dftd4 install

64_to_32 32_to_64 sngl_to_dbl dbl_to_sngl dependencies include_stamp:
	@true
