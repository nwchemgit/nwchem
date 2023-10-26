## $Id$

include ../../config/makefile.h
include ../DEFS

     CNFDIR := $(NWCHEM_TOP)/src/config
.PHONY:	 64_to_32 32_to_64
64_to_32:
	$(CNFDIR)/64_to_32 blas_lapack.h peigs_types.h

32_to_64:
	$(CNFDIR)/32_to_64 blas_lapack.h peigs_types.h
