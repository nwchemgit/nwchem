#
# $Id$
#

include ../../config/makefile.h
include ../DEFS

COBJ  = bortho_f.o \
        mxm25_f.o \
        mxm5x_f.o \
        mxm88_f.o \
        onenorm_f.o \
        ortho_f.o \
        pdspev_f.o \
        pdspevx_f.o \
        pdspgv_f.o \
        pdspgvx_f.o \
        resid_f.o \
        residual_f.o \
        sonenorm_f.o \
        tresid_f.o\
	choleski_f.o\
	inverse_f.o

#        pdsptri_f.o \

all:		$(COBJ)

#.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:		
		-rm -f *.o

 CNFDIR := $(NWCHEM_TOP)/src/config
.PHONY:	 64_to_32 32_to_64
64_to_32:
	$(CNFDIR)/64_to_32 "*c *f *F"

32_to_64:
	$(CNFDIR)/32_to_64 "*c *f *F"

.c.o:	 ;	$(peigs_CC) -c $<
 
