#
# $Id: GNUmakefile,v 1.5 1999-07-28 00:39:07 d3e129 Exp $
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

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:		
		-rm -f *.o

.c.o:	 ;	$(peigs_CC) -c $<
 
