#
# $Id: GNUmakefile,v 1.23 2004-12-16 18:47:22 bert Exp $
#


   OBJ = hnd_efgmap.o \
         hnd_elfmap.o \
         hnd_elpmap.o \
         hnd_eldmap.o \
         hnd_spind.o \
         hnd_elfcon.o \
         hnd_mtpcon.o \
         hnd_mtpole.o \
         prop_input.o \
         prop.o \
         prp.o \
         hnd_mulken.o \
         wnbofile.o \
         hnd_prp_get_dens.o \
         hnd_prp_vec_read.o \
         hnd_prop_dens_make.o \
         hnd_property.o \
         cosmo.o \
         giaofock.o \
         hnd_giaox.o \
         int_giao_1ega.o \
         hnd_diag.o \
         hnd_spinspin.o \
         hnd_sdfc.o \
         hnd_pso.o \
         hnd_dso.o \
         prop_atom.o \
         xc_tabcd_giao_nmr_rhs.o

   LIBRARY = libproperty.a
   LIB_DEFINES = 

    LIB_INCLUDES = -I../nwdft/include -I../ddscf

USES_BLAS = 

HEADERS = prop.fh atomdata.fh

include ../config/makefile.h
include ../config/makelib.h

#FDEBUG += -O0

