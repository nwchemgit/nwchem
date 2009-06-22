#
# $Id: GNUmakefile,v 1.32 2008-10-01 23:50:12 bylaska Exp $
#


   OBJ_OPTIMIZE = hnd_efgmap.o \
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
         hnd_psoten.o \
         hnd_dsoten.o \
         prop_atom.o  \
         aoresponse_driver.o \
         int_mpole_1ega.o \
         aoresponse.o \
         CalcPerturbedTDPmat1.o \
         aoresponse_giao_rhs.o \
         ao_1prdm_read.o

   LIBRARY = libproperty.a
   LIB_DEFINES = -DJEFF_OUTPUT

#
# If you add -DJEFF_OUTPUT to LIB_DEFINES, the dipole polarizability will be printed in the same format as the TCE.
# This is very useful (1) if you wrote scripts to parse that particular formatting and (2) if you want the anisotrophy.
#

   LIB_INCLUDES = -I../nwdft/include -I../ddscf

USES_BLAS = hnd_elfcon.F wnbofile.F hnd_spinspin.F hnd_mtpole.F hnd_mulken.F \
            hnd_elpmap.F hnd_elfmap.F hnd_mtpcon.F giaofock.F hnd_spinspin.F \
            aoresponse_driver.F

HEADERS = prop.fh atomdata.fh

include ../config/makefile.h
include ../config/makelib.h

#FDEBUG += -O0

