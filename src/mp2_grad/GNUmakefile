#$Id: GNUmakefile,v 1.11 1996-10-02 15:22:13 d3e129 Exp $

        OBJ = mp2_make_trhf.o \
	      mp2_make_tuhf.o \
	      mp2_grad.o \
              mp2_energy.o \
	      mp2_grad_doit.o \
	      mp2_read_tuhf.o \
	      mp2_read_trhf.o \
	      mp2_pijab_uhf.o \
	      mp2_wijab_uhf.o \
	      mp2_lai_uhf.o \
	      mp2_vvoo_uhf.o \
	      mp2_wai_uhf.o \
	      mp2_lai_fock_uhf.o \
	      mp2_make_cphf_rhs.o \
	      mp2_make_tot_dens.o \
	      mp2_make_tot_wdens.o \
	      mp2_make_dip_mom.o \
              mp2_back_transform.o

USES_BLAS = mp2_energy.F mp2_lai_uhf.F mp2_make_t.F mp2_make_trhf.F \
            mp2_make_tuhf.F mp2_pijab_uhf.F mp2_vooo_uhf.F \
            mp2_vvoo_uhf.F mp2_vvvo_uhf.F mp2_wai_uhf.F mp2_wijab_uhf.F

    LIBRARY = libmp2.a


LIB_DEFINES =  -DDEBUG_PRINT


# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */


include ../config/makefile.h
include ../config/makelib.h


