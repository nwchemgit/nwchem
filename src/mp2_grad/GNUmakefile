#$Id: GNUmakefile,v 1.7 1996-07-26 16:32:51 d3g681 Exp $

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
	      mp2_vvvo_uhf.o \
	      mp2_vvoo_uhf.o \
	      mp2_vooo_uhf.o \
	      mp2_wai_uhf.o \
	      mp2_lai_fock_uhf.o \

#ifneq ($(TARGET),CRAY-T3D)
#        OBJ += crayio.o
#endif
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


