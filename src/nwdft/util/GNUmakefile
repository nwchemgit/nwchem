#$Id$

	OBJ=	ga_chol.o ga_inverse.o ga_chol_seq.o ga_inv_seq.o dft_invio.o dft_print_mats.o dft_save_mats.o


	OBJ_OPTIMIZE=  \
	dft_gather.o \
	dft_header.o \
	dft_mpole.o \
	dft_scatter.o \
	dft_tstamp.o \
        dft_pstrace.o \
	get_col.o \
	movecs_log.o \
	nbf_to_mbf.o \
	dft_frozemb.o \
	dft_frozemb_fde.o \
	int_1e_fde.o \
	sca_lab_mat.o \
	dft_utilmap.o \
	dft_genutils.o


     LIBRARY = libnwdft.a

 LIB_DEFINES =

LIB_INCLUDES = -I../include


USES_BLAS = ga_chol_seq.F ga_inv_seq.F int_1e_fde.F

include ../../config/makefile.h
include ../../config/makelib.h

