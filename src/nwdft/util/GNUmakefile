#$Id: GNUmakefile,v 1.61 2001-10-18 20:10:20 edo Exp $

	OBJ=	ga_chol.o ga_inverse.o ga_chol_seq.o ga_inv_seq.o dft_invio.o


	OBJ_OPTIMIZE=  \
	dft_gather.o \
	dft_header.o \
	dft_mpole.o \
	dft_scatter.o \
	dft_second.o \
	dft_tstamp.o \
	dzero.o \
	get_col.o \
	movecs_log.o \
	nbf_to_mbf.o \
	new_eval_gb.o \
	sca_lab_mat.o dft_utilmap.o


     LIBRARY = libnwdft.a

# LIB_DEFINES = -DDFT_TIMER
 LIB_DEFINES =

LIB_INCLUDES = -I../include


USES_BLAS = dft_second.F ga_chol_seq.F ga_inv_seq.F new_eval_gb.F

include ../../config/makefile.h
include ../../config/makelib.h

