#$Id: GNUmakefile,v 1.66 2004-08-26 23:31:22 edo Exp $

	OBJ=	ga_chol.o ga_inverse.o ga_chol_seq.o ga_inv_seq.o dft_invio.o


	OBJ_OPTIMIZE=  \
	dft_gather.o \
	dft_header.o \
	dft_mpole.o \
	dft_scatter.o \
	dft_tstamp.o \
	get_col.o \
	movecs_log.o \
	nbf_to_mbf.o \
	sca_lab_mat.o dft_utilmap.o


     LIBRARY = libnwdft.a

 LIB_DEFINES =

LIB_INCLUDES = -I../include


USES_BLAS = ga_chol_seq.F ga_inv_seq.F 

include ../../config/makefile.h
include ../../config/makelib.h

