#$Id: GNUmakefile,v 1.58 1999-03-13 04:14:36 edo Exp $
# makefile,v 1.15 1994/12/14 02:13:59 og845 Exp
ifneq ($(TARGET),SGITFP)
ifeq ($(TARGET),CRAY-T3D)
	OBJ = ga_chol_seq.o ga_inv_seq.o dft_invio.o
	OBJ_OPTIMIZE = \
	new_eval_gb.o \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	dft_gather.o	dft_scatter.o \
	get_col.o	\
	dft_header.o	\
	pertbl.o	\
	sca_lab_mat.o	\
	dft_second.o	\
	dft_tstamp.o	\
	nbf_to_mbf.o \
	dft_mpole.o  movecs_log.o dft_utilmap.o
else
	OBJ_OPTIMIZE = \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	get_col.o	\
	dft_header.o	\
	dft_invio.o		\
	pertbl.o	\
	dft_gather.o	sca_lab_mat.o	dft_scatter.o	\
	dft_second.o	\
	dft_tstamp.o	\
	new_eval_gb.o \
	ga_chol_seq.o	\
	ga_inv_seq.o 	\
	nbf_to_mbf.o dft_mpole.o movecs_log.o dft_utilmap.o
endif
else
#	OBJ = nbf_to_mbf.o new_eval_gb.o         dft_eval_bas.o        dft_eval_fns.o   \
#        dft_eval_grd.o        dft_eval_hss.o dft_mpole.o
#
#
#	OBJ_OPTIMIZE=  \
#	dzero.o		\
#	ga_chol.o	\
#	ga_inverse.o		\
#	gacrea.o	\
#	get_col.o	\
#	dft_header.o	\
#	dft_invio.o		\
#	pertbl.o	\
#	dft_gather.o	sca_lab_mat.o	dft_scatter.o	\
#	dft_second.o	\
#	dft_tstamp.o	\
#	ga_chol_seq.o	\
#	ga_inv_seq.o 	movecs_log.o
#

	OBJ_OPTIMIZE=  \
	dft_gather.o \
	dft_header.o \
	dft_invio.o \
	dft_mpole.o \
	dft_scatter.o \
	dft_second.o \
	dft_tstamp.o \
	dzero.o \
	ga_chol.o \
	ga_chol_seq.o \
	ga_inv_seq.o \
	ga_inverse.o \
	get_col.o \
	movecs_log.o \
	nbf_to_mbf.o \
	new_eval_gb.o \
	pertbl.o \
	sca_lab_mat.o dft_utilmap.o

endif

     LIBRARY = libnwdft.a

# LIB_DEFINES = -DDFT_TIMER
 LIB_DEFINES =

LIB_INCLUDES = -I../include


USES_BLAS = dft_second.F ga_chol_seq.F ga_inv_seq.F new_eval_gb.F

include ../../config/makefile.h
include ../../config/makelib.h

