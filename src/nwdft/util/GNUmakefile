#$Id: GNUmakefile,v 1.43 1996-07-22 06:04:30 d3e129 Exp $
# makefile,v 1.15 1994/12/14 02:13:59 og845 Exp
ifneq ($(TARGET),SGITFP)
ifeq ($(TARGET),CRAY-T3D)
	OBJ = ga_chol_seq.o ga_inv_seq.o dft_invio.o
	OBJ_OPTIMIZE = \
	new_eval_gb.o \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	gacrea.o		\
	dft_gather.o	dft_scatter.o \
	get_col.o	\
	dft_header.o	\
	pertbl.o	\
	sca_lab_mat.o	\
	dft_second.o	\
	dft_tstamp.o	\
	bas_ncnmax.o       nbf_to_mbf.o \
	dft_mpole.o  movecs_log.o
else
	OBJ_OPTIMIZE = \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	gacrea.o		\
	get_col.o	\
	dft_header.o	\
	dft_invio.o		\
	pertbl.o	\
	dft_gather.o	sca_lab_mat.o	dft_scatter.o	\
	dft_second.o	\
	dft_tstamp.o	\
	new_eval_gb.o \
        dft_eval_bas.o        dft_eval_fns.o   \
        dft_eval_grd.o        dft_eval_hss.o    \
	bas_ncnmax.o		ga_chol_seq.o	\
	ga_inv_seq.o 	\
	nbf_to_mbf.o dft_mpole.o movecs_log.o
endif
else
	OBJ = bas_ncnmax.o 	nbf_to_mbf.o new_eval_gb.o         dft_eval_bas.o        dft_eval_fns.o   \
        dft_eval_grd.o        dft_eval_hss.o dft_mpole.o


	OBJ_OPTIMIZE=  \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	gacrea.o	\
	get_col.o	\
	dft_header.o	\
	dft_invio.o		\
	pertbl.o	\
	dft_gather.o	sca_lab_mat.o	dft_scatter.o	\
	dft_second.o	\
	dft_tstamp.o	\
	ga_chol_seq.o	\
	ga_inv_seq.o 	movecs_log.o
endif

     LIBRARY = libnwdft.a
LIB_INCLUDES = -I../include


USES_BLAS = dft_eval_fns.F dft_eval_grd.F dft_eval_hss.F dft_mpole.F dft_second.F ga_chol_seq.F ga_inv_seq.F new_eval_gb.F

include ../../config/makefile.h
include ../../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(movecs_log.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/cscfps.fh 
$(LIBRARY_PATH)(dft_eval_hss.o):	$(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_eval_grd.o):	$(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_eval_fns.o):	$(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_eval_bas.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(bas_ncnmax.o):	$(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_invio.o):	$(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(new_eval_gb.o):	$(INCDIR)/bas.fh 
$(LIBRARY_PATH)(ga_inverse.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_inv_seq.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_chol_seq.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(ga_chol.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_second.o):	../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_mpole.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_header.o):	$(INCDIR)/stdio.fh 
$(LIBRARY_PATH)(gacrea.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_tstamp.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/global.fh 
