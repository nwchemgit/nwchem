# $Id: GNUmakefile,v 1.33 1996-08-27 21:56:18 d3e129 Exp $
# Note: To switch on disk arrays, uncomment the DEFINE +=... after the
# includes.  Do NOT move the DEFINE ahead of the includes or it won't
# have the desired effect.


LIBRARY =	librimp2.a

OBJ =		rimp2_v_e2.o \
		mk_1d_shdist.o	ga_by_shells.o \
		xf3ci_step12.o	xf3ci_step3.o xform_3cints.o \
		occ_ranges.o	rimp2_xf.o \
		close_3ci.o	open_3ci.o	put_3cints.o \
		get_3cints.o	rimp2.o	\
		dgecop.o	lnblnk.o \
		dgewr.o		atoi.o		indint.o \
		max_v_diag.o	rimp2_parm.o \
		rimp2_timers.o  \
                mk_corr_map.o pr_corr_map.o      prog_report.o \
		xf3ci_mem12.o	iamax.o		rimp2_e_mem.o

OBJ_OPTIMIZE =      rimp2_denom.o   \
                sum_pairs.o     xf3ci_s12a.o \
		xf3ci_s12b.o

ifeq ($(TARGET),SGITFP)
   OBJ += eval_3ci.o
else
   OBJ_OPTIMIZE += eval_3ci.o
endif

USES_BLAS = eval_3ci.F rimp2_v_e2.F xf3ci_s12a.F xf3ci_s12b.F xf3ci_step3.F

include ../config/makefile.h
include ../config/makelib.h

# The default is to use standard Fortran direct access files.
# Uncomment this definition to switch to use disk arrays instead.
# You must also build and link against the DA library (see ../global/da).
#DEFINES += -DRIMP2_USE_DA

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(xform_3cints.o):	$(INCDIR)/util.fh $(INCDIR)/stdio.fh ./rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(xf3ci_s12b.o):	$(INCDIR)/schwarz.fh ./rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(xf3ci_s12a.o):	./rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(rimp2_parm.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(rimp2_e_mem.o):	$(INCDIR)/stdio.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rimp2.o):	$(INCDIR)/util.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh ./rimp2_timers.fh $(INCDIR)/pstat.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(rimp2_denom.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(open_3ci.o):	$(INCDIR)/util.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./ftnio.fh 
$(LIBRARY_PATH)(eval_3ci.o):	$(INCDIR)/schwarz.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(max_v_diag.o):	$(INCDIR)/stdio.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mk_v.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(get_3cints.o):	$(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/da.fh 
$(LIBRARY_PATH)(put_3cints.o):	$(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/da.fh 
$(LIBRARY_PATH)(xf3ci_step3.o):	$(INCDIR)/stdio.fh ./rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mk_1d_shdist.o):	$(INCDIR)/stdio.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(occ_ranges.o):	$(INCDIR)/numerical_constants.fh 
$(LIBRARY_PATH)(rimp2_timers.o):	./rimp2_timers.fh $(INCDIR)/pstat.fh 
$(LIBRARY_PATH)(rimp2_xf.o):	$(INCDIR)/util.fh $(INCDIR)/stdio.fh ./rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(sum_pairs.o):	$(INCDIR)/numerical_constants.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dgewr.o):	./output.h 
$(LIBRARY_PATH)(close_3ci.o):	$(INCDIR)/stdio.fh 
$(LIBRARY_PATH)(mk_corr_map.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh 
$(LIBRARY_PATH)(ga_by_shells.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(prog_report.o):	$(INCDIR)/global.fh 
$(LIBRARY_PATH)(xf3ci_step12.o):	$(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh ./rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(rimp2_v_e2.o):	$(INCDIR)/util.fh $(INCDIR)/stdio.fh ./rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
