# $Id: GNUmakefile,v 1.24 1995-10-30 22:01:45 gg502 Exp $
# Note: To switch on disk arrays, uncomment the DEFINE +=... after the
# includes.  Do NOT move the DEFINE ahead of the includes or it won't
# have the desired effect.


LIBRARY =	librimp2.a

OBJ =		rimp2_v_e2.o \
		mk_1d_shell_dist.o	ga_create_one_by_shells.o \
		xf3ci_step12.o	xf3ci_step3.o xform_3cints.o \
		occ_ranges.o	rimp2_xf_driver.o \
		close_3ci_file.o	open_3ci_file.o	put_3cints.o \
		get_3cints.o	rimp2.o	\
		dgecop.o	lnblnk.o \
		dgewr.o		atoi.o		indint.o \
		max_v_diagonal.o	rimp2_get_param.o \
		rimp2_timers.o  \
                mk_correl_map.o print_correl_map.o      progress_report.o \
		xf3ci_mem12.o	iamax.o		rimp2_e_mem.o

OBJ_OPTIMIZE =      rimp2_denom.o   \
                sum_pairs.o     xf3ci_step12a.o \
		xf3ci_step12b.o

ifeq ($(NWCHEM_TARGET),SGITFP)
   OBJ += eval_3ci_block.o
else
   OBJ_OPTIMIZE += eval_3ci_block.o
endif

USES_BLAS =	rimp2_v_e2.F    xf3ci_step12.F  xf3ci_step12a.F \
		xf3ci_step12b.F xform_3cints.F open_3ci_file.F \
		rimp2.F eval_3ci_block.F

include ../config/makefile.h
include ../config/makelib.h

# The default is to use standard Fortran direct access files.
# Uncomment this definition to switch to use disk arrays instead.
# You must also build and link against the DA library (see ../global/da).
#DEFINES += -DRIMP2_USE_DA

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
close_3ci_file.o :	 $(INCDIR)/stdio.fh 
dgewr.o :	 output.h 
eval_3ci_block.o :	 $(INCDIR)/schwarz.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/bas.fh 
ga_create_one_by_shells.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
get_3cints.o :	 $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
max_v_diagonal.o :	 $(INCDIR)/stdio.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
mk_1d_shell_dist.o :	 $(INCDIR)/stdio.fh $(INCDIR)/bas.fh 
mk_correl_map.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh 
mk_v.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
occ_ranges.o :	 $(INCDIR)/numerical_constants.fh 
open_3ci_file.o :	 $(INCDIR)/util.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ftnio.fh 
progress_report.o :	 $(INCDIR)/global.fh 
put_3cints.o :	 $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
rimp2.o :	 $(INCDIR)/util.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh rimp2_timers.fh $(INCDIR)/pstat.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
rimp2_denom.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
rimp2_e_mem.o :	 $(INCDIR)/stdio.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
rimp2_get_param.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh 
rimp2_timers.o :	 rimp2_timers.fh $(INCDIR)/pstat.fh 
rimp2_v_e2.o :	 $(INCDIR)/util.fh $(INCDIR)/stdio.fh rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
rimp2_xf_driver.o :	 $(INCDIR)/stdio.fh rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
sum_pairs.o :	 $(INCDIR)/numerical_constants.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
xf3ci_step12.o :	 $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
xf3ci_step12a.o :	 rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
xf3ci_step12b.o :	 $(INCDIR)/schwarz.fh rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
xf3ci_step3.o :	 $(INCDIR)/stdio.fh rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
xform_3cints.o :	 $(INCDIR)/util.fh $(INCDIR)/stdio.fh rimp2_timers.fh $(INCDIR)/numerical_constants.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
