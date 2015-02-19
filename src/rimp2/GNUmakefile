# $Id$
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
		max_v_diag.o	rimp2_parm.o \
		rimp2_timers.o  \
                pr_corr_map.o      prog_report.o \
		xf3ci_mem12.o	iamax.o		rimp2_e_mem.o

OBJ_OPTIMIZE =      rimp2_denom.o   \
                sum_pairs.o     xf3ci_s12a.o \
		xf3ci_s12b.o

   OBJ_OPTIMIZE += eval_3ci.o

USES_BLAS = dgecop.F eval_3ci.F rimp2_v_e2.F xf3ci_s12a.F xf3ci_s12b.F

ifdef ESSL_TRANSPOSE
   DEFINES += -DESSL_TRANSPOSE -DBLAS_COPY
else
   ifdef OSU_TRANSPOSE
      DEFINES += -DOSU_TRANSPOSE -DBLAS_COPY
      OBJ_OPTIMIZE += osu_transpose.o
   endif
endif

include ../config/makefile.h
include ../config/makelib.h

# The default is to use standard Fortran direct access files.
# Uncomment this definition to switch to use disk arrays instead.
# You must also build and link against the DA library (see ../global/da).
#DEFINES += -DRIMP2_USE_DA
#DEFINES += -DDEBUG
