# $Id$
# Note: To switch on disk arrays, uncomment the DEFINE +=... after the
# includes.  Do NOT move the DEFINE ahead of the includes or it won't
# have the desired effect.


LIBRARY =	librimp2g.a

OBJ =	rimp2g.o rimp2g_e_mem.o rimp2g_parm.o rimp2g_timers.o driver_e.o \
	driver_xf.o paire.o mk_t.o denom2.o util_detri.o asym_i.o e_ij.o \
	g_ij.o p2ab.o p2ij.o t_stats.o gamma3.o util_sqz_spc.o \
	ij_int.o e_final.o prep_dra.o driver_g.o mk_l.o l1.o l2.o \
	xform_3cintsg.o p2am.o mk_w2.o driver_p.o prep_sep.o clean_sep.o \
	nonsep2.o ns2_i3g3.o ns2_2e2c.o create_dra.o delete_dra.o

OBJ_OPTIMIZE = 

USES_BLAS = rimp2g_v_e2.F

include ../config/makefile.h
include ../config/makelib.h


