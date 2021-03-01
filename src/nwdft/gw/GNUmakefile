OBJ = gw_input.o gw_defaults.o task_gw.o dft_energy_gw.o\
	gw_init.o gw_finalize.o \
	gw_memory.o ga_chol_inv_seq.o gw_vxc.o \
	gw_cdgw_newton.o gw_cdgw_graph.o gw_cdgw_linear.o

USES_BLAS = gw_cdgw.F gw_ri_init.F gw_run.F ga_chol_inv_seq.F\
            gw_lu_solve.F gw_invdiag ga_chol_inv_seq.F\
	    gw_buildrn.o gw_buildiw.o gw_analytic_sigmac.o\
	    gw_analytic_wmn.o gw_minres.o gw_buildrn_minres.o\
	    cs_minres_qlp.o

OBJ_OPTIMIZE = gw_ri_init.o gw_cdgw.o gw_run.o gw_lu_solve.o\
               gw_invdiag.o gw_buildrn.o gw_buildin.o\
	       gw_buildiw.o gw_analytic_wmn.o gw_analytic_sigmac.o\
	       gw_analytic_addrpa.o gw_analytic_graph.o\
	       gw_analytic_linear.o gw_analytic_newton.o\
	       gw_analytic_ovscale.o gw_analytic_rpa.o\
	       gw_analytic_scalxpy.o gw_minres.o\
	       gw_buildrn_minres.o gw_buildiw_minres.o\
	       cs_minres_qlp.o gw_en.o


LIBRARY = libnwdft.a

HEADERS = 

LIB_DEFINES =
LIB_INCLUDES = -I../include -I../grid -I./

include ../../config/makefile.h
include ../../config/makelib.h
