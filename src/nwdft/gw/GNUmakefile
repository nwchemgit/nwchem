OBJ = gw_input.o gw_defaults.o task_gw.o dft_energy_gw.o\
	gw_init.o gw_finalize.o \
	gw_memory.o ga_chol_inv_seq.o gw_vxc.o \
	gw_cdgw_newton.o gw_cdgw_graph.o gw_cdgw_linear.o \
	gw_buildpi.o gw_donewton.o gw_tools.o

USES_BLAS = gw_cdgw.F gw_ri_init.F gw_run.F ga_chol_inv_seq.F\
            gw_lu_solve.F gw_invdiag.F \
	    gw_buildrn.F gw_buildiw.F \
	    gw_minres.F gw_buildrn_minres.F\
	    cs_minres_qlp.F \
	    gw_analytic_rpa.F gw_cdgw_newton.F gw_en.F gw_init.F gw_vxc.F \
	    gw_buildpi.F gw_buildrn_graph.F

OBJ_OPTIMIZE = gw_ri_init.o gw_cdgw.o gw_run.o gw_lu_solve.o\
               gw_invdiag.o gw_buildrn.o gw_buildin.o\
	       gw_buildiw.o gw_analytic_wmn.o gw_analytic_sigmac.o\
	       gw_analytic_addrpa.o gw_analytic_graph.o\
	       gw_analytic_linear.o gw_analytic_newton.o\
	       gw_analytic_ovscale.o gw_analytic_rpa.o\
	       gw_analytic_scalxpy.o gw_minres.o\
	       gw_buildrn_minres.o gw_buildiw_minres.o\
	       cs_minres_qlp.o gw_en.o gw_buildrn_graph.o\
	       gw_buildin_graph.o gw_analytic_sigmac_graph.o


LIBRARY = libnwdft.a

HEADERS = 

LIB_DEFINES =
LIB_INCLUDES = -I../include -I../grid -I./

include ../../config/makefile.h
include ../../config/makelib.h
