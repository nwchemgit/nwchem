OBJ = task_gw.o dft_energy_gw.o ga_chol_inv_seq.o gw_input.o \
      gw_defaults.o gw_init.o gw_finalize.o gw_memory.o \
      gw_vxc.o mycpu.o gw_donewton.o \
      gw_cdgw.o gw_cdgw_buildin.o gw_cdgw_newton.o gw_cdgw_graph.o \
      gw_analytic.o gw_analytic_newton.o gw_analytic_graph.o

USES_BLAS = gw_init.F gw_vxc.F gw_ri_init.F gw_invdiag.F \
            ga_chol_inv_seq.F gw_tools.F \
            gw_minres.F cs_minres_qlp.F \
	    gw_cdgw.F gw_cdgw_buildiw.F gw_cdgw_buildrn_minres.F \
	    gw_cdgw_buildpi.F gw_cdgw_buildrn_graph.F \
	    gw_analytic.F gw_analytic_rpa.F

OBJ_OPTIMIZE = gw_ri_init.o gw_invdiag.o gw_tools.o \
	       gw_analytic_wmn.o gw_analytic_sigmac.o \
	       gw_analytic_addrpa.o gw_analytic_ovscale.o \
               gw_analytic_rpa.o gw_analytic_scalxpy.o \
               gw_analytic_sigmac_graph.o \
	       gw_minres.o cs_minres_qlp.o \
	       gw_cdgw_buildrn_graph.o gw_cdgw_buildiw.o \
	       gw_cdgw_buildin_graph.o gw_cdgw_buildpi.o \
	       gw_cdgw_buildrn_minres.o


LIBRARY = libnwdft.a

HEADERS = 

LIB_DEFINES =
LIB_INCLUDES = -I../include -I../grid -I./

include ../../config/makefile.h
include ../../config/makelib.h
