OBJ = gwmol.o ga_chol_inv_seq.o gw_input.o \
      gw_defaults.o gw_init.o gw_finalize.o gw_memory.o \
      gw_vxc.o gw_donewton.o \
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


LIBRARY = libgwmol.a

HEADERS = ../nwdft/include/cdft.fh ../nwdft/grid/cgridfile.fh 

LIB_DEFINES =
LIB_INCLUDES = -I./

include ../config/makefile.h

ifeq ($(_USE_SCALAPACK),Y)
_USE_ELPA = $(shell ${GA_PATH}/bin/ga-config  --use_elpa| awk ' /1/ {print "Y"}')

ifeq ($(_USE_ELPA),Y)
 LIB_DEFINES += -DUSE_ELPA
endif
#assuming we have GAPZHEEVD
GOTGAPZHEEVD= Y
ifeq ($(GOTGAPZHEEVD),Y)
 LIB_DEFINES += -DUSE_GAPZHEEVD
endif
_PDSYEVR = $(shell [ -e ${NWCHEM_TOP}/src/tools/install/lib/libga.a ] &&  nm ${NWCHEM_TOP}/src/tools/install/lib/libga.a|egrep T |egrep  ga_pdsyevr |awk '/ga_pdsyevr/ {print "Y"}')
ifeq ($(_PDSYEVR),Y) 
 LIB_DEFINES += -DUSE_PDSYEVR
endif
endif
include ../config/makelib.h
