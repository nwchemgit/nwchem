#
# $Id$
#

 OBJ_OPTIMIZE =  esp.o\
                 esp_dalloc.o\
                 esp_denmat.o\
                 esp_elemom.o\
                 esp_esp.o\
                 esp_fastesp.o\
                 esp_getden.o\
                 esp_fit.o\
		 esp_getg.o\
                 esp_getx.o\
                 esp_grdcnt.o\
                 esp_grid.o\
                 esp_input.o\
                 esp_malloc.o\
                 esp_matinv.o\
                 esp_out.o\
                 esp_plt.o\
                 esp_rms.o\
                 esp_rtdbin.o\
                 esp_wrtinp.o

ifeq ($(FC),frt)
  OBJ = esp_jacobi.o
else
  OBJ_OPTIMIZE += esp_jacobi.o
endif


LIBRARY=libesp.a

HEADERS=

USES_BLAS=esp_esp.F esp_denmat.F

include ../config/makefile.h

include ../config/makelib.h
