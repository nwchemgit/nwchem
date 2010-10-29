#
# $Id$
#
# intf and aux break on IBM SP xlhpf with optimization ... why? 
 OBJ_OPTIMIZE = \
	      	detci_dens.o \
	      	detci_diag.o \
	      	detci_sigma.o \
		detci_spin.o


	  OBJ = detci_aux.o \
	      	detci_intf.o 

      LIBRARY = libmcscf.a

    USES_BLAS = detci_aux.F \
	 	detci_dens.F \
	        detci_diag.F \
		detci_sigma.F \
		detci_spin.F

  LIB_DEFINES = 


include ../../config/makefile.h
include ../../config/makelib.h

