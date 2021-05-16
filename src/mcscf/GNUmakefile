#
# $Id$
#
      SUBDIRS =	detci 


          OBJ = mcscf.o \
		mcscf_b.o \
		mcscf_aux.o \
		mcscf_ci.o \
		mcscf_etr.o \
		mcscf_fcore.o \
	        mcscf_fock.o \
		mcscf_hdiag.o \
		mcscf_ohv.o \
		mcscf_ohv1e.o \
		mcscf_ohv2e.o \
		mcscf_ohvao.o \
		mcscf_ohvmo.o \
		mcscf_orb.o \
		mcscf_prcnd.o \
		mcscf_debug.o \
                mcscf_input.o


      LIBRARY =	libmcscf.a

   USES_BLAS = hmat.F mcscf_aux.F mcscf_b.F mcscf_debug.F mcscf_fock.F mcscf_ohvmo.F mcscf_orb.F mcscf_prcnd.F mcscf.F

  LIB_INCLUDES = -I../ddscf -I./detci

  LIB_DEFINES = -DOLD_DAVIDSON




include ../config/makefile.h
include ../config/makelib.h









