      SUBDIRS =	detci 


          OBJ = mcscf.o \
		mcscf_aux.o \
		mcscf_ci.o \
		mcscf_etr.o \
		mcscf_fcore.o \
	        mcscf_fock.o \
		mcscf_fdiff.o \
		mcscf_hdiag.o \
		mcscf_ohv.o \
		mcscf_ohv1e.o \
		mcscf_ohv2e.o \
		mcscf_ohvao.o \
		mcscf_ohvmo.o \
		mcscf_orb.o \
		mcscf_prcnd.o \
		mcscf_debug.o \
		linesrch.o \
		hmat.o

      LIBRARY =	libmcscf.a

   USES_BLAS  =	mcscf.F \
		mcscf_aux.F \
		mcscf_dbg.F \
		mcscf_etr.F \
		mcscf_fock.F \
		mcscf_ohvmo.F \
		mcscf_orb.F \
		mcscf_prcnd.F \
		hmat.F


  LIB_INCLUDES = -I../ddscf -I./detci

  LIB_DEFINES =



include ../config/makefile.h
include ../config/makelib.h














