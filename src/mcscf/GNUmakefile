      SUBDIRS =	detci 


          OBJ = mcscf.o \
		mcscf_aux.o \
		mcscf_cisolve.o \
		mcscf_etrace.o \
		mcscf_fcore.o \
	        mcscf_fock.o \
		mcscf_fdiff.o \
		mcscf_getmovecs.o \
		mcscf_gfock2grad.o \
		mcscf_hdiag.o \
		mcscf_hdiag_scale.o \
		mcscf_hessv.o \
		mcscf_hessv_1e.o \
		mcscf_hessv_2e.o \
		mcscf_hessv_2e_ao.o \
		mcscf_hessv_2e_mo.o \
		mcscf_orb.o \
		mcscf_precond.o \
		mcscf_rohf_modens.o \
		linesrch.o \
		debug_crap.o \
		rrrrr.o \
		check_moblocked.o \
	        diagfill_patch.o

      LIBRARY =	libmcscf.a

   USES_BLAS  =	mcscf.F \
		mcscf_aux.F \
		mcscf_dbg.F \
		mcscf_dens2mo.F \
		mcscf_etrace.F \
		mcscf_fock.F \
		mcscf_gfock2grad.F \
		mcscf_hessv_2e_mo.F \
		mcscf_orb.F \
		mcscf_precond.F


  LIB_INCLUDES = -I../ddscf -I./detci

  LIB_DEFINES =



include ../config/makefile.h
include ../config/makelib.h














