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
		mcscf_debug.o \
		linesrch.o \
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














# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
check_moblocked.o :	 $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
diagfill_patch.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
linesrch.o :	 $(INCDIR)/global.fh 
mcscf.o :	 $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
mcscf_aux.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_cisolve.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/detciP.fh $(INCDIR)/detci.fh $(INCDIR)/bas.fh 
mcscf_debug.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_etrace.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
mcscf_fcore.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_fdiff.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
mcscf_fock.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
mcscf_getmovecs.o :	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
mcscf_gfock2grad.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_hdiag.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_hdiag_scale.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_hessv.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_hessv_1e.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_hessv_2e.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_hessv_2e_ao.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_hessv_2e_mo.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
mcscf_orb.o :	 $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
mcscf_precond.o :	 $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
