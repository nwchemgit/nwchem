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
		mcscf_debug.F \
		mcscf_etr.F \
		mcscf_fock.F \
		mcscf_ohvmo.F \
		mcscf_orb.F \
		mcscf_prcnd.F \
		linesrch.F \
		hmat.F


  LIB_INCLUDES = -I../ddscf -I./detci

  LIB_DEFINES = -DOLD_DAVIDSON



include ../config/makefile.h
include ../config/makelib.h









# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(linesrch.o):	 $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mcscf.o):	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh ./mcscfprof.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mcscf_aux.o):	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mcscf_ci.o):	 $(INCDIR)/rtdb.fh ./mcscfprof.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./detci/detciP.fh ./detci/detci.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mcscf_debug.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mcscf_etr.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mcscf_fcore.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mcscf_fdiff.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mcscf_fock.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh ./mcscfprof.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mcscf_hdiag.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mcscf_ohv.o):	 ./mcscfprof.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mcscf_ohv1e.o):	 ./mcscfprof.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mcscf_ohv2e.o):	 ./mcscfprof.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mcscf_ohvao.o):	 ./mcscfprof.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mcscf_ohvmo.o):	 ./mcscfprof.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mcscf_orb.o):	 $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh ./mcscfprof.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mcscf_prcnd.o):	 $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh ./mcscfprof.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
