
          OBJ = detci_aux.o \
	      	detci_dens.o \
	      	detci_dvddiag.o \
	      	detci_intface.o \
	      	detci_sigma.o \
		detci_spin.o

      LIBRARY = libmcscf.a

    USES_BLAS = detci_aux.F \
	        detci_dvddiag.F \
		detci_intface.F \
		detci_sigma.F \
		detci_spin.F


include ../../config/makefile.h
include ../../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
detci_aux.o :	 $(INCDIR)/global.fh detciP.fh detci.fh 
detci_dens.o :	 $(INCDIR)/global.fh detciP.fh detci.fh 
detci_driver.o :	 detciP.fh detci.fh 
detci_dvddiag.o :	 $(INCDIR)/global.fh 
detci_intface.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh detciP.fh detci.fh cdetcistats.fh cdetci.fh 
detci_sigma.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh detciP.fh detci.fh 
detci_spin.o :	 detciP.fh detci.fh 
