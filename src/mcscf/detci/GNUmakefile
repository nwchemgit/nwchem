
          OBJ = detci_aux.o \
	      	detci_dens.o \
	      	detci_diag.o \
	      	detci_intf.o \
	      	detci_sigma.o \
		detci_spin.o

      LIBRARY = libmcscf.a

    USES_BLAS = detci_aux.F \
	        detci_diag.F \
		detci_intf.F \
		detci_sigma.F \
		detci_spin.F


include ../../config/makefile.h
include ../../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(detci_aux.o):	 $(INCDIR)/global.fh ./detciP.fh ./detci.fh 
$(LIBRARY_PATH)(detci_dens.o):	 $(INCDIR)/global.fh ./detciP.fh ./detci.fh 
$(LIBRARY_PATH)(detci_diag.o):	 $(INCDIR)/util.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(detci_driver.o):	 $(INCDIR)/util.fh ./detciP.fh ./detci.fh 
$(LIBRARY_PATH)(detci_intf.o):	 $(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./detciP.fh ./detci.fh ./cdetcistats.fh ./cdetci.fh 
$(LIBRARY_PATH)(detci_sigma.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./detciP.fh ./detci.fh 
$(LIBRARY_PATH)(detci_spin.o):	 ./detciP.fh ./detci.fh 
