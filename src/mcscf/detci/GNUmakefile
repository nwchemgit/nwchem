
 OBJ_OPTIMIZE = detci_aux.o \
	      	detci_dens.o \
	      	detci_diag.o \
	      	detci_intf.o \
	      	detci_sigma.o \
		detci_spin.o

      LIBRARY = libmcscf.a

    USES_BLAS = detci_aux.F \
	 	detci_dens.F \
	        detci_diag.F \
		detci_sigma.F \
		detci_spin.F

  LIB_DEFINES = 


include ../../config/makefile.h
include ../../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(detci_spin.o):	$(INCDIR)/global.fh ./detciP.fh ./detci.fh 
$(LIBRARY_PATH)(detci_sigma.o):	$(INCDIR)/util.fh ./symmmul.fh ./symmdef.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./detciP.fh ./detci.fh ./cdetcistats.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(detci_intf.o):	$(INCDIR)/util.fh ./symmmul.fh ./symmdef.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./detciP.fh ./detci.fh ./cdetcistats.fh ./cdetci.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(detci_driver.o):	$(INCDIR)/util.fh ./detciP.fh ./detci.fh 
$(LIBRARY_PATH)(detci_diag.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(detci_dens.o):	$(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./detciP.fh ./detci.fh 
$(LIBRARY_PATH)(detci_aux.o):	./symmmul.fh ./symmdef.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./detciP.fh ./detci.fh $(INCDIR)/bitops.fh 
