# $Id: GNUmakefile,v 1.22 1996-08-19 15:22:25 d3g681 Exp $

 OBJ = gradients.o grad_force.o grad1.o scf_gradient.o grad_getdens.o \
       grad_dens.o grad_store.o grad_inp.o
 OBJ_OPTIMIZE = 

ifeq ($(TARGET),SGITFP)
 OBJ += grad2.o
else
 OBJ_OPTIMIZE += grad2.o
endif

   USES_BLAS = grad_dens.F

   LIBRARY = libgradients.a

include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(grad_inp.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(grad_dens.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(grad_store.o):	$(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(gradients.o):	$(INCDIR)/util.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(grad_force.o):	$(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(grad2.o):	$(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(grad1.o):	$(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
