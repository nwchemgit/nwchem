# $Id: GNUmakefile,v 1.14 1995-11-19 19:41:22 d3g681 Exp $

ifeq ($(NWCHEM_TARGET),SGITFP)
     OBJ = gradients.o grad_force.o grad1.o scf_gradient.o grad_get_dens.o \
           grad_dens.o grad2.o
else
     OBJ = gradients.o grad_force.o grad1.o scf_gradient.o grad_get_dens.o \
           grad_dens.o 
   OBJ_OPTIMIZE = grad2.o
endif

   USES_BLAS = grad_force.F

   LIBRARY = libgradients.a

include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(grad1.o):	 $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(grad2.o):	 $(INCDIR)/sym.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(grad_dens.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(grad_force.o):	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(grad_get_dens.o):	 $(INCDIR)/global.fh 
$(LIBRARY_PATH)(gradients.o):	 $(INCDIR)/util.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
