# $Id: GNUmakefile,v 1.11 1995-11-13 16:56:40 d3g681 Exp $

            OBJ = gradients.o grad_force.o grad1.o scf_gradient.o grad_get_dens.o \
			grad_dens.o grad2.o
   OBJ_OPTIMIZE = grad2.o

   USES_BLAS = grad_force.F

   LIBRARY = libgradients.a

include ../config/makefile.h
include ../config/makelib.h
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
grad1.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
grad2.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
grad_dens.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
grad_force.o :	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
grad_get_dens.o :	 $(INCDIR)/global.fh 
gradients.o :	 $(INCDIR)/util.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
