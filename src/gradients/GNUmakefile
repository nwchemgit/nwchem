# $Id: GNUmakefile,v 1.9 1995-11-07 20:00:07 vg038 Exp $

            OBJ = gradients.o grad_force.o grad1.o scf_gradient.o grad_get_dens.o \
			grad_dens.o grad2.o
#  OBJ_OPTIMIZE = grad2.o

   USES_BLAS = grad_force.F

   LIBRARY = libgradients.a

include ../config/makefile.h
include ../config/makelib.h
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
grad_force.o :	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
grad1.o :	 $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
grad2.o :	 $(INCDIR)/sym.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
gradients.o :	 $(INCDIR)/util.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
