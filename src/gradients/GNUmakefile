# $Id: GNUmakefile,v 1.8 1995-10-17 08:53:53 d3g681 Exp $

            OBJ = gradients.o force.o grad1.o scf_gradient.o get_dens.o
   OBJ_OPTIMIZE = grad2.o

   USES_BLAS = force.F

   LIBRARY = libgradients.a

include ../config/makefile.h
include ../config/makelib.h
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
force.o :	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
grad1.o :	 $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
grad2.o :	 $(INCDIR)/sym.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
gradients.o :	 $(INCDIR)/util.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
