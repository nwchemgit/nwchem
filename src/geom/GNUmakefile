#$Id: GNUmakefile,v 1.11 1995-11-16 20:38:20 d3g681 Exp $

 LIBRARY = libutil.a
     OBJ = geom.o geom_input.o geom_3d.o geom_num_core.o
 HEADERS = geom.fh geomP.fh

 USES_BLAS = geom.F geom_input.F

include ../config/makefile.h
include ../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
geom.o :	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh geomP.fh $(INCDIR)/context.fh 
geom_3d.o :	 $(INCDIR)/inp.fh geomP.fh geom.fh 
geom_input.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/inp.fh geomP.fh geom.fh 
geom_num_core.o :	 geomP.fh 
