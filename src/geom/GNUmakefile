#$Id: GNUmakefile,v 1.16 1996-07-14 18:17:35 gutowski Exp $

 LIBRARY = libutil.a
     OBJ = geom.o geom_input.o geom_3d.o geom_2d.o geom_1d.o geom_numcore.o
 HEADERS = geom.fh geomP.fh

 USES_BLAS = geom.F geom_input.F

include ../config/makefile.h
include ../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(geom.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh ./geomP.fh $(INCDIR)/context.fh 
$(LIBRARY_PATH)(geom_3d.o):	 $(INCDIR)/inp.fh ./geomP.fh ./geom.fh 
$(LIBRARY_PATH)(geom_input.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/inp.fh ./geomP.fh ./geom.fh 
$(LIBRARY_PATH)(geom_numcore.o):	 ./geomP.fh 
