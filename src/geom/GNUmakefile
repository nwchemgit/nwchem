#$Id: GNUmakefile,v 1.25 1996-11-14 22:54:04 mdupuis Exp $

 LIBRARY = libutil.a

     OBJ = geom.o geom_input.o geom_3d.o geom_2d.o geom_1d.o geom_numcore.o \
           geom_checksum.o

 HEADERS = geom.fh geomP.fh

 USES_BLAS = geom.F geom_input.F

include ../config/makefile.h
include ../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(geom_input.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgtypesf.h $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./geomP.fh ./geom.fh $(INCDIR)/context.fh 
$(LIBRARY_PATH)(geom_checksum.o):	$(INCDIR)/mafdecls.fh ./geomP.fh 
$(LIBRARY_PATH)(geom_2d.o):	$(INCDIR)/inp.fh ./geomP.fh ./geom.fh 
$(LIBRARY_PATH)(geom_1d.o):	$(INCDIR)/inp.fh ./geomP.fh ./geom.fh 
$(LIBRARY_PATH)(geom_3d.o):	$(INCDIR)/inp.fh ./geomP.fh ./geom.fh 
$(LIBRARY_PATH)(geom_numcore.o):	./geomP.fh 
$(LIBRARY_PATH)(geom.o):	$(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh ./geomP.fh $(INCDIR)/context.fh 
