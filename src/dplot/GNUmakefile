#$Id: GNUmakefile,v 1.4 1996-10-01 21:02:59 d3e129 Exp $
# makefile,v 1.8 1994/12/05 20:37:08 og845 Exp

   OBJ = dplot.o get_grid.o create_contour.o dmat.o errmem.o \
         get_dens.o rhogen.o dplot_input.o

   LIBRARY = libdplot.a
 USES_BLAS = dplot.F dmat.F create_contour.F rhogen.F


include ../config/makefile.h
include ../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(rhogen.o):	$(INCDIR)/bas.fh 
$(LIBRARY_PATH)(get_grid.o):	$(INCDIR)/stdio.fh 
$(LIBRARY_PATH)(get_dens.o):	$(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dplot_input.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(dplot.o):	$(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dmat.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(create_contour.o):	$(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
