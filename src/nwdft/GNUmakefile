#$Id: GNUmakefile,v 1.16 1995-11-19 19:42:01 d3g681 Exp $
# makefile,v 1.8 1994/12/05 20:37:08 og845 Exp

   OBJ = nwdft.o

   LIBRARY = libnwdft.a
   LIB_DEFINES = -DDFT_TIMER
   SUBDIRS = input_dft scf_dft coulomb xc util 3d
LIB_INCLUDES = -I./include


include ../config/makefile.h
include ../config/makelib.h

#

nwdft.o:  \
      include/timer.fh   \
      include/quad.fh  \
      include/cdft.fh  \
      include/dimensions.fh  

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(nwdft.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/cdft.fh $(INCDIR)/bas.fh 
