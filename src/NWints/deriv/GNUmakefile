# $Id: GNUmakefile,v 1.13 1997-01-21 07:40:56 d3e129 Exp $

      OBJ_OPTIMIZE = hf1d.o      hf2d.o     hfd2oi.o \
                     hfderi.o    hfdkei.o   hfdnai.o \
                     hfefi.o     hf2d_2a.o  hf2d_2b.o \
                     hf2d_2a2b.o hf2d_3a.o  hf2d_3a2b.o \
                     hf1PRtoAB.o hf2PRtoAB.o hf2QStoCD.o hfd3oi.o

     LIBRARY = libnwints.a

USES_BLAS = hf1PRtoAB.f hf1d.F hf2PRtoAB.f hf2QStoCD.f hf2d.F hf2d_2a.F hf2d_2a2b.F \
            hf2d_2b.F hf2d_3a.F hf2d_3a2b.F hfd3oi.F

include ../../config/makefile.h
include ../../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(hf2d.o):	$(INCDIR)/global.fh $(INCDIR)/apiP.fh 
