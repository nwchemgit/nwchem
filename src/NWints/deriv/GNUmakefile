# $Id: GNUmakefile,v 1.9 1995-10-17 08:52:34 d3g681 Exp $

      OBJ_OPTIMIZE = hf1d.o      hf2d.o     hfd2oi.o \
                     hfderi.o    hfdkei.o   hfdnai.o \
                     hfefi.o     hf2d_2a.o  hf2d_2b.o \
                     hf2d_2a2b.o hf2d_3a.o  hf2d_3a2b.o \
                     hf1PRtoAB.o hf2PRtoAB.o hf2QStoCD.o 

     LIBRARY = libnwints.a

   USES_BLAS = hf1d.F hf2d.F hf2d_2a.F hf2d_2b.F hf2d_2a2b.F hf2d_3a.F hf2d_3a2b.F \
   hf1PRtoAB.f  hf2PRtoAB.f  hf2QStoCD.f

include ../../config/makefile.h
include ../../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
hf1d.o :	 $(INCDIR)/pstat.fh $(INCDIR)/int_tim.fh 
hf2d.o :	 $(INCDIR)/global.fh $(INCDIR)/apiP.fh 
