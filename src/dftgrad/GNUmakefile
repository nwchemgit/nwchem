   OBJ = dft_grad_top.o  dftg_force.o \
         dftg_getxc.o       dftg_grad2.o \
         grad_dft.o      dftg_quadv0.o     \
         dftg_get_cden.o  dftg_get_xden.o  \
         dftg_bal_sum.o dftg_hf_x.o
   OBJ_OPTIMIZE = dftg_cdfit.o  dftg_quadv0_a.o  dftg_quadv0_b.o  
ifeq ($(NWCHEM_TARGET),SGITFP)
   OBJ = dftg_cdfit.o  dft_grad_top.o  dftg_force.o \
         dftg_getxc.o       dftg_grad2.o \
         grad_dft.o      dftg_quadv0.o     \
         dftg_get_cden.o  dftg_get_xden.o\
         dftg_bal_sum.o dftg_hf_x.o dftg_quadv0_a.o 
   OBJ_OPTIMIZE = dftg_quadv0_b.o 
endif
   USES_BLAS = dftg_force.F dftg_quadv0_a.F dftg_cdfit.F dftg_quadv0_b.F

   LIBRARY = libdftgrad.a

include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
dft_grad_top.o :	 $(INCDIR)/quad.fh 
dftg_bal_sum.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/global.fh 
dftg_cdfit.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
dftg_force.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
dftg_getxc.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
dftg_grad2.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
dftg_hf_x.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
dftg_quadv0.o :	 $(INCDIR)/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/quad.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
dftg_quadv0_b.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
grad_dft.o :	 $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
grad_xc.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
# DO NOT DELETE THIS LINE -- make depend depends on it.
