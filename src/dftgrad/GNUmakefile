   OBJ = dft_grad_top.o  dftg_force.o \
         dftg_getxc.o       dftg_grad2.o \
         grad_dft.o      dftg_quadv0.o     \
         dftg_getcden.o  dftg_getxden.o  \
         dftg_bal_sum.o dftg_hf_x.o dftg_beckew.o dftg_wderiv.o
   OBJ_OPTIMIZE = dftg_cdfit.o  dftg_quadv0a.o  dftg_quadv0b.o  
ifeq ($(NWCHEM_TARGET),SGITFP)
   OBJ = dftg_cdfit.o  dft_grad_top.o  dftg_force.o \
         dftg_getxc.o       dftg_grad2.o \
         grad_dft.o      dftg_quadv0.o     \
         dftg_get_cden.o  dftg_get_xden.o\
         dftg_bal_sum.o dftg_hf_x.o dftg_quadv0a.o dftg_beckew.o dftg_wderiv.o
   OBJ_OPTIMIZE = dftg_quadv0b.o 
endif
   USES_BLAS = dftg_force.F dftg_quadv0a.F dftg_cdfit.F dftg_quadv0b.F
LIB_INCLUDES = -I../nwdft/include

   LIBRARY = libdftgrad.a

include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(dftg_bal_sum.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dftg_beckew.o):	 $(INCDIR)/bas_staticP.fh 
$(LIBRARY_PATH)(dftg_cdfit.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dftg_force.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dftg_getxc.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dftg_grad2.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dftg_hf_x.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dftg_quadv0.o):	 ../nwdft/include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh ../nwdft/include/quad.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dftg_quadv0a.o):	 ../nwdft/include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh ../nwdft/include/quad.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dftg_quadv0b.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dftg_wderiv.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(grad_dft.o):	 $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(grad_xc.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
