#$Id: GNUmakefile,v 1.2 1996-03-20 22:13:29 mg141 Exp $

        OBJ = mp2_make_trhf.o \
	      mp2_make_tuhf.o \
	      mp2_grad.o \
              mp2_energy.o \
	      mp2_grad_doit.o \
	      

#ifneq ($(TARGET),CRAY-T3D)
#        OBJ += crayio.o
#endif


    LIBRARY = libmp2.a


LIB_DEFINES =  -DDEBUG_PRINT


# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */


include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(mp2_grad_utils.o):	$(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_uhf_doit.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_make_tuhf.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_energy.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mp2_rhf_doit.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_make_t.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_grad_doit.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_make_trhf.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_grad.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
