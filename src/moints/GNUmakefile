#$Id: GNUmakefile,v 1.57 1996-12-09 20:13:14 gg502 Exp $

        OBJ = atwtest.o \
	      direct_mp2.o \
	      mointdriver.o \
	      moints_1e.o \
	      moints_dbg.o \
	      moints_gsh1.o \
	      moints_mem.o \
	      moints_mp2.o \
	      moints_semi.o \
	      moints_stats.o \
              moints_full.o \
              sifs.o

ifneq ($(TARGET),CRAY-T3D)
        OBJ += crayio.o
endif

OBJ_OPTIMIZE = \
	      moints_2x.o \
              moints_6x.o \
              moints_aux2.o 

    LIBRARY = libmoints.a

USES_BLAS = moints_2x.F moints_6x.F moints_aux2.F moints_ccsd.F moints_full.F moints_semi.F sifs.F

LIB_DEFINES =  -DDEBUG_PRINT

    HEADERS = moints_stats.fh

#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#	       -DNOCOMMS              /* Turn off ALL communication calls */
#              -DBLOCK_TRANSF         /* Use block-transfer for comms calls */
#

include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(moints_stats.o):	$(INCDIR)/msgids.fh ./moints_stats.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(moints_semi.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/schwarz.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(moints_mp2.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(moints_mem.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(moints_gsh1.o):	$(INCDIR)/util.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(moints_full.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cmointsbuf.fh $(INCDIR)/bitops.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(moints_dbg.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(moints_ccsd.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/schwarz.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(moints_aux2.o):	$(INCDIR)/schwarz.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bitops.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(moints_6x.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/schwarz.fh $(INCDIR)/msgids.fh ./moints_stats.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(moints_2x.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/schwarz.fh $(INCDIR)/msgids.fh ./moints_stats.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(moints_1e.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mointdriver.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(get_JKelt.o):	$(INCDIR)/global.fh 
$(LIBRARY_PATH)(direct_mp2.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(atwtest.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
