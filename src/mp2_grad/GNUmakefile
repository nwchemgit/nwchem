#$Id: GNUmakefile,v 1.21 1997-02-13 01:21:18 d3g681 Exp $

        OBJ = mp2_make_tuhf.o \
	      mp2_grad.o \
              mp2_energy.o \
	      mp2_grad_doit.o \
	      mp2_pijab_uhf.o \
	      mp2_wijab_uhf.o \
	      mp2_lai_uhf.o \
	      mp2_lai_fock_uhf.o \
	      mp2_make_cphf_rhs.o \
	      mp2_make_tot_dens.o \
	      mp2_make_tot_wdens.o \
	      mp2_make_dip_mom.o \
              mp2_pstat.o mp2_read_tiajb.o \
	      mp2_memory.o

OBJ_OPTIMIZE = mp2_back_transform.o

USES_BLAS = mp2_energy.F mp2_lai_uhf.F \
            mp2_make_tuhf.F mp2_pijab_uhf.F \
            mp2_wai_uhf.F mp2_wijab_uhf.F mp2_back_transform.F

    LIBRARY = libmp2.a


LIB_DEFINES =  -DDEBUG_PRINT


# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */


include ../config/makefile.h
include ../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(mp2_energy.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cmp2ps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mp2_make_tot_dens.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_read_tiajb.o):	$(INCDIR)/util.fh $(INCDIR)/eaf.fh 
$(LIBRARY_PATH)(mp2_pstat.o):	$(INCDIR)/util.fh $(INCDIR)/pstat.fh ./cmp2ps.fh 
$(LIBRARY_PATH)(mp2_make_tot_wdens.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_back_transform.o):	$(INCDIR)/util.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/eaf.fh ./cmp2ps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mp2_wijab_uhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/eaf.fh 
$(LIBRARY_PATH)(mp2_wai_uhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_grad.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cmp2ps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mp2_make_tuhf.o):	$(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/eaf.fh 
$(LIBRARY_PATH)(mp2_pijab_uhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/eaf.fh 
$(LIBRARY_PATH)(mp2_make_dip_mom.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mp2_make_cphf_rhs.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_lai_uhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/eaf.fh 
$(LIBRARY_PATH)(mp2_lai_fock_uhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mp2_grad_utils.o):	$(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_memory.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_grad_doit.o):	$(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/eaf.fh ./cmp2ps.fh 
