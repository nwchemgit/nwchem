#$Id: GNUmakefile,v 1.10 1996-09-28 20:17:12 d3g681 Exp $

        OBJ = mp2_make_trhf.o \
	      mp2_make_tuhf.o \
	      mp2_grad.o \
              mp2_energy.o \
	      mp2_grad_doit.o \
	      mp2_read_tuhf.o \
	      mp2_read_trhf.o \
	      mp2_pijab_uhf.o \
	      mp2_wijab_uhf.o \
	      mp2_lai_uhf.o \
	      mp2_vvoo_uhf.o \
	      mp2_wai_uhf.o \
	      mp2_lai_fock_uhf.o \
	      mp2_make_cphf_rhs.o \
	      mp2_make_tot_dens.o \
	      mp2_make_tot_wdens.o \
	      mp2_make_dip_mom.o \
              mp2_back_transform.o

USES_BLAS = mp2_energy.F mp2_lai_uhf.F mp2_make_t.F mp2_make_trhf.F \
            mp2_make_tuhf.F mp2_pijab_uhf.F mp2_vooo_uhf.F \
            mp2_vvoo_uhf.F mp2_vvvo_uhf.F mp2_wai_uhf.F mp2_wijab_uhf.F

    LIBRARY = libmp2.a


LIB_DEFINES =  -DDEBUG_PRINT


# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */


include ../config/makefile.h
include ../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(mp2_make_tot_wdens.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_back_transform.o):	$(INCDIR)/util.fh $(INCDIR)/schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/eaf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mp2_wijab_uhf.o):	$(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_wai_uhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_vvoo_uhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_uhf_doit.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_rhf_doit.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_read_tuhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_read_trhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_pijab_uhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_make_tuhf.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_make_trhf.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_make_tot_dens.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_make_t.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_make_dip_mom.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_make_cphf_rhs.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_lai_uhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_lai_fock_uhf.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mp2_grad_utils.o):	$(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_grad_doit.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mp2_grad.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(mp2_energy.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
