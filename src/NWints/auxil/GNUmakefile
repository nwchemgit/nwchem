# $Id: GNUmakefile,v 1.23 1996-02-14 17:14:51 d3e129 Exp $

   OBJ_OPTIMIZE = int_norm.o nmcoef.o nmcoef2c.o int_norm_2c.o int_chk_sum.o 
            OBJ = do_chk_sum.o intb_chk_sum.o xlmcoeff.o int_spcart.o xlm_norm.o
LIBRARY = libnwints.a
LIB_TARGETS = spherical spherical.o testnat testnat.o

USES_BLAS = do_chk_sum.F

include ../../config/makefile.h
include ../../config/makelib.h

spherical:	spherical.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ spherical.o $(LIBS)
testnat:	testnat.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ testnat.o

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(do_chk_sum.o):	 $(INCDIR)/util.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_chk_sum.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_norm.o):	 $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/bas_exndcf_sfn.fh $(INCDIR)/bas_exndcf_dec.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_norm_2c.o):	 $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/bas_exndcf_sfn.fh $(INCDIR)/bas_exndcf_dec.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_spcart.o):	 $(INCDIR)/stdio.fh ./spcartP.fh $(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(intb_chk_sum.o):	 $(INCDIR)/sym.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(spherical.o):	 $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(xlm_norm.o):	 $(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(xlmcoeff.o):	 $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh 
