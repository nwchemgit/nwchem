# $Id: GNUmakefile,v 1.30 1996-10-23 05:48:23 d3g681 Exp $

    OBJ          = 
    OBJ_OPTIMIZE = hf3PEabc.o hf3mkr.o hf3pot.o \
                   defNxyz.o getNxyz.o matchNxyz.o hf1.o hf1mke3.o \
                   hf1mkr.o hf1set3.o hf2.o hf2mkr.o \
                   hf2oi.o hf3OIs.o hfabc.o hfctr3.o \
                   hferi.o hfkei.o hfmke.o hfmkr.o \
                   hfnai.o hfset.o igamma.o hf_print.o \
                   hfreord_pq.o hf1_tran.o
# rjh
# These take forever to compile, crash the SGI 7.0 compiler and are unused.
# hfmkr_9.o hfmkr_8.o hfmkr_7.o 
# hfmkr_6.o hfmkr_5.o hfmkr_4.o
# hfmkr_3.o hfmkr_2.o hfmkr_1.o 

# debug OBJ setting
#    OBJ          = hf3PEabc.o hf3mkr.o hf3pot.o hfmkr_9.o \
#                   hfmkr_8.o hfmkr_7.o \
#                   hfmkr_6.o hfmkr_5.o hfmkr_4.o \
#                   hfmkr_3.o hfmkr_2.o hfmkr_1.o \
#                   defNxyz.o getNxyz.o matchNxyz.o hf1.o hf1mke3.o \
#                   hf1mkr.o hf1set3.o hf2.o hf2mkr.o \
#                   hf2oi.o hf3OIs.o hfabc.o hfctr3.o \
#                   hferi.o hfkei.o hfmke.o hfmkr.o \
#                   hfnai.o hfset.o igamma.o hf_print.o \
#                   hfreord_pq.o hf1_tran.o

  USES_BLAS = hf2oi.f hfkei.f hfnai.f

  LIB_TARGETS = showxyz showxyz.o int_order.o int_order
  LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h

BIN_OTHER =	${BINDIR}/showxyz	${BINDIR}/int_order
bins:	$(BIN_OTHER)


write:
	f77 -g -Nl99   -c hferi.w.f -o hferi.o
	$(MAKE)
	(cd ../..;make link)

showxyz ${BINDIR}/showxyz:	showxyz.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o ${BINDIR}/showxyz showxyz.o $(LIBS)
	rm -f showxyz.o
int_order ${BINDIR}/int_order:	int_order.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o ${BINDIR}/int_order int_order.o $(LIBS) -lutil
	rm -f int_order.o
int_order.o:	int_order.F	int_order.fh

doc:	3ctrpot.tex
	latex 3ctrpot
	latex 3ctrpot
	dvips 3ctrpot -o
	ghostview 3ctrpot.ps

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(int_order.o):	./int_order.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(hfset.o):	$(INCDIR)/apiP.fh 
$(LIBRARY_PATH)(hf3pot.o):	$(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(hf1set3.o):	$(INCDIR)/apiP.fh 
$(LIBRARY_PATH)(hf1_tran.o):	$(INCDIR)/stdio.fh 
$(LIBRARY_PATH)(hf1.o):	$(INCDIR)/mafdecls.fh 
