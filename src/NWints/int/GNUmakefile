# $Id$

    OBJ          = 
    OBJ_OPTIMIZE = hf3PEabc.o hf3mkr.o hf3pot.o \
                   defNxyz.o getNxyz.o matchNxyz.o hf1.o hf1mke3.o \
                   hf1mkr.o hf1set3.o hf2.o hf2mkr.o \
                   hf2oi.o hf3OIs.o hfabc.o hfctr3.o \
                   hferi.o hferi_gen.o hfkei.o hfmke.o hfmkr.o \
                   hfnai.o hfset.o igamma.o hf_print.o \
                   hfreord_pq.o hf1_tran.o hf1_er.o hfnai_er.o case_impl.o

  USES_BLAS = hf2oi.f hfkei.f hfnai.f hfnai_er.F hf2.F igamma.F hferi.F hf1.F hf1mke3.f

  LIB_TARGETS = showxyz showxyz.o int_order.o int_order
  LIBRARY = libnwints.a

  HEADERS= sh_order.fh case.fh

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
	$(FC) $(FFLAGS) $(LDFLAGS) -o ${BINDIR}/int_order int_order.o $(LIBS) -lga -lutil
	rm -f int_order.o
int_order.o:	int_order.F	int_order.fh

doc:	3ctrpot.tex
	latex 3ctrpot
	latex 3ctrpot
	dvips 3ctrpot -o
	ghostview 3ctrpot.ps

