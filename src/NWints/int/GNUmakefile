# $Id: GNUmakefile,v 1.21 1996-05-20 07:16:37 d3e129 Exp $

    OBJ          = hf3PEabc.o hf3mkr.o hf3pot.o 
    OBJ_OPTIMIZE = hfmkr_9.o hfmkr_8.o hfmkr_7.o \
                   hfmkr_6.o hfmkr_5.o hfmkr_4.o \
                   hfmkr_3.o hfmkr_2.o hfmkr_1.o \
                   defNxyz.o getNxyz.o matchNxyz.o hf1.o hf1mke3.o \
                   hf1mkr.o hf1set3.o hf2.o hf2mkr.o \
                   hf2oi.o hf3OIs.o hfabc.o hfctr3.o \
                   hferi.o hfkei.o hfmke.o hfmkr.o \
                   hfnai.o hfset.o igamma.o hf_print.o \
                   hfreord_pq.o hf1_tran.o
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

   USES_BLAS     = hf2oi.f hfkei.f hfnai.f

  LIB_TARGETS = showxyz showxyz.o
  LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h

write:
	f77 -g -Nl99   -c hferi.w.f -o hferi.o
	$(MAKE)
	(cd ../..;make link)

showxyz:	showxyz.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ showxyz.o $(LIBS)

doc:	3ctrpot.tex
	latex 3ctrpot
	latex 3ctrpot
	dvips 3ctrpot -o
	ghostview 3ctrpot.ps

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(defNxyz.o):	 ./sh_order.fh 
$(LIBRARY_PATH)(getNxyz.o):	 ./sh_order.fh 
$(LIBRARY_PATH)(hf1.o):	 $(INCDIR)/pstat.fh $(INCDIR)/int_tim.fh 
$(LIBRARY_PATH)(hf1_tran.o):	 $(INCDIR)/stdio.fh 
$(LIBRARY_PATH)(hf1set3.o):	 $(INCDIR)/apiP.fh 
$(LIBRARY_PATH)(hf3pot.o):	 $(INCDIR)/stdio.fh 
$(LIBRARY_PATH)(hfmke.o):	 $(INCDIR)/pstat.fh $(INCDIR)/int_tim.fh 
$(LIBRARY_PATH)(hfmkr.o):	 $(INCDIR)/pstat.fh $(INCDIR)/int_tim.fh 
$(LIBRARY_PATH)(hfset.o):	 $(INCDIR)/apiP.fh 
$(LIBRARY_PATH)(igamma.o):	 $(INCDIR)/pstat.fh $(INCDIR)/int_tim.fh 
$(LIBRARY_PATH)(matchNxyz.o):	 ./sh_order.fh 
