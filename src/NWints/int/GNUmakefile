# $Id: GNUmakefile,v 1.14 1995-10-17 08:52:36 d3g681 Exp $

    OBJ_OPTIMIZE = hfmkr_9.o hfmkr_8.o hfmkr_7.o \
                   hfmkr_6.o hfmkr_5.o hfmkr_4.o \
                   hfmkr_3.o hfmkr_2.o hfmkr_1.o \
                   defNxyz.o getNxyz.o matchNxyz.o hf1.o hf1mke3.o \
                   hf1mkr.o hf1set3.o hf2.o hf2mkr.o \
                   hf2oi.o hf3OIs.o hfabc.o hfctr3.o \
                   hferi.o hfkei.o hfmke.o hfmkr.o \
                   hfnai.o hfset.o igamma.o hf_print.o \
                   hfreord_pq.o

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
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
defNxyz.o :	 sh_order.fh 
getNxyz.o :	 sh_order.fh 
hf1.o :	 $(INCDIR)/pstat.fh $(INCDIR)/int_tim.fh 
hf1set3.o :	 $(INCDIR)/apiP.fh 
hfmke.o :	 $(INCDIR)/pstat.fh $(INCDIR)/int_tim.fh 
hfmkr.o :	 $(INCDIR)/pstat.fh $(INCDIR)/int_tim.fh 
hfset.o :	 $(INCDIR)/apiP.fh 
igamma.o :	 $(INCDIR)/pstat.fh $(INCDIR)/int_tim.fh 
matchNxyz.o :	 sh_order.fh 
