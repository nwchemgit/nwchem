# $Id: GNUmakefile,v 1.5 1994-08-17 07:12:38 d3e129 Exp $

include ../../config/makefile.h

INT_OBJ = defNxyz.o getNxyz.o hf1.o hf1mke3.o \
          hf1mkr.o hf1set3.o hf2.o hf2mkr.o \
          hf2oi.o hf3OIs.o hfabc.o hfctr3.o \
          hferi.o hfkei.o hfmke.o hfmkr.o \
          hfnai.o hfset.o igamma.o hf_print.o \
          hfmkr_1.o hfmkr_2.o hfmkr_3.o \
          hfmkr_4.o hfmkr_5.o hfmkr_6.o \
          hfmkr_7.o hfmkr_8.o hfmkr_9.o

int_stamp:	$(INT_OBJ)
	/bin/rm -f ../nwints_stamp
	touch int_stamp

clean:
	/bin/rm -f $(INT_OBJ) int_stamp
