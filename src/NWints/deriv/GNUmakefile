# $Id: GNUmakefile,v 1.3 1994-08-17 07:12:36 d3e129 Exp $

include ../../config/makefile.h
DERIV_OBJ = hf1PRtoAB.o hf1d.o   hf2PRtoAB.o \
            hf2QStoCD.o hf2d.o   hfd2oi.o \
            hfderi.o    hfdkei.o hfdnai.o \
            hfefi.o

deriv_stamp:	$(DERIV_OBJ)
	/bin/rm -f ../nwints_stamp
	touch deriv_stamp 

clean:
	/bin/rm -f $(DERIV_OBJ) deriv_stamp 

