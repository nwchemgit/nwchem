# $Id: GNUmakefile,v 1.2 1994-07-29 22:42:55 d3e129 Exp $

include ../../config/makefile.h
DERIV_OBJ = hf1PRtoAB.o hf1d.o   hf2PRtoAB.o \
            hf2QStoCD.o hf2d.o   hfd2oi.o \
            hfderi.o    hfdkei.o hfdnai.o \
            hfefi.o

derivobj:	$(DERIV_OBJ)

#new:	
#	cp hf2d.F.new hf2d.F
#	$(MAKE) debug
#
#old:	
#	cp hf2d.F.org hf2d.F
#	$(MAKE) debug
#
#debug:	
#	$(MAKE) DEFINES="-DSUN -DHF2D_DEBUG"

clean:
	/bin/rm -f $(DERIV_OBJ)

