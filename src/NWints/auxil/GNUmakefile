# $Id: GNUmakefile,v 1.4 1994-08-17 07:12:34 d3e129 Exp $

include ../../config/makefile.h

AUX_OBJ = int_normalize.o nmcoef.o int_chk_sum.o

aux_stamp:	$(AUX_OBJ)
	/bin/rm -f ../nwints_stamp
	touch aux_stamp 

clean:
	/bin/rm -f $(AUX_OBJ) aux_stamp 
