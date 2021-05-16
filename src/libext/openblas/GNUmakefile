#$Id$
# makefile,v 1.8 1994/12/05 20:37:08 og845 Exp

# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../../config/makefile.h

../lib/libnwc_openblas.a:
	./build_openblas.sh 
LIB_TARGETS += OpenBLAS

include ../../config/makelib.h

clean:
	@rm -rf OpenBLAS*  ../lib/libnwc_openbla*



