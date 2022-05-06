#$Id$
# makefile,v 1.8 1994/12/05 20:37:08 og845 Exp

# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../../config/makefile.h


../lib/libplumed.a:
	./build_plumed.sh 
LIB_TARGETS += scalapack

include ../../config/makelib.h

clean:
	@rm -rf  plumed* ../lib/libplu* ../bin/plumed* ../include/plumed*



