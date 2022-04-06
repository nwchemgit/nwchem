#$Id$
# makefile,v 1.8 1994/12/05 20:37:08 og845 Exp

# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../../config/makefile.h

../lib/libnwc_tblite.a:
ifdef TBLITE_MESON
	./build_tblite.sh 
else
	./build_tblite_cmake.sh
endif

LIB_TARGETS += tblite


include ../../config/makelib.h

clean:
	@rm -rf tblite*  ../lib/lib*tblite*



