#$Id$
# makefile,v 1.8 1994/12/05 20:37:08 og845 Exp

# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../../config/makefile.h

install/lib/libxc.a:
	./build_libxc.sh $(LIBXC_VERSION)

LIB_TARGETS += libxc

include ../../config/makelib.h

clean:
	@rm -rf  libxc* install 



