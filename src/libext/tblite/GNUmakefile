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
	@rm -rf tblite*  ../lib/lib*tblite* ../include/tbli*  ../lib/cmake/tblite*
	@rm -rf  ../lib/pkgconfig/tbli*pc ../lib/pkgconfig/mctc-lib.pc ../lib/pkgconfig/mstore.pc ../lib/pkgconfig/toml-f.pc ../lib/pkgconfig/dftd4.pc ../lib/pkgconfig/multicharge.pc ../lib/pkgconfig/s-dftd3.pc
	@rm -rf   ../lib/libmctc-lib.a ../lib/libmstore.a ../lib/libtoml-f.a ../lib/libdftd4.a ../lib/libmulticharge.a ../lib/libs-dftd3.a
	@rm -rf   ../lib/cmake/mctc-lib ../lib/cmake/mstore ../lib/cmake/toml-f ../lib/cmake/dftd4 ../lib/cmake/multicharge ../lib/cmake/s-dftd3
