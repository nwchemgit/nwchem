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
	@rm -rf dmejiar* tblite*  ../lib/lib*tblite* ../include/tbli*  ../lib/cmake/tblite*
	@rm -rf  ../lib/pkgconfig/tbli*pc ../lib/pkgconfig/mctc-lib.pc ../lib/pkgconfig/mstore.pc ../lib/pkgconfig/toml-f.pc ../lib/pkgconfig/dftd4.pc ../lib/pkgconfig/multicharge.pc ../lib/pkgconfig/s-dftd3.pc ../lib/pkgconfig/test-drive.pc
	@rm -rf   ../lib/libmctc-lib.a ../lib/libmstore.a ../lib/libtoml-f.a ../lib/libdftd4.a ../lib/libmulticharge.a ../lib/libs-dftd3.a ../lib/libtest-drive.a
	@rm -rf   ../lib/cmake/mctc-lib ../lib/cmake/mstore ../lib/cmake/toml-f ../lib/cmake/dftd4 ../lib/cmake/multicharge ../lib/cmake/s-dftd3  ../lib/cmake/test-drive
	@rm -rf   ../share/licenses/mctc-lib ../share/licenses/mstore ../share/licenses/toml-f ../share/licenses/dftd4 ../share/licenses/multicharge ../share/licenses/s-dftd3  ../share/licenses/tblite ../share/licenses/test-drive
	@rm -rf   ../include/mctc-lib ../include/mstore ../include/toml-f ../include/dftd4 ../include/multicharge ../include/s-dftd3 ../include/test-drive ../include/dftd3.h ../include/dftd4.h ../include/s-dftd3.h
	@rm -rf   ../bin/mctc-convert ../bin/mstore-info ../bin/mstore-fortranize ../bin/dftd4 ../bin/multicharge ../bin/s-dftd3 ../bin/test-drive  ../bin/tblite
