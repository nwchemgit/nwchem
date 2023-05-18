#
# $Id$
#
#USE_I4FLAGS=1
include ../../config/makefile.h

ifdef NO_SIFS
OBJ = sifs_stubs.o
else

OBJ_OPTIMIZE = colib1.o colib2.o colib3.o colib4.o \
	       colib5.o colib6.o colib7.o colib8.o \
	       colib8a.o colib9.o blaswrapper.o izero.o \
	       cdcopy.o 

endif
      LIBRARY = libnwints.a


     LIB_INCLUDES = -I../../nwdft/grid -I../../nwdft/include \
	 -I../../ddscf

      LIB_DEFINES = 


include ../../config/makelib.h

izero.F:
	./build_colib.sh

blaswrapper.f cdcopy.f\
colib1.f colib2.f colib3.f colib4.f \
colib5.f colib6.f colib7.f colib8.f \
colib8a.f colib9.f: izero.F
	./build_colib.sh

MAKEFLAGS = -j 1 --no-print-directory
ALLF = $(filter-out ./sifs_stubs.F, $(shell find . -name '*.F'))


#bummer.F:
#	./build_colib.sh

cleanF:
	@rm -rf colu* $(ALLF)  dependencies sifs.patched

