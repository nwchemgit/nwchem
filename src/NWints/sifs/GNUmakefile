#
# $Id$
#
USE_I4FLAGS=1
include ../../config/makefile.h

ifdef NO_SIFS
OBJ = sifs_stubs.o
else
OBJ_OPTIMIZE = sifrd1.o sifr2.o sifr1x.o sifr1n.o  \
sifo2f.o siffw2.o siffw1.o sifew2.o \
sifew1.o sife2.o sife1.o sifd2.o sifd1.o sifcfg.o \
sif2w8.o sifzwh.o sifzh2.o sifwh.o sifw2.o \
sifskh.o sifsk1.o sifrh2.o \
sifrh1.o sifrd2.o aiopen.o aiwait.o airead.o \
aiwrit.o seqwbf.o seqrbf.o plab1.o plab8.o plab16.o plab32.o \
ulab1.o ulab8.o ulab16.o ulab32.o cdcopy.o izero.o

OBJ = bummer.o
endif
      LIBRARY = libnwints.a


     LIB_INCLUDES = -I../../nwdft/grid -I../../nwdft/include \
	 -I../../ddscf

      LIB_DEFINES = 


include ../../config/makelib.h

MAKEFLAGS = -j 1 --no-print-directory
ALLF = $(filter-out ./sifs_stubs.F, $(shell find . -name '*.F'))

bummer.F:
	./build_sifs.sh

cleanF:
	@rm -rf colu* $(ALLF)  dependencies sifs.patched

