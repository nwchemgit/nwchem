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

bummer.F \
aiopen.F aiwait.F airead.F aiwrit.F \
seqwbf.F seqrbf.F \
plab1.F plab8.F plab16.F plab32.F \
ulab1.F ulab8.F ulab16.F ulab32.F \
sif2w8.F sifcfg.F \
sifd1.F sifd2.F sife1.F sife2.F \
sifew1.F sifew2.F \
siffw1.F siffw2.F sifo2f.F \
sifr1n.F sifr2.F \
sifrd1.F sifrd2.F sifrh1.F sifrh2.F \
sifsk1.F sifskh.F \
sifw2.F izero.F\
sifwh.F sifzh2.F sifzwh.F:
	./build_sifs.sh

  MAKEFLAGS = -j 1 --no-print-directory

cleanF:
	@rm -rf colu* *F dependencies

