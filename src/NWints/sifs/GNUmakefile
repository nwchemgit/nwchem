#
# $Id$
#
ifdef NO_SIFS
OBJ = sifs_stubs.o
else
OBJ_OPTIMIZE = sifrd1.o sifr2.o sifr1x.o sifr1n.o sifpre.o \
sifo2f.o sifmge.o siffw2.o siffw1.o siffr2.o siffr1.o sifew2.o \
sifew1.o sife2.o sife1.o sifd2.o sifd1.o sifcfg.o sifc2f.o \
sif2w8.o sifzwh.o sifzh2.o sifwh.o sifw2.o sifw1x.o siftyp.o \
siftdy.o sifstv.o sifskh.o sifsk1.o sifsce.o sifrsh.o sifrh2.o \
sifrh1.o sifrd2.o

OBJ = bummer.o
endif
      LIBRARY = libnwints.a
ifndef NO_SIFS
    USES_BLAS = sifd1.F sifd2.F sife1.F sife2.F
endif


  LIB_DEFINES = 


include ../../config/makefile.h
include ../../config/makelib.h

bummer.F sif2w8.F sifc2f.F sifcfg.F \
sifd1.F sifd2.F sife1.F sife2.F \
sifew1.F sifew2.F siffr1.F siffr2.F \
siffw1.F siffw2.F sifmge.F sifo2f.F \
sifpre.F sifr1n.F sifr1x.F sifr2.F \
sifrd1.F sifrd2.F sifrh1.F sifrh2.F sifrsh.F \
sifsce.F sifsk1.F sifskh.F sifstv.F \
siftdy.F siftyp.F sifw1x.F sifw2.F \
sifwh.F sifzh2.F sifzwh.F:
	./build_sifs.sh

  MAKEFLAGS = -j 1 --no-print-directory

cleanF:
	@rm -rf colu* *F dependencies

