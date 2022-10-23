
# $Id

USE_I4FLAGS=1
include ../../config/makefile.h

OBJ_OPTIMIZE = \
 sif2w8.o sifc2f.o sifcfg.o sifd1.o sifd2.o sife1.o sife2.o \
 sifew1.o sifew2.o siffr1.o siffr2.o siffw1.o siffw2.o sifmge.o \
 sifo2f.o sifpre.o sifr1n.o sifr1x.o sifr2.o sifrd1.o sifrd2.o \
 sifrh1.o sifrh2.o sifrsh.o sifsce.o sifsk1.o sifskh.o sifstv.o \
 siftdy.o siftyp.o sifw1x.o sifw2.o sifwh.o sifzh2.o sifzwh.o \
 bummer.o seqwbf.o plab1.o plab8.o plab16.o plab32.o cdcopy.o \
 ulab16.o ulab1.o ulab32.o ulab8.o aiopen.o izero.o aiwait.o  aiwrit.o \
 seqrbf.o airead.o

LIBRARY = libnwints.a

LIB_INCLUDES = -I../../nwdft/grid -I../../nwdft/include \
 -I../../ddscf

include ../../config/makelib.h

iwfmt.x:  $(OJB_OPTIMIZE) iwfmt.o
	$(LINK.f) iwfmt.o $(LIBS) -o iwmft.x
