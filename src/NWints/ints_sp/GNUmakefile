# $Id: GNUmakefile,v 1.5 1994-08-22 00:02:54 d3g681 Exp $


      OBJ = genr70.o init70.o qout70.o \
            sinfo.o sp0000.o sp0001.o \
            sp0011.o sp0101.o sp0111.o \
            sp1111.o tabgen.o 

    LIBRARY = libnwints.a


  USES_BLAS = sp1111.F

include ../../config/makefile.h
include ../../config/makelib.h

$(OBJ):	auxvar.fh geom.fh maxc.fh  pgeom.fh qgeom.fh \
	const.fh  ginf.fh miscg.fh picon.fh shllfo.fh
