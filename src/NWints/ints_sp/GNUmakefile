# $Id: GNUmakefile,v 1.8 1995-11-16 20:38:04 d3g681 Exp $


      OBJ_OPTIMIZE = genr70.o init70.o qout70.o \
                     sinfo.o sp0000.o sp0001.o \
                     sp0011.o sp0101.o sp0111.o \
                     sp1111.o tabgen.o 

    LIBRARY = libnwints.a


  USES_BLAS = sp1111.F

include ../../config/makefile.h
include ../../config/makelib.h

$(OBJ):	auxvar.fh geom.fh maxc.fh  pgeom.fh qgeom.fh \
	const.fh  ginf.fh miscg.fh picon.fh shllfo.fh


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
genr70.o :	 shllfo.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
init70.o :	 picon.fh maxc.fh auxvar.fh 
sinfo.o :	 shllfo.fh qgeom.fh picon.fh pgeom.fh miscg.fh maxc.fh ginf.fh geom.fh const.fh 
sp0000.o :	 shllfo.fh qgeom.fh pgeom.fh miscg.fh maxc.fh ginf.fh geom.fh auxvar.fh 
sp0001.o :	 shllfo.fh qgeom.fh pgeom.fh miscg.fh maxc.fh ginf.fh geom.fh auxvar.fh 
sp0011.o :	 shllfo.fh qgeom.fh pgeom.fh miscg.fh maxc.fh ginf.fh geom.fh auxvar.fh 
sp0101.o :	 shllfo.fh qgeom.fh pgeom.fh miscg.fh maxc.fh ginf.fh geom.fh const.fh auxvar.fh 
sp0111.o :	 shllfo.fh qgeom.fh pgeom.fh miscg.fh maxc.fh ginf.fh geom.fh const.fh auxvar.fh 
sp1111.o :	 shllfo.fh qgeom.fh pgeom.fh miscg.fh maxc.fh ginf.fh geom.fh const.fh auxvar.fh 
