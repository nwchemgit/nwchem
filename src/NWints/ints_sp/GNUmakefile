# $Id: GNUmakefile,v 1.4 1994-08-21 22:00:36 d3e129 Exp $

include ../../config/makefile.h

INTSP_OBJ = genr70.o init70.o qout70.o \
            sinfo.o sp0000.o sp0001.o \
            sp0011.o sp0101.o sp0111.o \
            sp1111.o tabgen.o 


int_sp_obj:	$(INTSP_OBJ)

clean:
	/bin/rm -f $(INTSP_OBJ) 

$(INTSP_OBJ):	auxvar.fh geom.fh maxc.fh  pgeom.fh qgeom.fh \
	const.fh  ginf.fh miscg.fh picon.fh shllfo.fh
