# $Id: GNUmakefile,v 1.12 1996-01-23 00:26:41 d3e129 Exp $


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
$(LIBRARY_PATH)(genr70.o):	 ./shllfo.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/bas_ibs_sfn.fh $(INCDIR)/bas_ibs_dec.fh $(INCDIR)/bas_exndcf_sfn.fh $(INCDIR)/bas_exndcf_dec.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(init70.o):	 ./picon.fh ./maxc.fh ./auxvar.fh 
$(LIBRARY_PATH)(sinfo.o):	 ./spgeom.fh ./shllfo.fh ./qgeom.fh ./picon.fh ./pgeom.fh ./miscg.fh ./maxc.fh ./ginf.fh ./const.fh 
$(LIBRARY_PATH)(sp0000.o):	 ./spgeom.fh ./shllfo.fh ./qgeom.fh ./pgeom.fh ./miscg.fh ./maxc.fh ./ginf.fh ./auxvar.fh 
$(LIBRARY_PATH)(sp0001.o):	 ./spgeom.fh ./shllfo.fh ./qgeom.fh ./pgeom.fh ./miscg.fh ./maxc.fh ./ginf.fh ./auxvar.fh 
$(LIBRARY_PATH)(sp0011.o):	 ./spgeom.fh ./shllfo.fh ./qgeom.fh ./pgeom.fh ./miscg.fh ./maxc.fh ./ginf.fh ./auxvar.fh 
$(LIBRARY_PATH)(sp0101.o):	 ./spgeom.fh ./shllfo.fh ./qgeom.fh ./pgeom.fh ./miscg.fh ./maxc.fh ./ginf.fh ./const.fh ./auxvar.fh 
$(LIBRARY_PATH)(sp0111.o):	 ./spgeom.fh ./shllfo.fh ./qgeom.fh ./pgeom.fh ./miscg.fh ./maxc.fh ./ginf.fh ./const.fh ./auxvar.fh 
$(LIBRARY_PATH)(sp1111.o):	 ./spgeom.fh ./shllfo.fh ./qgeom.fh ./pgeom.fh ./miscg.fh ./maxc.fh ./ginf.fh ./const.fh ./auxvar.fh 
