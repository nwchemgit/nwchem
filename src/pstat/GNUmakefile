# $Id: GNUmakefile,v 1.7 1995-10-17 08:55:40 d3g681 Exp $


#LIBRARY = libpstat.a
LIBRARY =      libutil.a

OBJ =		pstat_init.o	pstat_terminate.o	pstat_allocate.o \
		pstat_on.o	pstat_off.o		pstat_print_all.o \
                pstat_free.o	pstat_print_handles.o	pstat_print_detail.o

HEADERS =	pstat.fh	pstat_consts.fh


include ../config/makefile.h
include ../config/makelib.h

$(OBJ): pstatP.fh pstat_consts.fh pstat.fh
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
pstat_allocate.o :	 pstatP.fh 
pstat_free.o :	 pstatP.fh 
pstat_init.o :	 pstatP.fh 
pstat_off.o :	 $(INCDIR)/util.fh pstatP.fh 
pstat_on.o :	 $(INCDIR)/util.fh pstatP.fh 
pstat_print_all.o :	 $(INCDIR)/util.fh pstatP.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
pstat_print_detail.o :	 pstatP.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
pstat_print_handles.o :	 pstatP.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
pstat_terminate.o :	 pstatP.fh 
