# $Id: GNUmakefile,v 1.11 1995-11-27 02:18:25 d3g681 Exp $


#LIBRARY = libpstat.a
LIBRARY =      libutil.a

OBJ =		pstat_init.o	pstat_terminate.o	pstat_allocate.o \
		pstat_on.o	pstat_off.o		pstat_print_all.o \
                pstat_free.o	pstat_print_handles.o	pstat_print_detail.o

HEADERS =	pstat.fh	pstat_consts.fh


include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(pstat_allocate.o):	 ./pstatP.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(pstat_free.o):	 ./pstatP.fh 
$(LIBRARY_PATH)(pstat_init.o):	 ./pstatP.fh 
$(LIBRARY_PATH)(pstat_off.o):	 $(INCDIR)/util.fh ./pstatP.fh 
$(LIBRARY_PATH)(pstat_on.o):	 $(INCDIR)/util.fh ./pstatP.fh 
$(LIBRARY_PATH)(pstat_print_all.o):	 $(INCDIR)/util.fh ./pstatP.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(pstat_print_detail.o):	 ./pstatP.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(pstat_print_handles.o):	 ./pstatP.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(pstat_terminate.o):	 ./pstatP.fh 
