# $Id: GNUmakefile,v 1.12 1995-12-13 01:36:28 d3g681 Exp $


#LIBRARY = libpstat.a
LIBRARY =      libutil.a

OBJ =		pstat_init.o	pstat_term.o	pstat_alloc.o \
		pstat_on.o	pstat_off.o		pstat_pr_all.o \
                pstat_free.o	pstat_pr_han.o	pstat_pr_det.o

HEADERS =	pstat.fh	pstat_consts.fh


include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(pstat_alloc.o):	 ./pstatP.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(pstat_free.o):	 ./pstatP.fh 
$(LIBRARY_PATH)(pstat_init.o):	 ./pstatP.fh 
$(LIBRARY_PATH)(pstat_off.o):	 $(INCDIR)/util.fh ./pstatP.fh 
$(LIBRARY_PATH)(pstat_on.o):	 $(INCDIR)/util.fh ./pstatP.fh 
$(LIBRARY_PATH)(pstat_pr_all.o):	 $(INCDIR)/util.fh ./pstatP.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(pstat_pr_det.o):	 ./pstatP.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(pstat_pr_han.o):	 ./pstatP.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(pstat_term.o):	 ./pstatP.fh 
