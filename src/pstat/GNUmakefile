# $Id: GNUmakefile,v 1.4 1994-08-22 03:11:00 d3g681 Exp $


#LIBRARY = libpstat.a
LIBRARY =      libutil.a

OBJ =		pstat_init.o	pstat_terminate.o	pstat_allocate.o \
		pstat_on.o	pstat_off.o		pstat_print_all.o pstat_free.o

HEADERS =	pstat.fh	pstat_consts.fh


include ../config/makefile.h
include ../config/makelib.h

$(OBJ): pstatP.fh pstat_consts.fh pstat.fh
