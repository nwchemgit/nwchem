# $Id: GNUmakefile,v 1.3 1994-08-20 00:26:10 d3g681 Exp $


OBJ =		pstat_init.o	pstat_terminate.o	pstat_allocate.o \
		pstat_on.o	pstat_off.o		pstat_print_all.o pstat_free.o

HEADERS =	pstat.fh	pstat_consts.fh

LIBRARY =	libpstat.a

include ../config/makefile.h
include ../config/makelib.h

$(OBJ): pstatP.fh pstat_consts.fh pstat.fh
