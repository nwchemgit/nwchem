# $Id: GNUmakefile,v 1.5 1994-09-25 07:54:49 d3e129 Exp $


#LIBRARY = libpstat.a
LIBRARY =      libutil.a

OBJ =		pstat_init.o	pstat_terminate.o	pstat_allocate.o \
		pstat_on.o	pstat_off.o		pstat_print_all.o \
                pstat_free.o	pstat_print_handles.o

HEADERS =	pstat.fh	pstat_consts.fh


include ../config/makefile.h
include ../config/makelib.h

$(OBJ): pstatP.fh pstat_consts.fh pstat.fh
