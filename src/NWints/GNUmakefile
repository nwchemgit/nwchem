
# $Id: GNUmakefile,v 1.21 1995-10-17 08:52:27 d3g681 Exp $

    HEADERS = api/int_tim.fh api/apiP.fh
    LIBRARY = libnwints.a
    SUBDIRS = api int aux ints_sp deriv texas

include ../config/makefile.h
include ../config/makelib.h


tasks tasks.ps:	tasks.tex
	latex tasks
	latex tasks
	dvips tasks -o
	rm -f tasks.aux tasks.dvi tasks.log
view tasks.v:	tasks.ps
	ghostview tasks.ps
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
