
# $Id: GNUmakefile,v 1.22 1995-10-31 09:21:30 d3e129 Exp $

    HEADERS = api/int_tim.fh api/apiP.fh
    LIBRARY = libnwints.a
    SUBDIRS = api int aux ints_sp deriv texas

include ../config/makefile.h
include ../config/makelib.h


nolib:
	rm -i $(NWCHEM_TOP)/lib/$(NWCHEM_TARGET)/$(LIBRARY)
tasks tasks.ps:	tasks.tex
	latex tasks
	latex tasks
	dvips tasks -o
	rm -f tasks.aux tasks.dvi tasks.log
view tasks.v:	tasks.ps
	ghostview tasks.ps
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
