
# $Id: GNUmakefile,v 1.27 1995-12-20 22:44:43 d3e129 Exp $

    HEADERS = api/int_tim.fh api/apiP.fh
    LIBRARY = libnwints.a
    SUBDIRS = api int aux ints_sp deriv texas

    USE_TEXAS = NOPE

include ../config/makefile.h
include ../config/makelib.h

txs:
	@touch api/cando_sp.F texas/texas_stubs.F
	$(MAKE) "USE_TEXAS=YEP"
fix:
	@(cd texas; $(MAKE) "USE_TEXAS=YEP" removeobj)
	@touch api/cando_sp.F texas/*.F texas/*.f
	$(MAKE)
nolib:
	rm -i $(LIBRARY_PATH)

tasks tasks.ps:	tasks.tex
	latex tasks
	latex tasks
	dvips tasks -o
	rm -f tasks.aux tasks.dvi tasks.log

view tasks.v:	tasks.ps
	ghostview tasks.ps

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
