
# $Id: GNUmakefile,v 1.25 1995-11-21 19:35:34 gg502 Exp $

    HEADERS = api/int_tim.fh api/apiP.fh
    LIBRARY = libnwints.a
    SUBDIRS = api int aux ints_sp deriv texas

    USE_TEXAS = NOPE
include ../config/makefile.h
include ../config/makelib.h

txs:
	rm -f $(LIBRARY_PATH)
	rm -f api/cando_sp.o texas/texas_stubs.o
	$(MAKE) "USE_TEXAS=YEP"
fix:
	rm -f $(LIBRARY_PATH)
	rm -f api/cando_sp.o texas/*.o
	$(MAKE)

nocan:
	rm -f api/cando_sp.o

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
