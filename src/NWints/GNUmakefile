
# $Id: GNUmakefile,v 1.23 1995-11-01 17:33:49 d3e129 Exp $

    HEADERS = api/int_tim.fh api/apiP.fh
    LIBRARY = libnwints.a
    SUBDIRS = api int aux ints_sp deriv texas

    USE_TEXAS = NOPE
include ../config/makefile.h
include ../config/makelib.h

txs:
	rm -f $(NWCHEM_TOP)/lib/$(NWCHEM_TARGET)/$(LIBRARY)
	rm -f api/cando_sp.o texas/texas_stubs.o
	$(MAKE) "USE_TEXAS=YEP"
fix:
	rm -f $(NWCHEM_TOP)/lib/$(NWCHEM_TARGET)/$(LIBRARY)
	rm -f api/cando_sp.o texas/*.o
	$(MAKE)

nocan:
	rm -f api/cando_sp.o
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
