
# $Id: GNUmakefile,v 1.34 1996-06-27 22:33:24 d3g681 Exp $

    HEADERS = api/int_tim.fh api/apiP.fh seint/se_defs.fh
    LIBRARY = libnwints.a
    SUBDIRS = api int aux ints_sp deriv texas

    USE_TEXAS = YEP


include ../config/makefile.h
include ../config/makelib.h


semi:
	cp seint/se_defs.fh $(INCDIR)
	(cd seint; $(MAKE) )
txs:
	@touch api/cando_sp.F texas/texas_stubs.F api/int_2e2c.F api/int_2e3c.F
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
