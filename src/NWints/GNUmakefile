
# $Id: GNUmakefile,v 1.44 1999-11-02 23:07:14 bert Exp $

    HEADERS = api/int_tim.fh api/apiP.fh seint/se_defs.fh api/int_nbf.fh aux/spcartP.fh hondo/hnd_whermt.fh hondo/hnd_rys.fh hondo/hnd_tol.fh
    LIBRARY = libnwints.a
    SUBDIRS = api int aux ints_sp deriv texas ecp hondo rel dk
    OBJ = int_input.o
  USE_TEXAS = YEP
export USE_TEXAS


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

