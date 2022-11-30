# $Id$

    HEADERS = api/int_tim.fh api/apiP.fh api/int_nbf.fh auxil/spcartP.fh hondo/hnd_whermt.fh hondo/hnd_rys.fh hondo/hnd_tol.fh ../property/prop.fh
    LIBRARY = libnwints.a
    SUBDIRS = simint api int auxil ints_sp deriv texas ecp hondo rel dk sifs
    OBJ = intgrl_input.o
  USE_TEXAS = YEP
export USE_TEXAS


include ../config/makefile.h
include ../config/makelib.h


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

