# $Id: GNUmakefile,v 1.4 1997-05-28 07:07:48 d3e129 Exp $

OBJ = seetex.o

LIBRARY = libseetex.a
LIB_TARGETS = seetex $(OBJ)

include ../../src/config/makefile.h

target:
	$(MAKE) seetex
	rm -f sample.tex
	seetex sample.th sample.F sample.tex

include ../../src/config/makelib.h

install:
	$(MAKE) seetex
	mv seetex $(BINDIR)/seetex
