# $Id: GNUmakefile,v 1.3 1997-02-24 20:26:29 d3e129 Exp $

OBJ = seetex.o

LIBRARY = libseetex.a
LIB_TARGETS = seetex $(OBJ)

include ../../src/config/makefile.h

target:
	$(MAKE) seetex
	seetex sample.th sample.F sample.tex

include ../../src/config/makelib.h

install:
	$(MAKE) seetex
	mv seetex $(BINDIR)/seetex
