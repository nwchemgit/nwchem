#
# $Id: GNUmakefile,v 1.8 1997-10-31 23:51:45 d3e129 Exp $
#

include ../../config/makefile.h

 OBJ_OPTIMIZE = hdbm.o

  LIB_TARGETS = test.o words.o 

  HEADERS = hdbm.h

  LIBRARY = libutil.a

include ../../config/makelib.h


test:	test.o hdbm.o
	$(CC) $(CFLAGS) -o $@ $^

words:	words.o hdbm.o
	$(CC) $(CFLAGS) -o $@ $^

test.o word.o hdbm.o:	hdbm.h

