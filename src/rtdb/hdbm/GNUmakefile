#
# $Id: GNUmakefile,v 1.9 2001-08-28 00:54:38 edo Exp $
#

include ../../config/makefile.h

 OBJ_OPTIMIZE = hdbm.o

  LIB_TARGETS = test.o words.o 

  HEADERS = hdbm.h

  LIBRARY = libnwcutil.a

include ../../config/makelib.h


test:	test.o hdbm.o
	$(CC) $(CFLAGS) -o $@ $^

words:	words.o hdbm.o
	$(CC) $(CFLAGS) -o $@ $^

test.o word.o hdbm.o:	hdbm.h

