#
# $Id: GNUmakefile,v 1.10 2004-02-11 02:27:36 edo Exp $
#

include ../../config/makefile.h

 OBJ_OPTIMIZE = hdbm.o

  LIB_TARGETS = test.o words.o 

  HEADERS = hdbm.h

  LIBRARY = libnwcutil.a

#ifdef USE_FFIO
  OBJ_OPTIMIZE += stdiof2ffio.o

  HEADERS += stdiof2ffio.h
#endif


include ../../config/makelib.h


test:	test.o hdbm.o
	$(CC) $(CFLAGS) -o $@ $^

words:	words.o hdbm.o
	$(CC) $(CFLAGS) -o $@ $^

test.o word.o hdbm.o:	hdbm.h

