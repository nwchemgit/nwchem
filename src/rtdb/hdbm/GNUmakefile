#
# $Id$
#

include ../../config/makefile.h

 OBJ_OPTIMIZE = hdbm.o

  LIB_TARGETS = test.o words.o 

  HEADERS = hdbm.h

  LIBRARY = libnwcutil.a

ifdef USE_FFIO
  OBJ_OPTIMIZE += stdiof2ffio.o

  HEADERS += stdiof2ffio.h
endif


include ../../config/makelib.h

ifeq ($(TARGET),CATAMOUNT)
   LIB_DEFINES += -I$(IOBUF_INC)
#   COMM_LIBS += $(IOBUF_POST_LINK_OPTS)
endif

test:	test.o hdbm.o
	$(CC) $(CFLAGS) -o $@ $^

words:	words.o hdbm.o
	$(CC) $(CFLAGS) -o $@ $^

test.o word.o hdbm.o:	hdbm.h

