
include ../../config/makefile.h

  OBJ = hdbm.o

  LIB_TARGETS = test.o words.o 

  HEADERS = hdbm.h

  LIBRARY = libutil.a

include ../../config/makelib.h


test:	test.o hdbm.o
	$(CC) $(CFLAGS) -o $@ $^

words:	words.o hdbm.o
	$(CC) $(CFLAGS) -o $@ $^

test.o word.o hdbm.o:	hdbm.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(hdbm.o):	 ./hdbm.h 
$(LIBRARY_PATH)(test.o):	 ./hdbm.h 
$(LIBRARY_PATH)(words.o):	 ./hdbm.h 
