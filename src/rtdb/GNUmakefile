# $Id: GNUmakefile,v 1.1.1.1 1994-03-29 06:44:27 d3g681 Exp $

include ../config/makefile.h

             OBJ = rtdb_f2c.o rtdb.o rtdb_par.o rtdb_par_f2c.o \
                   context.o context_f2c.o
         LIBRARY = librtdb.a

         HEADERS = context.h rtdb.h rtdb.fh context.fh

     LIB_TARGETS = test test.o rtdbtest rtdbtest.o interact context davetest \
                   context.o davetest.o interact.o rtdb_par_f2c.o \
                   cntx.o cntx

    LIB_INCLUDES = -I../db

       TEST_LIBS = $(LIBRARY) $(LIBS) 

include ../config/makelib.h

davetest:	davetest.o $(LIBRARY)
	$(FC) $(FFLAGS) -o $@ davetest.o librtdb.a ../db/libdb.a ../util/libutil.a ../ma/libma.a

rtdb_par_f2c.c:	rtdb_f2c.c
	sed 's/rtdb_/rtdb_par_/g' rtdb_f2c.c > rtdb_par_f2c.c

cntx:	cntx.o $(LIBRARY) 
	$(FC) $(FFLAGS) -o $@ cntx.o $(LIBS)

interact:	interact.o $(LIBRARY)
	$(CC) $(CFLAGS) -o $@ interact.o $(TEST_LIBS)

test:	test.o $(LIBRARY)
	$(FC) $(FFLAGS) -o $@ test.o $(TEST_LIBS)

rtdbtest:	rtdbtest.o $(LIBRARY)
	$(CC) $(CFLAGS) -o $@ rtdbtest.o $(TEST_LIBS)

rtdbpartest:	rtdb_par_test.o $(LIBRARY)
	$(CC) $(CFLAGS) -o $@ rtdb_par_test.o $(TEST_LIBS)

rtdb.o rtdbf2c.o rtdbtest.o:	misc.h
rtdb.o rtdbf2c.o rtdbtest.o:	rtdb.h
test.o:	rtdb.fh

context:	context.o $(LIBRARY)
	$(CC) $(CFLAGS) -o $@ context.o $(TEST_LIBS)

