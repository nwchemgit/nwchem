# $Id: GNUmakefile,v 1.4 1994-08-31 19:25:41 d3g681 Exp $

include ../config/makefile.h

#        LIBRARY = librtdb.a
         LIBRARY = libutil.a

             OBJ = rtdb_f2c.o rtdb.o rtdb_seq.o \
                   context.o context_f2c.o

         HEADERS = context.h rtdb.h rtdb.fh context.fh

     LIB_TARGETS = test test.o rtdbtest rtdbtest.o interact context davetest \
                   context.o davetest.o interact.o rtdb_par_f2c.o \
                   cntx.o cntx

    LIB_INCLUDES = -I../db

       TEST_LIBS = $(LIBRARY) $(LIBS) 

include ../config/makelib.h

davetest:	davetest.o $(LIBRARY)
	$(FC) $(FFLAGS) -o $@ davetest.o $(TEST_LIBS) 

cntx:	cntx.o $(LIBRARY) 
	$(FC) $(FFLAGS) -o $@ cntx.o $(LIBS)

interact:	interact.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS) -lm

test:	test.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS) 

rtdbtest:	rtdbtest.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

rtdbpartest:	rtdb_par_test.o $(LIBRARY)
	$(CC) $(CFLAGS) -o $@ rtdb_par_test.o $(TEST_LIBS)

rtdb.o rtdbf2c.o rtdbtest.o:	misc.h
rtdb.o rtdbf2c.o rtdbtest.o:	rtdb.h
test.o:	rtdb.fh

context:	context.o $(LIBRARY)
	$(CC) $(CFLAGS) -o $@ context.o $(TEST_LIBS)

