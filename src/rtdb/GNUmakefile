# $Id: GNUmakefile,v 1.19 1999-07-29 01:03:13 d3e129 Exp $

include ../config/makefile.h

         LIBRARY = libutil.a


    OBJ_OPTIMIZE = rtdb_f2c.o rtdb.o rtdb_seq.o \
                   context.o context_f2c.o
ifeq ($(TARGET),LINUX)
    OBJ = rtdb_seq.o
    OBJ_OPTIMIZE = rtdb_f2c.o rtdb.o context.o context_f2c.o
endif
         HEADERS = context.h rtdb.h rtdb.fh context.fh rtdb.cray.h

     LIB_TARGETS = test test.o rtdbtest rtdbtest.o interact context davetest \
                   context.o davetest.o interact.o rtdb_par_f2c.o \
                   cntx.o cntx

       TEST_LIBS = $(LIBRARY) $(LIBS) 

#define this to use HDBM ... otherwise DB will be used
#                            DB is broken-ish on 64 bit machines
        USE_HDBM = 1

ifdef USE_HDBM
     LIB_DEFINES = -DUSE_HDBM
         SUBDIRS = hdbm
else
    LIB_INCLUDES = -I./db
     LIB_DEFINES = -DUSE_DB
         SUBDIRS = db
endif

include ../config/makelib.h

davetest:	davetest.o $(LIBRARY)
	$(FC) $(FFLAGS) -o $@ davetest.o $(TEST_LIBS) 

cntx:	cntx.o $(LIBRARY) 
	$(FC) $(FFLAGS) -o $@ cntx.o $(LIBS)

interact:	interact.o $(LIBRARY_PATH)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ interact.o $(LIBRARY_PATH) -lglobal -ltcgmsg -L/msrc/apps/f771.4/SC1.0 -lF77 -lV77 -lm

test:	test.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS) 

rtdbtest:	rtdbtest.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ -lmalloc $(LIBS)

rtdbpartest:	rtdb_par_test.o $(LIBRARY)
	$(CC) $(CFLAGS) -o $@ rtdb_par_test.o $(TEST_LIBS)

context:	context.o $(LIBRARY)
	$(CC) $(CFLAGS) -o $@ context.o $(TEST_LIBS)

