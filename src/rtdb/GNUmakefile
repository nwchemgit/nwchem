# $Id: GNUmakefile,v 1.11 1995-11-19 19:43:16 d3g681 Exp $

include ../config/makefile.h

         LIBRARY = libutil.a

             OBJ = rtdb_f2c.o rtdb.o rtdb_seq.o \
                   context.o context_f2c.o

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

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(cntx.o):	 rtdb.fh $(INCDIR)/mafdecls.fh context.fh 
$(LIBRARY_PATH)(context.o):	 rtdb.h misc.h $(INCDIR)/macdecls.h 
$(LIBRARY_PATH)(context_f2c.o):	 context.h 
$(LIBRARY_PATH)(davetest.o):	 rtdb.fh $(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(interact.o):	 rtdb.h $(INCDIR)/macdecls.h 
$(LIBRARY_PATH)(rtdb.o):	 $(INCDIR)/sndrcv.h rtdb.h $(INCDIR)/macdecls.h 
$(LIBRARY_PATH)(rtdb.save.o):	 $(INCDIR)/sndrcv.h rtdb.h misc.h $(INCDIR)/macdecls.h 
$(LIBRARY_PATH)(rtdb_f2c.o):	 rtdb.h $(INCDIR)/macdecls.h 
$(LIBRARY_PATH)(rtdb_seq.o):	 $(INCDIR)/sndrcv.h rtdb.h misc.h $(INCDIR)/macdecls.h $(INCDIR)/hdbm.h 
$(LIBRARY_PATH)(rtdbtest.o):	 $(INCDIR)/sndrcv.h rtdb.h misc.h $(INCDIR)/macdecls.h 
$(LIBRARY_PATH)(test.o):	 rtdb.fh $(INCDIR)/mafdecls.fh 
