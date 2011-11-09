# $Id$

include ../config/makefile.h

         LIBRARY = libnwcutil.a

    OBJ_OPTIMIZE = rtdb_f2c.o rtdb.o rtdb_seq.o context.o context_f2c.o

# BEGIN_WIN32_IGNORE

# replace TARGET with NWCHEM_TARGET so ARMCI build system headers are not required
ifeq ($(NWCHEM_TARGET),$(findstring $(NWCHEM_TARGET),BGL BGP LINUX))
    OBJ = rtdb_seq.o
    OBJ_OPTIMIZE = rtdb_f2c.o rtdb.o context.o context_f2c.o
endif

# END_WIN32_IGNORE

         HEADERS = context.h rtdb.h rtdb.fh context.fh rtdb.cray.h

     LIB_TARGETS = test test.o rtdbtest rtdbtest.o interact context davetest \
                   context.o davetest.o interact.o rtdb_par_f2c.o \
                   cntx.o cntx testgr.o testgr

       TEST_LIBS = $(LIBRARY) $(LIBS) 

#define this to use HDBM ... otherwise DB will be used
#                            DB is broken-ish on 64 bit machines
        USE_HDBM = 1

ifdef USE_HDBM
     LIB_DEFINES = -DUSE_HDBM
         SUBDIRS = hdbm
else
# BEGIN_WIN32_IGNORE
    LIB_INCLUDES = -I./db
     LIB_DEFINES = -DUSE_DB
         SUBDIRS = db
# END_WIN32_IGNORE
endif

include ../config/makelib.h

ifeq ($(TARGET),CATAMOUNT)
   LIB_DEFINES += -I$(IOBUF_INC) -I/opt/xt-libc/default/amd64/include
#   COMM_LIBS += $(IOBUF_POST_LINK_OPTS)
endif


davetest:	davetest.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ davetest.o $(LIBS)

cntx:	cntx.o  $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^  $(LIBS)

interact:	interact.o $(LIBRARY_PATH)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ interact.o $(LIBRARY_PATH) -lglobal -ltcgmsg -L/msrc/apps/f771.4/SC1.0 -lF77 -lV77 -lm

test:	test.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS) 
testgr:	testgr.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS) 

rtdbtest:	rtdbtest.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ -lm $(LIBS)

rtdbpartest:	rtdb_par_test.o $(LIBRARY)
	$(CC) $(CFLAGS) -o $@ rtdb_par_test.o $(TEST_LIBS)

context:	context.o $(LIBRARY)
	$(CC) $(CFLAGS) -o $@ context.o $(TEST_LIBS)

