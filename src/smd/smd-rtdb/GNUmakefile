# $Id: GNUmakefile,v 1.1 2008-04-21 18:19:27 marat Exp $

include ../../config/makefile.h


      LIBRARY = libsmd.a
    OBJ_OPTIMIZE = srtdb_f2c.o srtdb.o srtdb_seq.o \
                   scontext.o scontext_f2c.o
      HEADERS = srtdb.fh srtdb.h
# BEGIN_WIN32_IGNORE
ifeq ($(TARGET),$(findstring $(TARGET),BGL LINUX))
    OBJ = srtdb_seq.o
    OBJ_OPTIMIZE = srtdb_f2c.o srtdb.o scontext.o scontext_f2c.o
endif
# END_WIN32_IGNORE

#define this to use HDBM ... otherwise DB will be used
#                            DB is broken-ish on 64 bit machines
        USE_HDBM = 1
        DB_INMEM = 1

ifdef USE_HDBM
     LIB_DEFINES = -DUSE_HDBM
else
# BEGIN_WIN32_IGNORE
    LIB_INCLUDES = -I./db
     LIB_DEFINES = -DUSE_DB
         SUBDIRS = db
# END_WIN32_IGNORE
endif

include ../../config/makelib.h

ifeq ($(TARGET),CATAMOUNT)
   LIB_DEFINES += -I$(IOBUF_INC) -I/opt/xt-libc/default/amd64/include
#   COMM_LIBS += $(IOBUF_POST_LINK_OPTS)
endif

test:	test.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS) 
