# $Id$

include ../../tools/armci/config/makecoms.h
include ../../config/makefile.h

         LIBRARY = libsmd.a


    OBJ_OPTIMIZE = srtdb_f2c.o srtdb.o srtdb_seq.o \
                   scontext.o scontext_f2c.o
         HEADERS = srtdb.h srtdb.fh 

     LIB_TARGETS = test test.o 

       TEST_LIBS = $(LIBRARY) $(LIBS) 

     LIB_DEFINES = -DUSE_HDBM

include ../../config/makelib.h

test:	test.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS) 


