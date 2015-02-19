#$Id$

include ../../config/makefile.h

# LIBRARY = libdb.a

  LIBRARY = libnwcutil.a


  LIB_TARGETS = rjhtest rjhtest.o

   OBJ_MOST = db.o    hash_bigkey.o   hash_func.o  hash_page.o \
            hash.o  hash_buf.o      hash_log2.o  mkstemp.o 

# hsearch.o

  EXTRA = memmove.o realloc.o

ifeq ($(TARGET),IBM)
  EXTRA = realloc.o
endif
ifeq ($(TARGET),SP)
  EXTRA = realloc.o
endif

# default definition
ifndef LIB_DEFINES
  LIB_DEFINES = -DBYTE_ORDER=BIG_ENDIAN -DMEMMOVE
endif

      OBJ = $(OBJ_MOST) $(EXTRA)


include ../../config/makelib.h

rjhtest:	rjhtest.o $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $^


db.c:		db.h
hash.c:		db.h
hash.c:		hash.h
hash.c:		page.h
hash.c:		extern.h
hash_bigkey.c:	db.h
hash_bigkey.c:	hash.h
hash_bigkey.c:	page.h
hash_bigkey.c:	extern.h
hash_buf.c:	db.h
hash_buf.c:	hash.h
hash_buf.c:	page.h
hash_buf.c:	extern.h
hash_func.c:	db.h
hash_func.c:	hash.h
hash_func.c:	page.h
hash_func.c:	extern.h
hash_page.c:	db.h
hash_page.c:	hash.h
hash_page.c:	page.h
hash_page.c:	extern.h
hsearch.c:	db.h
hsearch.c:	search.h
rjhtest.c:	db.h
rjhtest.c:	rjhmisc.h
db.h:		cdefs.h compat.h
