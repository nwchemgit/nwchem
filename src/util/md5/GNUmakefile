#
# $Id: GNUmakefile,v 1.5 1997-10-31 20:45:42 d3e129 Exp $
#


           OBJ = md5wrap.o md5.o
       LIBRARY = libutil.a
   LIB_TARGETS = ftest mddriver
       HEADERS = checksum.h

include ../../config/makefile.h
include ../../config/makelib.h

ftest:	ftest.o $(OBJ)
	$(LINK.f) -o $@ $^ $(LIBS) $(LIBS) $(LIBS) -lm

mddriver:	mddriver.o $(OBJ)
	$(CC) -o $@ $^ -lm

