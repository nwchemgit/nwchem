
           OBJ = md5wrap.o md5.o
       LIBRARY = libutil.a
   LIB_TARGETS = ftest mddriver
       HEADERS = checksum.h

include ../../config/makefile.h
include ../../config/makelib.h

ftest:	ftest.o $(OBJ)
	$(F77) -o $@ $(FFLAGS) $^ -lm

mddriver:	mddriver.o $(OBJ)
	$(CC) -o $@ $^ -lm

