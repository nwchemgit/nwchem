
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

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(mddriver.o):	./mdglobal.h 
$(LIBRARY_PATH)(md5wrap.o):	./md5.h 
$(LIBRARY_PATH)(md5.o):	./md5.h 
