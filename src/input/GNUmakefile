#$Id: GNUmakefile,v 1.35 1997-02-06 18:56:24 mdupuis Exp $

# Put into task library so that it occurs at the front of the link path

 LIBRARY = libtask.a

     OBJ = input_parse.o input_mem.o memory_input.o input_set.o \
           input_title.o \
           get_inp_file.o \
	   memory_def.o \
	   speech_input.o

     HEADERS = ../global/src/types.f2c.h

include ../config/makefile.h
include ../config/makelib.h

ifeq ($(TARGET),SGITFP)
$(LIBRARY_PATH)(get_inp_file.o):	get_inp_file.F
	f77 -g -d8 $(INCLUDES) $(DEFINES)  -c get_inp_file.F
endif

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(get_inp_file.o):	$(INCDIR)/inp.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(input_parse.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(speech_input.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(input_set.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(memory_input.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/context.fh 
$(LIBRARY_PATH)(memory_def.o):	$(INCDIR)/inp.fh 
$(LIBRARY_PATH)(input_title.o):	$(INCDIR)/rtdb.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(input_mem.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh 
