#$Id: GNUmakefile,v 1.37 1997-03-21 23:52:20 d3e129 Exp $

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
	f77 -g -d8 -mips4 -64 $(INCLUDES) $(DEFINES)  -c get_inp_file.F
endif

