#$Id: GNUmakefile,v 1.38 1998-08-04 19:13:51 d3g681 Exp $

# Put into task library so that it occurs at the front of the link path

 LIBRARY = libtask.a

     OBJ = input_parse.o input_mem.o memory_input.o input_set.o \
           input_title.o \
           get_inp_file.o \
	   memory_def.o \
	   speech_input.o

     HEADERS = ../global/src/types.f2c.h

include ../config/makefile.h

HOSTNAME := $(shell hostname)

ifeq ($(HOSTNAME),mcurie)
  LIB_DEFINES += -DDFLT_TOT_MEM=26214400
endif
ifeq ($(TARGET),SGITFP)
$(LIBRARY_PATH)(get_inp_file.o):	get_inp_file.F
	f77 -g -d8 -mips4 -64 $(INCLUDES) $(DEFINES)  -c get_inp_file.F
endif

include ../config/makelib.h


