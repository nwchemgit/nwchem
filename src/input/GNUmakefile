#$Id: GNUmakefile,v 1.13 1995-02-02 23:15:20 d3g681 Exp $
 LIBRARY = libutil.a
#LIBRARY = libinput.a

     OBJ = input_parse.o input_mem_size.o memory_input.o input_set.o \
           input_start_opt.o input_title.o task_input.o scf_input.o \
           get_input_filename.o vectors_input.o stepper_input.o input_echo.o \
	   dft_input.o md_input.o gradients_input.o

include ../config/makefile.h
include ../config/makelib.h

ifeq ($(TARGET),SGITFP)
get_input_filename.o:	get_input_filename.F
	f77 -g -d8 $(INCLUDES) $(DEFINES)  -c get_input_filename.F -o get_input_filename.o	
endif
