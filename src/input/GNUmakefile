#$Id: GNUmakefile,v 1.43 1999-06-08 00:47:02 d3h325 Exp $

# Put into task library so that it occurs at the front of the link path

 LIBRARY = libnwctask.a

     OBJ = input_parse.o input_mem.o memory_input.o input_set.o \
           input_title.o \
           get_inp_file.o \
	   memory_def.o \
	   speech_input.o

#     HEADERS = ../global/src/types.f2c.h

include ../config/makefile.h

HOSTNAME := $(shell hostname)

HOSTSUFFIX := $(shell hostname | sed -e 's/et....//' | sed -e 's/^\.//')

ifeq ($(HOSTNAME),mcurie)
# 180 MB              = 23592960
# 200 MB on NERSC T3E = 26214400
  LIB_DEFINES += -DDFLT_TOT_MEM=23592960
endif

ifeq ($(HOSTSUFFIX),nwmpp1.emsl.pnl.gov)
# 450 MB on EMSL IBMSP
  LIB_DEFINES += -DDFLT_TOT_MEM=58982400
endif


include ../config/makelib.h

ifeq ($(TARGET),SGITFP)
$(LIBRARY_PATH)(get_inp_file.o):	get_inp_file.F
	f77 -g -d8 -mips4 -64 $(INCLUDES) $(DEFINES)  -c get_inp_file.F
endif


