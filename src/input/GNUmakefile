#$Id: GNUmakefile,v 1.47 2000-08-02 01:25:17 bjohnson Exp $

# Put into task library so that it occurs at the front of the link path

 LIBRARY = libnwctask.a

     OBJ = input_parse.o input_mem.o memory_input.o input_set.o \
           input_title.o \
           get_inp_file.o \
	   memory_def.o \
	   speech_input.o  

#     HEADERS = ../global/src/types.f2c.h

include ../config/makefile.h

# BEGIN_WIN32_IGNORE

HOSTNAME := $(shell hostname)

HOSTSUFFIX := $(shell hostname | sed -e 's/et....//' | sed -e 's/^\.//')

ifeq ($(HOSTNAME),mcurie)
# 180 MB              = 23592960
# 200 MB on NERSC T3E = 26214400
  LIB_DEFINES += -DDFLT_TOT_MEM=23592960
endif

ifeq ($(HOSTSUFFIX),nwmpp1.emsl.pnl.gov)
# 384 MB on EMSL IBMSP
  LIB_DEFINES += -DDFLT_TOT_MEM=50331648
endif

# END_WIN32_IGNORE

include ../config/makelib.h
