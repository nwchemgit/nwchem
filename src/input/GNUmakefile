#$Id: GNUmakefile,v 1.50 2003-04-17 04:16:16 edo Exp $

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
# 170 MB              = 22282240
# 180 MB              = 23592960
# 200 MB on NERSC T3E = 26214400
  LIB_DEFINES += -DDFLT_TOT_MEM=22282240
endif
ifeq ($(HOSTNAME),ng126)
# 1800 MB = 700 + 700 + 400
 LIB_DEFINES += -DDFLT_TOT_MEM=235929600 -DELAN3HACK
endif
ifeq ($(HOSTNAME),ng127)
# 1800 MB = 700 + 700 + 400
 LIB_DEFINES += -DDFLT_TOT_MEM=235929600 -DELAN3HACK
endif

# END_WIN32_IGNORE

include ../config/makelib.h
