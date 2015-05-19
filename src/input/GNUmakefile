#$Id$

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

ifeq ($(HOSTNAME),ng126)
# 1800 MB = 700 + 700 + 400
 LIB_DEFINES += -DDFLT_TOT_MEM=235929600 -DELAN3HACK
endif
ifeq ($(HOSTNAME),ng127)
# 1800 MB = 700 + 700 + 400
 LIB_DEFINES += -DDFLT_TOT_MEM=235929600 -DELAN3HACK
endif

ifeq ($(HOSTNAME),cu0login1)
  LIB_DEFINES += -DDFLT_TOT_MEM=419430400
endif
ifeq ($(HOSTNAME),cu0login2)
  LIB_DEFINES += -DDFLT_TOT_MEM=419430400
endif
ifeq ($(HOSTNAME),cu0login3)
  LIB_DEFINES += -DDFLT_TOT_MEM=419430400
endif
ifeq ($(HOSTNAME),cu0login4)
  LIB_DEFINES += -DDFLT_TOT_MEM=419430400
endif
ifeq ($(HOSTNAME),cu0login5)
  LIB_DEFINES += -DDFLT_TOT_MEM=419430400
endif

ifdef USE_FORTRAN2003
  LIB_DEFINES += -DUSE_FORTRAN2003
endif

# END_WIN32_IGNORE

ifdef MRCC_METHODS
#ckbn mrcc-r-1 -3
      LIB_DEFINES += -DMRCC_METHODS
endif

include ../config/makelib.h
