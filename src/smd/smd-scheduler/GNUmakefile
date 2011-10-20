#
# $Id$
#

          OBJ = sched_groups.o sched_perf.o 
  LIB_DEFINES = -DMPI
  ifneq ($(MPI_INCLUDE),)
    LIB_DEFINES += $(patsubst -I-I%,-I%,-I$(MPI_INCLUDE))
  endif

LIB_INCLUDES = -I../include


      LIBRARY = libsmd.a

include ../../config/makefile.h
include ../../config/makelib.h









