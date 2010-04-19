#
# $Id: GNUmakefile,v 1.14 2008/05/20 16:59:46 marat Exp $
#

          OBJ = sched_groups.o sched_perf.o 
  LIB_DEFINES = -DMPI
  ifneq ($(MPI_INCLUDE),)
    LIB_DEFINES += -I$(MPI_INCLUDE)
    endif

LIB_INCLUDES = -I../include


      LIBRARY = libsmd.a

include ../../config/makefile.h
include ../../config/makelib.h









