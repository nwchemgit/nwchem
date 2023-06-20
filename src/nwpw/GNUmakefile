#
# $Id$
#

#      SUBDIRS = pspw band nwpwlib libraryps
       SUBDIRS = nwpwlib libraryps pspw paw band ofpw

          OBJ = nwpw_input.o \
                nwpw_default.o \
                nwpw_vectors_input.o


      LIBRARY = libnwpw.a
$(info USE_OPENMP $(origin USE_OPENMP))

ifdef USE_OPENMP
# QA issues
ifeq ($(FC),ifx)
override undefine USE_OPENMP
endif
endif

$(info USE_OPENMP $(origin USE_OPENMP))

  LIB_DEFINES =

      HEADERS = 




include ../config/makefile.h
include ../config/makelib.h









