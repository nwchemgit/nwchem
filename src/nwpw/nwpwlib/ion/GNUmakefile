#
# $Id$
#
 SUBDIRS = shake

 OBJ_OPTIMIZE = \
            center.o \
			ion_FixIon.o \
            incell1.o \
            incell2.o \
            incell3.o \
			xyz.o \
            ion.o \
			katm_list.o \
            fcoord.o \
			seperate_molpsp.o \
			seperate_pointcharge.o


	  OBJ = 

      LIBRARY = libnwpwlib.a

    USES_BLAS = ion.F fcoord.f  incell1.f incell3.f seperate_molpsp.F seperate_pointcharge.F incell2.f ion_FixIon.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

