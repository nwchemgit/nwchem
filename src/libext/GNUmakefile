#$Id$
# makefile,v 1.8 1994/12/05 20:37:08 og845 Exp

# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../config/makefile.h


SUBDIRS=
ifdef BUILD_OPENBLAS
    SUBDIRS += openblas
endif
ifdef BUILD_MPICH
    SUBDIRS += mpich
endif
ifdef BUILD_SCALAPACK
    SUBDIRS += scalapack
endif
ifdef BUILD_ELPA
    SUBDIRS += elpa
endif
ifdef USE_LIBXC
    SUBDIRS += libxc
endif
ifdef USE_TBLITE
    SUBDIRS += tblite
endif
ifdef BUILD_PLUMED
    SUBDIRS += plumed
endif

include ../config/makelib.h

64_to_32 32_to_64:
	@echo No blas conversion done here





