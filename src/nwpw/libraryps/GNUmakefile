#
# $Id: GNUmakefile,v 1.4 2002-05-03 21:20:25 edo Exp $
#
 OBJ_OPTIMIZE = nwpw_library.o nwpw_libfile.o


	  OBJ = 

      LIBRARY = libnwpwlib.a

    USES_BLAS = 

ifdef NWCHEM_NWPW_LIBRARY_PATH
  LIB_DEFINES =  -DNWPW_LIBRARY="'$(NWCHEM_NWPW_LIBRARY_PATH)'"
else
  LIB_DEFINES =  -DNWPW_LIBRARY="'$(SRCDIR)/nwpw/libraryps/'"
endif

  LIB_INCLUDES = 

include ../../config/makefile.h
include ../../config/makelib.h









