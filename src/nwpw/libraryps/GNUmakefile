#
# $Id: GNUmakefile,v 1.5 2003-12-02 19:16:53 bylaska Exp $
#
 OBJ_OPTIMIZE = nwpw_library.o nwpw_libfile.o nwpw_set_library.o


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









