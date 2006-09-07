#
# $Id: GNUmakefile,v 1.7 2006-09-07 16:55:09 d3p307 Exp $
#
 OBJ_OPTIMIZE = nwpw_library.o nwpw_libfile.o nwpw_set_library.o


	  OBJ = 

      LIBRARY = libnwpwlib.a

    USES_BLAS = 

   LONGNWTOP = $(shell if [ `echo $(SRCDIR) |wc -m` -gt 64 ]; then echo "Y"; fi  )
ifeq ($(LONGNWTOP),Y)
errorlongpwd:
	@echo " "
	@echo "The directory name chosen for NWCHEM_TOP is longer than "
	@echo "the maximum allowed value of 65 characters"
	@echo "      current NWCHEM_TOP="$(SRCDIR)" equal to" $(shell echo $(BSRCDIR) |wc -m) "characters"
	@echo "please chose a directory with a shorter name"
	@echo "E.g., setenv NWCHEM_TOP /home/nwchem_user/nwchem"
	@exit 1
endif
  LIB_DEFINES =  -DNWPW_LIBRARY="'$(SRCDIR)'"

  LIB_INCLUDES = 

include ../../config/makefile.h
include ../../config/makelib.h









