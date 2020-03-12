#
# $Id$
#

LIBRARY = libpeigs.a

include ../config/makefile.h
include ./DEFS

#peigs.lib:; $(MAKE) -f makefile.peigs lib

peigs_stamp:	include_stamp
		($(MAKE) lib ) || exit 1
		touch peigs_stamp.$(TARGET)

# no include files needed. 
# The tie to an include file is necessary for timing issues in parallel makes RAK 2/16/98

include_stamp:	./h/blas_lapack.h
	touch include_stamp

#.EXPORT_ALL_VARIABLES:

.PHONY :        clean	cleanF	realclean

cleanF:
	@echo " not done here"
cleanDEP:
	@echo " not done here"
clean:
	(cd ./src/c; $(MAKE) clean)
	(cd ./src/f77; $(MAKE) clean)
	(cd ./comm; $(MAKE) clean)
	(cd ./ctof; $(MAKE) clean)

realclean: clean
	rm -f $(LIBDIR)/$(LIBRARY)
	rm -f peigs_stamp.$(TARGET) include_stamp

sngl_to_dbl dbl_to_sngl:
	@echo No conversion necessary

64_to_32:
	(cd ./src/c; $(MAKE) 64_to_32)
	(cd ./src/f77; $(MAKE) 64_to_32)
	(cd ./comm; $(MAKE) 64_to_32)
	(cd ./ctof; $(MAKE) 64_to_32)

32_to_64:
	(cd ./src/c; $(MAKE) 32_to_64)
	(cd ./src/f77; $(MAKE) 32_to_64)
	(cd ./comm; $(MAKE) 32_to_64)
	(cd ./ctof; $(MAKE) 32_to_64)

depend:
	@echo No dependencies in include


source:
	@echo Nothing to be done for $@

lib:
	(cd ./src/c; $(MAKE))
	(cd ./src/f77; $(MAKE))
	(cd ./comm; $(MAKE))
	(cd ./ctof; $(MAKE))
	-rm -f libpeigs.a
	$(peigs_AR) $(peigs_LIB) comm/*.o ctof/*.o src/c/*.o src/f77/*.o
	$(peigs_RANLIB) $(peigs_LIB)
#	($(MAKE) clean)

examples:
	( cd example; $(MAKE) all )

all:
	$(MAKE) peigs.a
	$(MAKE) examples

doc:
	( cd doc; $(MAKE) )





