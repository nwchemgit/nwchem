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

.EXPORT_ALL_VARIABLES:

.PHONY :        clean

clean:
	(cd ./src/c; make clean)
	(cd ./src/f77; make clean)
	(cd ./comm; make clean)
	(cd ./ctof; make clean)
	rm -f peigs_stamp.$(TARGET) include_stamp

realclean: clean
	rm -f $(LIBDIR)/$(LIBRARY)
	rm -f peigs_stamp.$(TARGET) include_stamp


sngl_to_dbl dbl_to_sngl:
	@echo No conversion necessary

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
	$(peigs_AR) libpeigs.a comm/*.o ctof/*.o src/c/*.o src/f77/*.o
	$(peigs_RANLIB) libpeigs.a
	mv libpeigs.a $(LIBDIR)

examples:
	( cd example; $(MAKE) all )

all:
	$(MAKE) peigs.a
	$(MAKE) examples

doc:
	( cd doc; $(MAKE) )





