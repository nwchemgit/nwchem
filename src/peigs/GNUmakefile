LIBRARY = libpeigs.a

include ../config/makefile.h
include ./DEFS

#peigs.lib:; $(MAKE) -f makefile.peigs lib

peigs_stamp:	$(LIBDIR)/$(LIBRARY)
	($(MAKE) -f makefile.peigs lib ) || exit 1
	touch peigs_stamp


# no include files needed
include_stamp:
	touch include_stamp

clean:
	$(MAKE) -f makefile.peigs clean
	rm -f peigs_stamp include_stamp
realclean:	clean
	$(MAKE) -f makefile.peigs realclean
	rm -f $(LIBDIR)/$(LIBRARY)
	rm -f peigs_stamp include_stamp


sngl_to_dbl dbl_to_sngl:
	@echo No conversion necessary

depend:
	@echo No dependencies in include


source:
	@echo Nothing to be done for $@
