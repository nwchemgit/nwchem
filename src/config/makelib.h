# $Id: makelib.h,v 1.5 1994-08-22 02:18:00 d3g681 Exp $

#
# A makefile for a library should
#
# 1) include ../config/makefile.h ... amoung other things this will
#    define TARGET from which any machine dependent actions are driven
# 2) define LIBRARY as the name of the library to be made
# 3) define OBJ as the list of object files to be made
# 4) define HEADERS as the list of header/include files to be exported
#    into the common include directory
# 5) optionally define LIB_TARGETS as any additional files made in
#    this subdirectory that may need cleaning up
# 6) optionally define LIB_DEFINES as any additional defines for
#    the C preprocessor
# 7) optionally define LIB_INCLUDES as any additional includes 
# 8) optionally define SUBDIRS as any subdirectories to build
# 9) include ../config/makelib.h
# 10) define any additional targets (e.g., test programs)
# 11)if the directory contains references to fortran BLAS
#    define USES_BLAS to be the list of FORTRAN files that
#    need converting (e.g., ddot -> sdot)
#
# Note that the library is now put directly into the LIBDIR directory.
#
# E.g.
#
# include ../config/makefile.h
#
#          OBJ = a.o b.o c.o
#      LIBRARY = libsimple.a
#      HEADERS = simple.h
#  LIB_TARGETS = test.o test.x
#  LIB_DEFINES = -DGOODBYE="\"Have a nice day\""
# LIB_INCLUDES = -I../testdir
#    USES_BLAS = test.f
#
# include ../config/makelib.h
#
# test: test.o $(LIBRARY)
#       $(CC) -o $@ $^
#
# a.o b.o c.o test.o: simple.h
#

# Problem with library target is that adding into the library needs
# to be single threaded but parallel make is important
# LIBOBJ = $(patsubst %,$(LIBDIR)/$(LIBRARY)(%),$(OBJ))
#
#$(LIBDIR)/$(LIBRARY):	$(LIBOBJ)
#	$(RANLIB) $@

$(LIBDIR)/$(LIBRARY):	$(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	$(RANLIB) $@

ifdef SUBDIRS
$(LIBDIR)/$(LIBRARY):	subdirs

subdirs:        
	for dir in $(SUBDIRS); do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir || exit 1 ;  \
        done
endif


ifdef HEADERS
include_stamp:	$(HEADERS)
ifdef SUBDIRS
	for dir in $(SUBDIRS); do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir $@ || exit 1 ;  \
	done
endif
	cp -p $(HEADERS) $(INCDIR)
	touch include_stamp

$(OBJ):	$(HEADERS)

else
include_stamp:
ifdef SUBDIRS
	for dir in $(SUBDIRS); do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir $@ || exit 1 ;  \
	done
endif
	touch include_stamp
endif

ifdef USES_BLAS
sngl_to_dbl:
ifdef SUBDIRS
	for dir in $(SUBDIRS); do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir $@ || exit 1 ;  \
	done
endif
	$(CNFDIR)/sngl_to_dbl $(USES_BLAS)
dbl_to_sngl:
ifdef SUBDIRS
	for dir in $(SUBDIRS); do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir $@ || exit 1 ;  \
	done
endif
	$(CNFDIR)/dbl_to_sngl $(USES_BLAS)
else
sngl_to_dbl dbl_to_sngl:
ifdef SUBDIRS
	for dir in $(SUBDIRS); do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir $@ || exit 1 ;  \
	done
endif
	echo $@ : no conversion necessary
endif

clean:
ifdef SUBDIRS
	for dir in $(SUBDIRS); do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir $@ || exit 1 ;  \
	done
endif
	/bin/rm -f *.o *.a core include_stamp $(LIB_TARGETS)


realclean:	clean
ifdef SUBDIRS
	for dir in $(SUBDIRS); do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir $@ || exit 1 ;  \
	done
endif
	/bin/rm -f *~ \#*\#

