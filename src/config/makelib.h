# $Id: makelib.h,v 1.1.1.1 1994-03-29 06:44:24 d3g681 Exp $

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
# 8) include ../config/makelib.h
# 9) define any additional targets (e.g., test programs)
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
#
# include ../config/makelib.h
#
# test: test.o $(LIBRARY)
#       $(CC) -o $@ $^
#
# a.o b.o c.o test.o: simple.h
#

$(LIBRARY):	$(OBJ)
	/bin/rm -f $@
	$(AR) $(ARFLAGS) $@ $(OBJ)
	$(RANLIB) $@
	cp -p $(LIBRARY) $(LIBDIR)

ifdef HEADERS
include_stamp:	$(HEADERS)
	cp -p $(HEADERS) $(INCDIR)
	touch include_stamp
else
include_stamp:
	touch include_stamp
endif

clean:
	/bin/rm -f $(LIBRARY) $(OBJ) core include_stamp $(LIB_TARGETS)


realclean:	clean
	/bin/rm -f *~ \#*\#
