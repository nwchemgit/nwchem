# $Id: makelib.h,v 1.10 1994-08-27 20:01:16 d3g681 Exp $

#
# A makefile for a library should
#
# 1) include ../config/makefile.h ... amoung other things this will
#    define TARGET from which any machine dependent actions are driven
#    (if you don't need to use TARGET then it's best to include this
#     file at the same point that makelib.h is included).
# 2) define LIBRARY as the name of the library to be made
# 3) optionally define OBJ as the list of object files to be made 
#    without optimization
# 4) optionally define OBJ_OPTIMIZE as the list of object files to 
#    be made with optimization.  It is good practice to keep this list 
#    short so as to minimize exposure to possible compiler errors.
# 5) optionally define HEADERS as the list of header/include files to be
#    copied into the common include directory
# 6) optionally define LIB_TARGETS as any additional files made in
#    this subdirectory that may need cleaning up
# 7) optionally define LIB_DEFINES as any additional defines for
#    the C preprocessor (for both Fortran and C)
# 8) optionally define LIB_INCLUDES as any additional includes 
# 9) optionally define SUBDIRS as any subdirectories to build
#    (note that makefiles in subdirectories will need to modify
#     the paths to the include files)
# 10) optionally define USES_BLAS to be the list of FORTRAN files that
#     need BLAS names converting between single and double (e.g., ddot -> sdot)
# 11) include ../config/makelib.h.  The first rule in this file
#     builds the library so there should be NO targets before this.
# 12) define any additional targets (e.g., test programs)
#
#
# Notes 
# -----

# A)  The library is now put directly into the LIBDIR directory.
#
# B)  To modify the optimization being used specify on the command
#     line C/FDEBUG or C/FOPTIMIZE to override the flags for the 
#     OBJ and OBJ_OPTIMIZE files respectively
#
#         e.g.  make FDEBUG="-g -O1"
#               make FOPTIMIZE="-O3 -Superfast -bugs" FDEBUG="-O2"
#
#
# Sample makefile
# ---------------
#
#          OBJ := a.o b.o c.o
# OBJ_OPTIMIZE := d.o
#      LIBRARY := libsimple.a
#      HEADERS := simple.h
#  LIB_TARGETS := test.o test.x
#  LIB_DEFINES := -DGOODBYE="\"Have a nice day\""
# LIB_INCLUDES := -I../testdir
#    USES_BLAS := test.f
#
# include ../config/makefile.h
# include ../config/makelib.h
#
# test: test.o $(LIBRARY)
#       $(CC) -o $@ $^
#
# a.o b.o c.o test.o: simple.h
#

####################################################################
#
# We use here the default dependency chain (%.o) <- %.o <- %.c/F/...
# However, the default rule for this inserts each .o into the
# library separately which does not work in parallel.  Hence,
# modify the default rule to do nothing but keep the dependency.
# The rule where we usually just do ranlib is also used to build the 
# archive from the objects.  The dependency on the OBJECTS is thus
# also needed since the (%.o): %.o rule no longer puts the .o file
# into the library.
#
####################################################################

OBJECTS := $(OBJ) $(OBJ_OPTIMIZE)

 LIBOBJ := $(patsubst %,$(LIBDIR)/$(LIBRARY)(%),$(OBJECTS))

$(LIBDIR)/$(LIBRARY):   $(LIBOBJ) $(OBJECTS)
ifneq ($(strip $(OBJECTS)),)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)
	$(RANLIB) $@
endif

crap:
	echo $(LIBOBJ) 
	echo $(OBJECTS) 

(%.o):  %.o
	
# Need a tab on the preceding line to make an empty rule

####################################################################
# 
# The definition of C/FFLAGS in makefile.h changes if the variable
# OPTIMIZE is defined.  Without it C/FDEBUG are used.  With it
# C/FOPTIMIZE are used.  To compile the files from OBJ_OPTIMIZE
# with optimization the rule below simply does a make in this
# directory with OPTIMIZE set.  The ifndef OPTIMIZE is to eliminate
# the infinite loop.  The ifdef OBJ_OPTIMIZE is to eliminate the
# empty rule if there are no files to optimize.
# 
####################################################################


ifdef OBJ_OPTIMIZE
ifndef OPTIMIZE
$(OBJ_OPTIMIZE):	optimized
	
# Need a tab on the preceding line to make an empty rule
.PHONY:	optimized
optimized:
	@$(MAKE) OPTIMIZE="Yes" $(OBJ_OPTIMIZE)
endif
endif

####################################################################
#
# Subdirectories are handled with an explicit shell loop propagating
# any make arguments and overrides
#
####################################################################


MAKESUBDIRS = for dir in $(SUBDIRS); do $(MAKE)	 -C $$dir $@ || exit 1 ; done

ifdef SUBDIRS

$(LIBDIR)/$(LIBRARY):	subdirs

.PHONY:	subdirs
subdirs:        
	for dir in $(SUBDIRS); do \
		$(MAKE)	 -C $$dir || exit 1 ;  \
        done
endif


ifdef HEADERS
include_stamp:	$(HEADERS)
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	cp -p $(HEADERS) $(INCDIR)
	touch include_stamp

$(OBJECTS):	$(HEADERS)

else
include_stamp:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	touch include_stamp
endif

ifdef USES_BLAS
.PHONY:	sngl_to_dbl dbl_to_sngl
sngl_to_dbl:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	$(CNFDIR)/sngl_to_dbl $(USES_BLAS)

dbl_to_sngl:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	$(CNFDIR)/dbl_to_sngl $(USES_BLAS)

else

sngl_to_dbl dbl_to_sngl:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	@echo $@ : no conversion necessary
endif

.PHONY:	clean
clean:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	-/bin/rm -f *.o *.a core include_stamp $(LIB_TARGETS)


.PHONY:	realclean
realclean:	clean
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	-/bin/rm -f *~ \#*\#

