# $Id: makelib.h,v 1.27 1995-11-21 19:33:49 gg502 Exp $

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
#     OBJ and OBJ_OPTIMIZE files respectively.
#
#         e.g.  make FDEBUG="-g -O1"
#               make FOPTIMIZE="-O3 -Superfast -bugs" FDEBUG="-O2"
#
# C)  The makelib.h defines the macro LIBRARY_PATH to be the full
#     path of the library being built.
#
# D)  The objectfiles are now deleted
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
# test: test.o $(LIBRARY_PATH)
#       $(CC) -o $@ $^
#
# a.o b.o c.o test.o: simple.h
#

####################################################################
#
# We use here the default dependency chain (%.o) <- %.o <- %.c/F/...
# However, the default rule for this inserts each .o into the
# library separately which does not work in parallel.  
# So we modify the (%o) <- %.o %.c %.F %.f rules so that they do
# not insert directly into the archive.  Putting them in one at
# a time with explicit locking is possible but slow and error 
# prone.
#
# The definition of C/FFLAGS in makefile.h changes if the variable
# OPTIMIZE is defined.  Without it C/FDEBUG are used.  With it
# C/FOPTIMIZE are used.  To compile the files from OBJ_OPTIMIZE
# with optimization the rule below simply does a make in this
# directory with OPTIMIZE set.  The ifndef OPTIMIZE is to eliminate
# the infinite loop.  The ifdef OBJ_OPTIMIZE is to eliminate the
# empty rule if there are no files to optimize.  We don't want to
# ranlib twice so we also modify the default (%.o) rule so that
# it touches .doranlib if a ranlib is necessary.
# 
####################################################################

LIBRARY_PATH := $(LIBDIR)/$(LIBRARY)

.PRECIOUS:	$(LIBRARY_PATH) 

OBJECTS := $(OBJ) $(OBJ_OPTIMIZE)

ifndef OPTIMIZE
ifdef OBJ_OPTIMIZE
 OPT_TARGET = optimized
else
 OPT_TARGET = 
endif

 LIBOBJ := $(patsubst %,$(LIBRARY_PATH)(%),$(OBJ))
$(LIBRARY_PATH):       dummy include_stamp $(LIBOBJ) $(OPT_TARGET)
	@$(MAKE) update_archive

ifdef OBJ_OPTIMIZE
.phony:	optimized
optimized:	
	@$(MAKE) OPTIMIZE=Yes
endif

dummy:
	
else

  LIBOBJ := $(patsubst %,$(LIBRARY_PATH)(%),$(OBJ_OPTIMIZE))
$(LIBRARY_PATH):	$(LIBOBJ)
	
# Previous line must contain tab for empty command
endif

# This puts any floating object files into the library and then
# deletes them. The .notthere is just in case there are no objects
# to be made in a top level directory.

.phony:	update_archive
update_archive:	
	@( list=`for file in $(OBJECTS) .notthere; do if [ -f $$file ] ; then echo $$file; fi ; done`; \
	  if [ "$$list" ] ; then \
		echo $(AR) $(ARFLAGS) $(LIBRARY_PATH) $$list ; \
		     $(AR) $(ARFLAGS) $(LIBRARY_PATH) $$list 2>&1 | \
				grep -v truncated ; \
		/bin/rm -f $$list ; \
		echo ranlib $(LIBRARY_PATH) ; $(RANLIB) $(LIBRARY_PATH) ; \
	  fi; )

# Explicit rules to avoid infinite recursion, to get dependencies, and
# for efficiency

(%.o):	%.F
ifdef EXPLICITF
	@echo Converting $< '->' $*.f
	@$(FCONVERT)
	$(FC) -c $(FFLAGS) -o $% $*.f
	/bin/rm -f $*.f
else
	$(FC) -c $(FFLAGS) $(CPPFLAGS) $<
endif

(%.o):	%.f
	$(FC) -c $(FFLAGS) -o $% $<

(%.o):	%.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $% $<

(%.o):  %.o
	
# Preceding line has a tab to make an empty rule

####################################################################
#
# Subdirectories are handled with an explicit shell loop propagating
# any make arguments and overrides
#
####################################################################


MAKESUBDIRS = @for dir in $(SUBDIRS); do \
			echo Making $@ in $$dir; \
			$(MAKE)	-C $$dir $@ || exit 1 ; done

ifdef SUBDIRS

$(LIBRARY_PATH):	subdirs

.PHONY:	subdirs
subdirs:        
	@for dir in $(SUBDIRS); do \
		echo Making all in $$dir; \
		$(MAKE)	 -C $$dir || exit 1 ;  \
        done
endif

.PHONY:	depend
depend:	
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	@echo Making depend in `pwd`
	$(CNFDIR)/makedepend $(LIB_INCLUDES)


ifdef HEADERS
include_stamp:	$(HEADERS)
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	cp -p $(HEADERS) $(INCDIR)
	touch include_stamp

# The below dependency has now been superceded by the complete
# dependency analysis of headers of makedepend/dependencies
#$(OBJECTS):	$(HEADERS)

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
	-$(RM) *.o *.a core include_stamp $(LIB_TARGETS)


.PHONY:	cleanF
cleanF:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	@for file in *F; do \
		body=`basename $$file .F` ; \
		if [ -f $$body.f ] ; then \
		  echo $$file and $$body.f both exist ... deleting $$body.f; \
		  /bin/rm -f $$body.f ; \
		fi ; \
        done
	


.PHONY:	realclean
realclean:	clean
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	-$(RM) *~ \#*\# makefile.bak $(LIBRARY_PATH)

