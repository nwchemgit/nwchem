# $Id$

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
ifdef USE_SHARED
LIBRARY_SO := $(shell echo $(LIBRARY) | sed s/\\.a/\\.so/g )
endif

.PRECIOUS:	$(LIBRARY_PATH) 

OBJECTS := $(OBJ) $(OBJ_OPTIMIZE)

######################################################################
# makefile in each directory might define error message for undefined 
# symbols etc. When error message is defined, it should be displayed 
# and then make processing aborted.
# This comes from GA, TCGMSG-MPI & DA
######################################################################

define print_error
	@echo $(ERRMSG)
	exit 1
	@echo
endef

ifdef ERRMSG
error:
	$(print_error)
endif

# Make sure that nothing gets compiled in case of error
#
ifdef ERRMSG
      CC = $(print_error)
      FC = $(print_error)
endif

######################################################################
ifndef OPTIMIZE
ifdef OBJ_OPTIMIZE
 OPT_TARGET = optimized
else
 OPT_TARGET = 
endif

all:	$(LIBRARY_PATH) $(LIB_ALSO_BUILD)


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
# to be made in a top level directory.  The locking is necessary
# because concurrent makes in different directories may attempt
# to update the same archive.  The locking mechanism works on a
# coarse grain like this but is too slow for finer grain locking

    LOCKFILE := $(LIBRARY_PATH:.a=.lock)

.phony: update_archive
update_archive: 
ifdef USE_SHARED
	@( list=`for file in $(OBJECTS) .notthere; do if [ -f $$file ] ; then echo $$file; fi ; done`; \
	  if [ "$$list" ] ; then \
	        $(CNFDIR)/lockfile -steal $(LOCKFILE) || exit 1 ; \
	        echo $(CC) -shared -o $(LIBDIR)/$(LIBRARY_SO) $$list ;\
                $(CC) -shared -o $(LIBDIR)/$(LIBRARY_SO) $$list ;\
	        /bin/rm -f $(LOCKFILE) ; \
	  fi; )
endif
ifdef SUMO
	@( list=`for file in $(OBJECTS) .notthere; do if [ -f $$file ] ; then echo $$file; fi ; done`; \
	  if [ "$$list" ] ; then \
	        $(CNFDIR)/lockfile -steal $(LOCKFILE) || exit 1 ; \
	        if test ! -e $(LIBDIR)/objs; then mkdir $(LIBDIR)/objs; fi;\
	        cp $$list $(LIBDIR)/objs/;\
	        /bin/rm -f $(LOCKFILE) ; \
	  fi; )
endif
	@( list=`for file in $(OBJECTS) .notthere; do if [ -f $$file ] ; then echo $$file; fi ; done`; \
	  if [ "$$list" ] ; then \
	        $(CNFDIR)/lockfile -steal $(LOCKFILE) || exit 1 ; \
                echo  $(AR) $(ARFLAGS) $(LIBRARY_PATH) $$list ; \
  		echo $$list|grep -v truncated | xargs -n 300 $(AR) $(ARFLAGS) $(LIBRARY_PATH) 2>&1  ; \
	        /bin/rm -f $$list ; \
	        echo $(RANLIB) $(LIBRARY_PATH) ; $(RANLIB) $(LIBRARY_PATH) ; \
	        /bin/rm -f $(LOCKFILE) ; \
	   fi; )


#
# The explict rules for (%.o), .o, .f and .F have been moved
# back into makefile.h since there is machine dependence generated
# at least by the T3D
#

####################################################################
#
# Subdirectories are handled with an explicit shell loop propagating
# any make arguments and overrides
#
####################################################################


MAKESUBDIRS = +@for dir in $(SUBDIRS); do \
			echo Making $@ in $$dir; \
			$(MAKE)	-C $$dir $@ || exit 1 ; done

ifdef SUBDIRS
ifndef OPTIMIZE

$(LIBRARY_PATH):	subdirs

.PHONY:	subdirs
subdirs:        
	@for dir in $(SUBDIRS); do \
		echo Making all in $$dir; \
		$(MAKE)	 -C $$dir || exit 1 ;  \
        done
endif
endif


ifdef HEADERS
include_stamp:	$(HEADERS)
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
ifdef USE_CPPRESERVE
	cp --preserve=timestamp $(HEADERS) $(INCDIR)
else
	cp -p $(HEADERS) $(INCDIR)
endif

	touch include_stamp

# The below dependency has now been superceded by the complete
# dependency analysis of headers of makedepend/dependencies (now depend.x)
#$(OBJECTS):	$(HEADERS)

else
include_stamp:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	touch include_stamp
endif

ifdef USES_BLAS
.PHONY:	sngl_to_dbl dbl_to_sngl 64_to_32 32_to_64
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

64_to_32:
	$(CNFDIR)/64_to_32 $(USES_BLAS)
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	@/bin/rm -f dependencies

32_to_64:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	$(CNFDIR)/32_to_64 $(USES_BLAS)
	@/bin/rm -f dependencies

else

sngl_to_dbl dbl_to_sngl 64_to_32 32_to_64:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	@/bin/rm -f dependencies
	@echo $@ : no conversion necessary
endif

.PHONY:	clean
clean:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	-$(RM) -f *.o *.a *.mod *__genmod.f90 *core *stamp *trace mputil.mp* *events* *ipo *optrpt $(LIB_TARGETS)
	if [ -f $(LIBRARY_PATH) ] ; then \
  		echo $(OBJ) $(OBJ_OPTIMIZE)| xargs -n 300 $(AR) d $(LIBRARY_PATH) ; \
		if [ `$(AR) t $(LIBRARY_PATH) | wc | awk ' {print $$1;}'` -eq 0 ] ; then \
			$(RM) -f $(LIBRARY_PATH) ; \
		fi ; \
	fi ;
	@-$(RM) dependencies

#
# This is a convenience target that will make the TAGS file for current 
# checked out source tree.  This is only useful if you know something 
# about emacs.  Note: find, grep and etags must be in your path.
#
.PHONY: tags_file
tags_file:
	find . \( -name "*.[cfFh]" -o -name "*.fh" \) -print | grep -v "\./include" | grep -v "\./tools/include" | grep -v "NWints/seint" | etags -
.PHONY: tags_clean
tags_clean:
	find . -name TAGS -print -exec rm -f "{}" ";"
	
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
realclean:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	-$(RM) -f *~ \#*\#  makefile.bak $(LIBRARY_PATH)
	$(MAKE) clean
	-$(RM) dependencies
ifdef SUMO
	-$(RM) -rf $(LIBDIR)/objs
endif


.PHONY:	source
source:
	@/bin/rm -f source
	@touch source
	@for file in nonexistent $(OBJ) $(OBJ_OPTIMIZE); do \
		body=`basename $$file .o` ; \
		if [ -f $$body.f ] ; then \
		  echo $$body.f; \
		  cat $$body.f >> source; \
		fi ; \
		if [ -f $$body.F ] ; then \
		  echo $$body.F; \
		  cat $$body.F >> source; \
		fi ; \
        done
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif

.PHONY:	cleanDEP
cleanDEP:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	@/bin/rm -f dependencies

#
# If make cannot find the dependencies file it will generate it
# using this rule.  We also need to make sure that the program
# that makes the dependencies has been built.
#
$(BINDIR)/depend.x:	
	( cd $(CNFDIR); $(MAKE) $@ ; )

dependencies:	$(wildcard *.c) $(wildcard *.F) $(BINDIR)/depend.x
	$(BINDIR)/depend.x $(LIB_INCLUDES) $(INCPATH) > dependencies


-include dependencies

