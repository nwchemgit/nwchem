
# $Id: GNUmakefile,v 1.47 1997-01-17 01:45:13 d3e129 Exp $

include config/makefile.h

  SUBDIRS =  $(NWSUBDIRS)

  LIB_DEFINES =  -DCOMPILATION_DATE="'`date`'" \
                 -DCOMPILATION_DIR="'$(NWCHEM_TOP)'"

#
# This should be the first target so that a simple make builds
# the entire world
#

.PHONY:	nwchem_config

.PRECIOUS:	$(BINDIR)/nwchem nwchem
all $(BINDIR)/nwchem nwchem:	directories libraries
	@/bin/rm -f nwchem.o stubs.o
	$(MAKE) nwchem.o stubs.o
ifeq ($(TARGET),SOLARIS)
	/bin/rm -f $(BINDIR)/nwchem
endif
	$(LINK.f) -o $(BINDIR)/nwchem nwchem.o stubs.o -ltask $(LIBS)
	/bin/rm -f nwchem.o stubs.o
	$(INSTALL)
#
# This to link only and avoid tedious subdir makes
#

link nwchem_link:    
	@/bin/rm -f nwchem.o stubs.o
	$(MAKE) nwchem.o stubs.o
ifeq ($(TARGET),SOLARIS)
	/bin/rm -f $(BINDIR)/nwchem
endif
	$(LINK.f) -o $(BINDIR)/nwchem nwchem.o stubs.o -ltask $(LIBS)
	/bin/rm -f nwchem.o stubs.o
	$(INSTALL)
#
# This for linking version with profiling
#
prof profile:    
	@/bin/rm -f nwchem.o stubs.o
	$(MAKE) nwchem.o stubs.o
ifeq ($(TARGET),SOLARIS)
	/bin/rm -f nwchem_prof
endif
	$(LINK.f) -p -o nwchem_prof nwchem.o stubs.o -ltask $(LIBS)
	/bin/rm -f nwchem.o stubs.o
	$(INSTALL)
#
# This for linking test version
#
test:	nwchem.o stubs.o
ifeq ($(TARGET),SOLARIS)
	/bin/rm -f nwchem_test
endif
	$(LINK.f) -o nwchem_test nwchem.o stubs.o -ltask $(LIBS)
# removed because if testing don't recompile these:	/bin/rm -f nwchem.o stubs.o
	$(INSTALL)

# This dependency so that includes which are made by libraries rule
# are updated before nwchem is compiled

libraries:
	@for dir in $(SUBDIRS); do \
		echo Making $@ in $$dir; \
		$(MAKE)	-C $$dir || exit 1 ;  \
	done

showdir: 
	@echo "subdirectories in this build <$(SUBDIRS)>"

clean:
	-$(RM) nwchem core *.o *.a *.trace
	@for dir in $(SUBDIRS) config ../tests; do \
		echo Making $@ in $$dir; \
		$(MAKE)	-C $$dir clean || exit 1 ;  \
	done

realclean:	
	-$(RM) nwchem core *.o *.a *.trace makefile.bak
	-$(RM) *~ \#*\#
	@for dir in $(SUBDIRS) config ../tests; do \
		echo Making $@ in $$dir; \
		$(MAKE)	-C $$dir realclean || exit 1 ;  \
	done
	-$(RM) $(LIBDIR)/*.a $(BINDIR)/*
	-$(RM) *.db movecs fort.41 fort.44

.PHONY:	cleanF
cleanF:
	@for dir in $(SUBDIRS) ../tests; do \
		echo Making $@ in $$dir; \
		$(MAKE)	-C $$dir cleanF || exit 1 ;  \
	done
	@for file in *F; do \
		body=`basename $$file .F` ; \
		if [ -f $$body.f ] ; then \
		  echo $$file and $$body.f both exist ... deleting $$body.f; \
		  /bin/rm -f $$body.f ; \
		fi ; \
        done

sngl_to_dbl dbl_to_sngl:
	for dir in $(SUBDIRS); do \
		$(MAKE)	-C $$dir $@ || exit 1 ;  \
	done
        
#
# This rule to make necessary subdirectories in lib and bin (normally)
#

directories:
	test -d $(LIBDIR) || $(MKDIR) -p $(LIBDIR)
	test -d $(BINDIR) || $(MKDIR) -p $(BINDIR)

# This target will change config/nwchem_config.h include file which
# controls what we're going to build.
nwchem_config:
ifeq ($(NWCHEM_MODULES),)
	@echo
	@echo -----
	@echo You must set NWCHEM_MODULES to the list of modules to be \
		included in NWChem builds.
	@$(MAKE) -C config -f make_nwchem_config list-options
else
	$(MAKE) -C config -f make_nwchem_config NWCHEM_MODULES="$(NWCHEM_MODULES)"
endif

# This is a convenience target that will make sure all subdirectories
# listed in $(SUBDIRS) exist. Useful because CVS will only update
# directories that already exist in your working tree, so this is a
# way ot make sure all that are relevant to the code exist.

make_source_dirs:
ifneq ($(SUBDIRS),)
	for dir in $(SUBDIRS); do \
	   test -d $$dir  || $(MKDIR) $$dir; \
	done
endif

# This is a convenience target to remove source directories that are
# configured out of the program.  Useful to save disk space.

remove_unused_dirs:
ifneq ($(EXCLUDED_SUBDIRS),)
	for dir in $(EXCLUDED_SUBDIRS); do \
	   $(RM) -r $$dir; \
	done
endif

.PHONY:	depend
depend:	
	( cd config; $(MAKE) depend.x ; )
	config/stripdepend
	for dir in $(SUBDIRS); do \
		$(MAKE)	-C $$dir $@ || exit 1 ;  \
	done
	config/depend.x

.PHONY:	src.F
src.F:	
	cat stubs.f nwchem.F > src.F
ifdef SUBDIRS
	for dir in $(SUBDIRS); do \
		$(MAKE)	-C $$dir $@ || exit 1 ;  \
	done
endif


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(nwchem.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh 
