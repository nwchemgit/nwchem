#
# $Id: GNUmakefile,v 1.2 1999-08-10 23:05:15 d3j191 Exp $
#

OBJ_CLASS = nwmd.class GUI_Nwchem.class RUN_Nwchem.class Graph.class

include ../config/makefile.h
include ../config/makelib.h

CLASS_PATH = .:./ptplot.zip:$(TOPDIR)/lib/classes/ptplot.zip:$(SWING)

.SUFFIXES:	.f .o .java .class

.java.class:	
		(export CLASSPATH ; CLASSPATH=${CLASS_PATH}; javac $<)

all:	
	$(MAKE) javacheck
	$(MAKE) classes

CLASSDIR = $(TOPDIR)/lib/classes

classes:
	-test -d $(CLASSDIR) || $(MKDIR) -p $(CLASSDIR) 
	$(MAKE) $(OBJ_CLASS)
	-mv *.class $(CLASSDIR)

javacheck:	whichjava.class
	java whichjava || exit 911
