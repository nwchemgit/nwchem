# $Id: GNUmakefile,v 1.2 1998-05-29 21:27:09 d3e129 Exp $


include ../../src/config/makefile.h
include ../../src/config/makelib.h

.SUFFIXES:	.f .o .java .class

.java.class:	
		javac $<

all:	
	$(MAKE) javacheck.class
	$(MAKE) classes

CLASSDIR = $(TOPDIR)/lib/classes

classes:
	-test -d $(CLASSDIR) || $(MKDIR) -p $(CLASSDIR) 
	$(MAKE) NWChem.class
	-mv *.class $(CLASSDIR)

NWChem.class: 	NWChem.java Email.java submitted.java notsent.java waiting.java
	javac NWChem.java

Email.class:	Email.java
	javac Email.java

submitted.class:	submitted.java
	javac submitted.java

waiting.class:	waiting.java
	javac waiting.java

notsent.class:	notsent.java
	javac notsent.java

javacheck.class:
realclean clean:
	rm -f *.class
