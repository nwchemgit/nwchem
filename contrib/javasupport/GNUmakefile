# $Id$


include ../../src/config/makefile.h
include ../../src/config/makelib.h

.SUFFIXES:	.f .o .java .class

.java.class:	
		javac $<

all:	
	$(MAKE) javacheck.class
	$(MAKE) classes

GIF_FILES = pc22.gif
CLASSDIR = $(TOPDIR)/lib/classes

classes:
	-test -d $(CLASSDIR) || $(MKDIR) -p $(CLASSDIR) 
	$(MAKE) NWChem.class
	-mv *.class  $(CLASSDIR)
	-cp $(GIF_FILES) $(CLASSDIR)

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
