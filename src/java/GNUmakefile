#
# $Id: GNUmakefile,v 1.7 1999-11-16 04:47:46 d3j191 Exp $
#

OBJ_CLASS = nwchem.class nwchem_Main.class nwchem_Filter.class \
	nwchem_Task.class nwchem_Job.class nwchem_MD.class \
	nwchem_MD_energy.class nwchem_MD_optimize.class \
	nwchem_MD_dynamics.class nwchem_MD_thermo.class \
	nwchem_Synch.class nwchem_Input.class nwchem_Times.class\
	nwchem_MD_control.class nwchem_Property.class nwchem_Rasmol.class\
        nwchem_NWChem.class \
        nwchem_Prepare_topology.class nwchem_Prepare_restart.class\
        nwchem_Prepare.class \
	IntegerField.class DoubleField.class Graph.class

CLASS_PATH = .:$(PTOLEMY):$(SWING)

.SUFFIXES:	.f .o .java .class

.java.class:	
		(export CLASSPATH ; CLASSPATH=${CLASS_PATH}; javac $<)

all:	
	$(MAKE) javacheck
	$(MAKE) classes

CLASSDIR = $(NWCHEM_TOP)/lib/classes

classes:
	-test -d $(CLASSDIR) || $(MKDIR) -p $(CLASSDIR) 
	$(MAKE) $(OBJ_CLASS)
	-mv *.class $(CLASSDIR)

javacheck:	whichjava.class
	java whichjava || exit 911
