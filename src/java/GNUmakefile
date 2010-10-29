#
# $Id$
#

OBJ_CLASS = PrintfFormat.class PrintfWriter.class PrintfStream.class\
 	nwchem.class nwchem_Main.class nwchem_Filter.class \
	nwchem_Task.class nwchem_Job.class nwchem_MD.class \
	nwchem_MD_energy.class nwchem_MD_optimize.class \
	nwchem_MD_dynamics.class nwchem_MD_thermo.class \
        nwchem_Free.class nwchem_Analysis.class nwchem_RMS.class\
	nwchem_Synch.class nwchem_Input.class nwchem_Times.class nwchem_Timing.class\
	nwchem_MD_control.class nwchem_Property.class nwchem_Rasmol.class\
        nwchem_NWChem.class nwchem_Socket.class\
        nwchem_Prepare_topology.class nwchem_Prepare_restart.class\
        nwchem_Prepare.class nwchem_Rama.class nwchem_Param.class nwchem_Segment.class nwchem_Fragment.class\
	IntegerField.class DoubleField.class Graph.class \
	Atom.class Cros.class Bond.class Angle.class Torsion.class Improper.class Rule.class Segment.class Fragment.class TimingSet.class\
	AtomDefinition.class BondDefinition.class AngleDefinition.class TorsionDefinition.class ImproperDefinition.class

CLASS_PATH = .:${MYCLASSES}:$(PTPLOT_APPLICATION)

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
	-cd $(CLASSDIR)
	-chmod a+rx *

javacheck:	whichjava.class
	java whichjava || exit 911
